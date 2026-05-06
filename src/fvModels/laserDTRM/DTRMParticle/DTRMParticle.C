/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "DTRMParticle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(DTRMParticle, 0);

    label DTRMParticle::nParticles = 0;

    scalar DTRMParticle::qLost = 0.0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DTRMParticle::DTRMParticle
(
    const meshSearch& searchEngine,
    const vector& position,
    const label cellI,
    label& nLocateBoundaryHits,
    const vector& direction,
    const scalar power,
    const bool transmissive
)
:
    particle(searchEngine, position, cellI, nLocateBoundaryHits),
    q0_(power),
    trackIndex_(nParticles++),
    q_(power),
    d_(direction),
    transmissive_(transmissive)
{
    this->reset(0.0);
}

Foam::DTRMParticle::DTRMParticle(const DTRMParticle& p)
:
    particle(p),
    q0_(p.q0_),
    trackIndex_(p.trackIndex_),
    q_(p.q_),
    d_(p.d_),
    transmissive_(p.transmissive_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::DTRMParticle::move
(
    lagrangian::Cloud<DTRMParticle>& cloud,
    trackingData& td
)
{
    td.keepParticle = true;
    td.sendToProc = -1;

    // Initial position
    td.append
    (
        this->position(td.mesh),
        trackIndex_,
        q_
    );

    while
    (
        q_ > 0.01 * q0_ &&
        td.keepParticle && td.sendToProc == -1 && stepFraction() < 1
    )
    {
        // Initial data
        const tetIndices oldTetIs = this->currentTetIndices(td.mesh);
        const scalar oldAlpha = td.alphaInterp().interpolate
            (
                this->coordinates(),
                oldTetIs
            );
        const scalar oldAbsorp = td.absorpInterp().interpolate
            (
                this->coordinates(),
                oldTetIs
            );
        const label oldCell = this->cell();
        const vector oldPos = this->position(td.mesh);

        const scalar TOL = 1e-3;
        if (transmissive_ && oldAlpha < 0.5 - TOL)
        {
            DebugInfo<< "Transmissive ray inside material detected"<<endl;
            qLost += q_;
            td.keepParticle = false;
            continue;
        }

        // Track to new face and cell
        trackToAndHitFace(d_, 1.0, cloud, td);

        // New data
        const tetIndices tetIs = this->currentTetIndices(td.mesh);
        const scalar alpha = td.alphaInterp().interpolate
            (
                this->coordinates(),
                tetIs
            );
        const vector pos = this->position(td.mesh);

        // Distance traveled by ray
        const vector dx = pos - oldPos;
        const scalar ds = mag(dx);

        // Check for reflection
        if
        (
            oldAlpha >= 0.5 && alpha < 0.5 && transmissive_
        )
        {
            const meshSearch& searchEngine = td.searchEngine();
            label nLocateBoundaryHits = 0;

            DTRMParticle* pPtr = new
                DTRMParticle
                (
                    searchEngine,
                    oldPos,
                    oldCell,
                    nLocateBoundaryHits,
                    d_,
                    q0_,
                    true
                );

            // Initial Bounds
            scalar a = 0, fa = oldAlpha - 0.5;
            scalar b = 1, fb = alpha - 0.5;
            scalar c = 0, fc = fa;

            // Track to interface
            label i = 0, MAX = 20;
            while (mag(fc) > TOL && i < MAX)
            {
                i++;

                c = (a*fb - b*fa) / (fb - fa);
                const vector p = oldPos + c * dx;

                pPtr->locate(searchEngine, p, oldCell);
                fc = td.alphaInterp().interpolate
                    (
                        pPtr->coordinates(),
                        pPtr->currentTetIndices(td.mesh)
                    ) - 0.5;

                if (fc * fa > 0)
                {
                    a = c;
                    fa = fc;
                }
                else
                {
                    b = c;
                    fb = fc;
                }
            }
            if (i >= MAX)
            {
                DebugInfo<<"Failed search: "
                         << oldAlpha << " -> " << fc+0.5
                         << endl;
            }

            const vector gradAlpha = td.gradAlphaInterp().interpolate
                (
                    pPtr->coordinates(),
                    pPtr->currentTetIndices(td.mesh)
                );
            const vector nHat = gradAlpha / mag(gradAlpha);
            scalar cosTheta = -nHat & d_;
            if (cosTheta < 0.0)
            {
                DebugInfo<< "False reflection direction" << endl;
                qLost += q_;
                td.keepParticle = false;
                delete pPtr;
            }
            else
            {
                pPtr->d_ = td.reflection().R(d_, nHat);
                pPtr->q_ = td.reflection().rho(cosTheta) * q_;
                //DebugInfo<< cosTheta << " " << (pPtr->q_ / q_)<<endl;

                q_ -= pPtr->q_;
                transmissive_ = false;

                if (pPtr->cell() != -1)
                {
                    cloud.addParticle(pPtr);
                }
            }
        }

        // Laser power absorption in ray
        if (!transmissive_)
        {
            const scalar qAbsorped = max(min(ds * oldAbsorp, 1.0), 0.0) * q_;
            td.Q(oldCell) += qAbsorped;
            q_ -= qAbsorped;
        }

        // New position
        td.append
        (
            pos,
            trackIndex_,
            q_
        );
    }

    return td.keepParticle;
}

void Foam::DTRMParticle::hitWallPatch
(
    lagrangian::Cloud<DTRMParticle>& cloud,
    trackingData& td
)
{
    td.keepParticle = false;
}


// ************************************************************************* //