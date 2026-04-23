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
    trackIndex_(nParticles++),
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

    //const label origProc = this->origProc();

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
        const scalar oldAlpha = td.alphaInterp().interpolate
            (
                this->coordinates(),
                this->currentTetIndices(td.mesh)
            );
        const scalar oldAbsorp = td.absorpInterp().interpolate
            (
                this->coordinates(),
                this->currentTetIndices(td.mesh)
            );
        const label oldCell = this->cell();
        const vector oldPos = this->position(td.mesh);

        // Track to new face and cell
        trackToFace(td.mesh, d_, 1.0);

        // New data
        const scalar alpha = td.alphaInterp().interpolate
            (
                this->coordinates(),
                this->currentTetIndices(td.mesh)
            );
        const vector pos = this->position(td.mesh);

        // Distance traveled by ray
        const scalar ds = mag(pos - oldPos);

        // Handle reflection
        const scalar reflectivity = 0.5;
        const scalar TOL = 1e-3;
        if (oldAlpha >= 0.5 && alpha <= 0.5 && transmissive_)
        {
            DTRMParticle* pPtr(new DTRMParticle(*this));

            // Bounds
            vector lBound = oldPos;
            vector rBound = pos;
            scalar al = oldAlpha;
            scalar ar = alpha;
            DebugInfo<< "Reflection detected: " << al << " - " << ar;

            // Initial interface guess
            scalar t = (0.5 - ar) / (al - ar + TOL);
            vector p = t * lBound + (1-t) * rBound;

            const meshSearch& searchEngine = meshSearch::New(td.mesh);
            pPtr->locate(searchEngine, p, oldCell);
            scalar a = td.alphaInterp().interpolate
                (
                    pPtr->coordinates(),
                    pPtr->currentTetIndices(td.mesh)
                );
            DebugInfo<< " : " << a;

            // Track interface position
            label i = 0;
            while (mag(a - 0.5) > TOL)
            {
                if (a > 0.5)
                {
                    lBound = p;
                    al = a;
                }
                else
                {
                    rBound = p;
                    ar = a;
                }

                t = (0.5 - ar) / (al - ar + TOL);
                p = t * lBound + (1-t) * rBound;

                pPtr->locate(searchEngine, p, oldCell);
                a = td.alphaInterp().interpolate
                (
                    pPtr->coordinates(),
                    pPtr->currentTetIndices(td.mesh)
                );

                i++;
                if (i > 20)
                {
                    break;
                }
            }
            DebugInfo<< " -> " << a << " i = " << i << endl;

            vector nHat = td.nHatInterp().interpolate
                (
                    pPtr->coordinates(),
                    pPtr->currentTetIndices(td.mesh)
                );
            pPtr->d_ = normalised(d_ - 2.0 * (nHat & d_) * nHat);
            pPtr->q_ = reflectivity * q_;

            //if (origProc == Pstream::myProcNo())
            //{
                cloud.addParticle(pPtr);
            //}
            //DebugInfo<<pPtr<<endl;

            q_ -= pPtr->q_;
            transmissive_ = false;
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

        // Patch interactions
        hitFace(d_, 1.0, cloud, td);
    }

    return td.keepParticle;
}


void Foam::DTRMParticle::hitProcessorPatch
(
    lagrangian::Cloud<DTRMParticle>& cloud,
    trackingData& td
)
{
    particle::hitProcessorPatch(cloud, td);

    td.append
    (
        this->position(td.mesh),
        trackIndex_,
        q_
    );
}

void Foam::DTRMParticle::hitWallPatch
(
    lagrangian::Cloud<DTRMParticle>& cloud,
    trackingData& td
)
{
    td.keepParticle = false;

    td.append
    (
        this->position(td.mesh),
        trackIndex_,
        q_
    );
}


// ************************************************************************* //