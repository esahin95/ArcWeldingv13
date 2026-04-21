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
        // Initial values
        const scalar alphaOld = td.alphaInterp().interpolate
            (
                this->coordinates(),
                this->currentTetIndices(td.mesh)
            );
        const scalar absorpOld = td.absorpInterp().interpolate
            (
                this->coordinates(),
                this->currentTetIndices(td.mesh)
            );
        const label origProcOld = this->origProc();
        const label cellIOld = this->cell();
        const vector posOld = this->position(td.mesh);

        // Track to new face and cell
        trackToAndHitFace(d_, 1.0, cloud, td);

        // New values
        const scalar alpha = td.alphaInterp().interpolate
            (
                this->coordinates(),
                this->currentTetIndices(td.mesh)
            );
        const vector pos = this->position(td.mesh);

        const scalar ds = mag(pos - posOld);

        // Handle reflection
        const scalar reflectivity = 0.5;
        if (alphaOld > 0.5 && alpha <= 0.5)
        {
            // Bounds
            vector lBound = posOld;
            vector rBound = pos;
            scalar al = alphaOld;
            scalar ar = alpha;
            DebugInfo<< "Reflection detected: " << al << " - " << ar;

            // Initial interface guess
            scalar t = (0.5 - ar) / (al - ar);
            vector p = t * lBound + (1-t) * rBound;
            scalar a = td.alphaInterp().interpolate(p, cellIOld);
            DebugInfo<< " : " << a;

            // Track interface position
            label i = 0;
            while (mag(a - 0.5) > 0.01)
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

                t = (0.5 - ar) / (al - ar);
                p = t * lBound + (1-t) * rBound;
                a = td.alphaInterp().interpolate(p, cellIOld);

                i++;
                if (i > 20)
                {
                    break;
                }
            }
            DebugInfo<< " -> " << a << " i = " << i << endl;

            // New particle at interface
            vector nHat = td.nHatInterp().interpolate(p, cellIOld);
            nHat /= mag(nHat);
            const scalar qReflected = reflectivity * q_;
            if (origProcOld == Pstream::myProcNo())
            {
                cloud.addParticle(new DTRMParticle
                    (
                        td.searchEngine,
                        p,
                        cellIOld,
                        td.nLocateBoundaryHits,
                        normalised(d_ - 2.0 * (nHat & d_) * nHat),
                        qReflected,
                        true
                    )
                );
            }

            q_ -= qReflected;

            transmissive_ = false;
        }

        // Laser power absorption in ray
        if (!transmissive_)
        {
            const scalar qAbsorped = max(min(ds * absorpOld, 1.0), 0.0) * q_;
            td.Q(cellIOld) += qAbsorped;
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