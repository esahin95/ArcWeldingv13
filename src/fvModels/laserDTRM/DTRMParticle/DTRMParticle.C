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
    const scalar absorption
)
:
    particle(searchEngine, position, cellI, nLocateBoundaryHits),
    q0_(power),
    trackIndex_(nParticles++),
    a_(absorption),
    q_(power),
    d_(direction)
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

    const scalar trackTime = td.mesh.time().deltaTValue();

    // Initial position
    td.append
    (
        this->position(td.mesh),
        trackIndex_,
        q_
    );

    const scalar limitAlpha = 0.5 + 0.001;

    while (td.keepParticle && td.sendToProc == -1 && stepFraction() < 1)
    {
        if (debug)
        {
            Info<< "Time = " << td.mesh.time().name()
                << " trackTime = " << trackTime
                << " stepFraction() = " << stepFraction()
                << " trackIndex = "<< trackIndex_ << endl;
        }

        scalar alpha = td.alphaInterp().interpolate
            (
                this->coordinates(),
                this->currentTetIndices(td.mesh)
            );

        DebugInfo<<"Tracking alpha: "<< alpha;
        while
        (
            alpha > limitAlpha &&
            td.keepParticle && td.sendToProc == -1 && stepFraction() < 1
        )
        {
            trackToAndHitFace(d_, 1.0, cloud, td);

            alpha = td.alphaInterp().interpolate
                (
                    this->coordinates(),
                    this->currentTetIndices(td.mesh)
                );

            DebugInfo<< " " << alpha;
        }
        DebugInfo<<endl;

        if (alpha <= limitAlpha)
        {
            // Final position
            td.append
            (
                this->position(td.mesh),
                trackIndex_,
                q_
            );

            // Reflection
            td.Q(this->cell()) += q_;
            q_ = Zero;
            stepFraction() = 1.0;
        }
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