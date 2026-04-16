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
    const scalar absorption,
    const label trackIndex
)
:
    particle(searchEngine, position, cellI, nLocateBoundaryHits),
    q0_(power),
    trackIndex_(trackIndex),
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
            alpha > 0.5 &&
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

        if (alpha <= 0.5)
        {
            td.Q(this->cell()) += q_;
            q_ = Zero;
            stepFraction() = 1.0;
        }

        /*
        const scalar sfrac = stepFraction();

        const scalar f = 1 - stepFraction();
        trackToAndHitFace(f*trackTime*U_, f, cloud, td);

        const scalar dt = (stepFraction() - sfrac)*trackTime;

        const tetIndices tetIs = this->currentTetIndices(td.mesh);
        scalar rhoc = td.rhoInterp().interpolate(this->coordinates(), tetIs);
        vector Uc = td.UInterp().interpolate(this->coordinates(), tetIs);
        scalar nuc = td.nuInterp().interpolate(this->coordinates(), tetIs);

        scalar rhop = cloud.rhop();
        scalar magUr = mag(Uc - U_);

        scalar ReFunc = 1.0;
        scalar Re = magUr*d_/nuc;

        if (Re > 0.01)
        {
            ReFunc += 0.15*pow(Re, 0.687);
        }

        scalar Dc = (24.0*nuc/d_)*ReFunc*(3.0/4.0)*(rhoc/(d_*rhop));

        U_ = (U_ + dt*(Dc*Uc + (1.0 - rhoc/rhop)*td.g()))/(1.0 + dt*Dc);
        */
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