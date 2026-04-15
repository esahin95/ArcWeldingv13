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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::DTRMParticle::trackingData::trackingData
(
    const DTRMParticleCloud& spc,
    const interpolationCellPoint<scalar>& rhoInterp,
    const interpolationCellPoint<vector>& UInterp,
    const interpolationCellPoint<scalar>& nuInterp,
    const vector& g
)
:
    particle::trackingData(spc),
    rhoInterp_(rhoInterp),
    UInterp_(UInterp),
    nuInterp_(nuInterp),
    g_(g)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::DTRMParticle::trackingData::rhoInterp() const
{
    return rhoInterp_;
}


inline const Foam::interpolationCellPoint<Foam::vector>&
Foam::DTRMParticle::trackingData::UInterp() const
{
    return UInterp_;
}


inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::DTRMParticle::trackingData::nuInterp() const
{
    return nuInterp_;
}

inline const Foam::vector& Foam::DTRMParticle::trackingData::g() const
{
    return g_;
}


inline Foam::scalar Foam::DTRMParticle::d() const
{
    return d_;
}


inline const Foam::vector& Foam::DTRMParticle::U() const
{
    return U_;
}


// ************************************************************************* //
