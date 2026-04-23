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
    const lagrangian::Cloud<DTRMParticle>& cloud,
    const interpolationCellPoint<scalar>& alphaInterp,
    const interpolationCellPoint<scalar>& absorpInterp,
    const interpolationCellPoint<vector>& nHatInterp,
    volScalarField& Q,
    DynamicField<vector>& allPositions,
    DynamicField<label>& allTracks,
    DynamicField<scalar>& allPowers
)
:
    particle::trackingData(cloud),
    alphaInterp_(alphaInterp),
    absorpInterp_(absorpInterp),
    nHatInterp_(nHatInterp),
    Q_(Q),
    allPositions_(allPositions),
    allTracks_(allTracks),
    allPowers_(allPowers)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::DTRMParticle::trackingData::alphaInterp() const
{
    return alphaInterp_;
}

inline const Foam::interpolationCellPoint<Foam::scalar>&
Foam::DTRMParticle::trackingData::absorpInterp() const
{
    return absorpInterp_;
}

inline const Foam::interpolationCellPoint<Foam::vector>&
Foam::DTRMParticle::trackingData::nHatInterp() const
{
    return nHatInterp_;
}

inline Foam::scalar& Foam::DTRMParticle::trackingData::Q(label cellI)
{
    return Q_[cellI];
}

inline void Foam::DTRMParticle::trackingData::append
(
    vector position,
    label trackIndex,
    scalar power
)
{
    if (mesh.time().writeTime())
    {
        allPositions_.append(position);
        allTracks_.append(trackIndex);
        allPowers_.append(power);
    }
}

// ************************************************************************* //
