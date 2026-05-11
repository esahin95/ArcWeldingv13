/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

#include "phaseChangeVoFphase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeVoFphase::phaseChangeVoFphase
(
    const word& name,
    const fvMesh& mesh,
    const volScalarField& T
)
:
    compressibleVoFphase(name, mesh, T),
    solidFraction_
    (
        IOobject
        (
            IOobject::groupName("solid", name),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
Foam::autoPtr<Foam::VoFphase> Foam::phaseChangeVoFphase::clone() const
{
    NotImplemented;
    return autoPtr<VoFphase>(nullptr);
}


void Foam::phaseChangeVoFphase::correct
(
    const volScalarField& p,
    const volScalarField& T
)
{
    thermo_->he() = thermo_->he(p, T);
    thermo_->correct();
}
*/

// ************************************************************************* //
