/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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

#include "noSolidificationModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidificationModels
{
    defineTypeNameAndDebug(none, 0);

    addToRunTimeSelectionTable
    (
        solidificationModel,
        none,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidificationModels::none::none
(
    const fvMesh& mesh,
    const word& group
)
:
    solidificationModel(mesh, group)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidificationModels::none::correct(const volScalarField& T)
{}


Foam::tmp<Foam::fvScalarMatrix> Foam::solidificationModels::none::hSource
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const surfaceScalarField& alphaRhoPhi,
    const volScalarField& T
) const
{
    return tmp<fvScalarMatrix>(nullptr);
}


Foam::tmp<Foam::fvVectorMatrix> Foam::solidificationModels::none::USource
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U
) const
{
    return tmp<fvVectorMatrix>(nullptr);
}

// ************************************************************************* //
