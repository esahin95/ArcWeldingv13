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

#include "linearExplicitSolidificationModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidificationModels
{
    defineTypeNameAndDebug(linearExplicit, 0);

    addToRunTimeSelectionTable
    (
        solidificationModel,
        linearExplicit,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidificationModels::linearExplicit::linearExplicit
(
    const fvMesh& mesh,
    const word& group
)
:
    solidificationModel(mesh, group),

    Lm_
    (
        "Lm",
        dimEnergy/dimMass,
        dict.lookup<scalar>("Lm")
    ),

    Tliq_
    (
        "Tliq",
        dimTemperature,
        dict.lookup<scalar>("Tliq")
    ),

    Tsol_
    (
        "Tsol",
        dimTemperature,
        dict.lookup<scalar>("Tsol")
    ),

    Cu_
    (
        "Cu",
        dimDensity/dimTime,
        dict.lookupOrDefault<scalar>("Cu", 1e5)
    ),

    q_(dict.lookupOrDefault<scalar>("q", 0.001)),

    relax_(dict.lookupOrDefault<scalar>("relax", 1.0))
{
    Info<<"   Lm   = " << Lm_ << endl;
    Info<<"   Tliq = " << Tliq_ << endl;
    Info<<"   Tsol = " << Tsol_ << endl;
    Info<<"   Cu   = " << Cu_ << endl;
    Info<<"   q    = " << q_ << endl;
    Info<<"   relax= " << relax_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidificationModels::linearExplicit::correct
(
    const volScalarField& T
)
{}


Foam::tmp<Foam::fvScalarMatrix>
Foam::solidificationModels::linearExplicit::hSource
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const surfaceScalarField& alphaRhoPhi,
    const volScalarField& T
) const
{
    return tmp<fvScalarMatrix>(nullptr);
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::solidificationModels::linearExplicit::USource
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volVectorField& U
) const
{
    return tmp<fvVectorMatrix>(nullptr);
}

// ************************************************************************* //
