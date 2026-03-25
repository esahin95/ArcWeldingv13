/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2025 OpenFOAM Foundation
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

#include "compressibleReactingVoF.H"
#include "localEulerDdtScheme.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(compressibleReactingVoF, 0);
    addToRunTimeSelectionTable(solver, compressibleReactingVoF, fvMesh);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::compressibleReactingVoF::compressibleReactingVoF(fvMesh& mesh)
:
    compressibleVoF(mesh),
    
    dict_
    (
        IOobject
        (
            "phaseChangeProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    phase_(dict_.lookup<word>("phase")),

    alphaSolidT_
    (
        Function1<scalar>::New
        (
            "alphaSolidT",
            dimTemperature,
            unitFraction,
            dict_
        )
    ),

    alphaSolid_
    (
        IOobject
        (
            IOobject::groupName(phase_, "solid"),
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0.0),
        zeroGradientFvPatchField<scalar>::typeName
    ),

    alphaVoF_
    (
        mesh_.lookupObjectRef<volScalarField>
        (
            IOobject::groupName("alpha", phase_)
        )
    ),

    L_
    (
        "L",
        dimEnergy/dimMass,
        dict_.lookup<scalar>("L")
    ),

    relax_(dict_.lookupOrDefault<scalar>("relax", 1.0)),

    Cu_
    (
        "Cu",
        dimDensity/dimTime,
        dict_.lookupOrDefault<scalar>("Cu", 1e5)
    ), 

    q_(dict_.lookupOrDefault<scalar>("q", 0.001))
{

    // Initialize solid fraction field
    const volScalarField& T = mixture_.T();
    alphaSolid_.primitiveFieldRef() = 
        min 
        (
            alphaVoF_,
            alphaSolidT_->value(T.primitiveField())
        );
    alphaSolid_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::compressibleReactingVoF::~compressibleReactingVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleReactingVoF::prePredictor()
{
    compressibleVoF::prePredictor();
}

/*
void Foam::solvers::compressibleReactingVoF::momentumTransportPredictor()
{
    momentumTransport.predict();
}


void Foam::solvers::compressibleReactingVoF::thermophysicalTransportPredictor()
{
    thermophysicalTransport.predict();
}


void Foam::solvers::compressibleReactingVoF::momentumTransportCorrector()
{
    momentumTransport.correct();
}


void Foam::solvers::compressibleReactingVoF::thermophysicalTransportCorrector()
{
    thermophysicalTransport.correct();
}
*/

// ************************************************************************* //
