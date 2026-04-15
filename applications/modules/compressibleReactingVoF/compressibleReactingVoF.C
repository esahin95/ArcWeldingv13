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

    alphaVoF_
    (
        mesh_.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", phase_)
        )
    ),

    solidFraction_
    (
        IOobject
        (
            IOobject::groupName(alphaVoF_.name(), "solidFraction"),
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0.0)
    ),

    alphaSolid_
    (
        IOobject
        (
            IOobject::groupName(alphaVoF_.name(), "solid"),
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0.0)
    ),

    L_
    (
        "L",
        dimEnergy/dimMass,
        dict_.lookup<scalar>("L")
    ),

    Tsol_
    (
        "Tsol",
        dimTemperature,
        dict_.lookup<scalar>("Tsol")
    ),

    Tliq_
    (
        "Tliq",
        dimTemperature,
        dict_.lookup<scalar>("Tliq")
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

    // Initialize solid phase fraction field
    const volScalarField& T = mixture_.T();

    solidFraction_ =
        max
        (
            min
            (
                (Tliq_ - T) / (Tliq_ - Tsol_),
                1.0
            ),
            0.0
        );
    alphaSolid_ = alphaVoF_ * solidFraction_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::compressibleReactingVoF::~compressibleReactingVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::compressibleReactingVoF::prePredictor()
{
    compressibleVoF::prePredictor();
}

// ************************************************************************* //
