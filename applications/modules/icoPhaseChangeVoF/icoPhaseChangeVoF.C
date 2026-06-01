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

#include "icoPhaseChangeVoF.H"
#include "localEulerDdtScheme.H"
#include "fvCorrectPhi.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(icoPhaseChangeVoF, 0);
    addToRunTimeSelectionTable(solver, icoPhaseChangeVoF, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::icoPhaseChangeVoF::icoPhaseChangeVoF(fvMesh& mesh)
:
    icoThermoVoF(mesh),

    phaseChangeDict_
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

    solidFraction_
    (
        IOobject
        (
            IOobject::groupName(alpha1.name(), "solidFraction"),
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
            IOobject::groupName(alpha1.name(), "solid"),
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0.0)
    ),

    Lm_
    (
        "Lm",
        dimEnergy/dimMass,
        phaseChangeDict_.lookup<scalar>("Lm")
    ),

    Tsol_
    (
        "Tsol",
        dimTemperature,
        phaseChangeDict_.lookup<scalar>("Tsol")
    ),

    Tliq_
    (
        "Tliq",
        dimTemperature,
        phaseChangeDict_.lookup<scalar>("Tliq")
    ),

    relax_(phaseChangeDict_.lookupOrDefault<scalar>("relax", 1.0)),

    Cu_
    (
        "Cu",
        dimDensity/dimTime,
        phaseChangeDict_.lookupOrDefault<scalar>("Cu", 1e5)
    ),

    q_(phaseChangeDict_.lookupOrDefault<scalar>("q", 0.001))
{
    const volScalarField& T = mixture_.T();

    // Initialize solid phase fraction field
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
    alphaSolid_ = alpha1 * solidFraction_;
    alphaSolid_.write();

    Info<<Tsol_ << " " <<Tliq_ << " " <<Lm_ << " " <<q_ << " " <<Cu_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::icoPhaseChangeVoF::~icoPhaseChangeVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //




// ************************************************************************* //
