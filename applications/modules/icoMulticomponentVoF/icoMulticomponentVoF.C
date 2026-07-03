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

#include "icoMulticomponentVoF.H"
#include "geometricZeroField.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(icoMulticomponentVoF, 0);
    addToRunTimeSelectionTable(solver, icoMulticomponentVoF, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::icoMulticomponentVoF::icoMulticomponentVoF
(
    fvMesh& mesh
)
:
    multiphaseVoFSolver
    (
        mesh,
        autoPtr<multiphaseVoFMixture>
        (
            new multicomponentVoFMixture(mesh)
        )
    ),

    mixture
    (
        refCast<multicomponentVoFMixture>
        (
            multiphaseVoFSolver::mixture
        )
    ),

    phases(mixture.phases()),

    p(mixture.p()),

    rhoCp(mixture.rhoCp()),

    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict(),
        false
    ),

    momentumTransport_
    (
        compressible::momentumTransportModel::New
        (
            rho,
            U,
            rhoPhi,
            mixture
        )
    ),

    momentumTransport(momentumTransport_()),

    contErr
    (
        new volScalarField
        (
            IOobject
            (
                "contErr",
                mesh.time().name(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimDensity/dimTime, 0.0),
            zeroGradientFvPatchField<scalar>::typeName
        )
    )
{
    if (correctPhi || mesh.topoChanging())
    {
        rAU = new volScalarField
        (
            IOobject
            (
                "rAU",
                runTime.name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime/dimDensity, 1)
        );
    }

    if (!incompressible())
    {
        FatalErrorInFunction
            << "At least one phase is compressible!"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::icoMulticomponentVoF::~icoMulticomponentVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::icoMulticomponentVoF::prePredictor()
{
    multiphaseVoFSolver::prePredictor();

    contErr.ref() = fvc::ddt(rho)() + fvc::div(rhoPhi)();

    forAll(mixture.phases(), phasei)
    {
        const volScalarField& rho = phases[phasei].thermo().rho();
        contErr.ref() -= (fvModels().source(phases[phasei], rho)&rho)();
    }

    contErr.ref().correctBoundaryConditions();
}


void Foam::solvers::icoMulticomponentVoF::momentumTransportPredictor()
{
    momentumTransport.predict();
}


void Foam::solvers::icoMulticomponentVoF::
thermophysicalTransportPredictor()
{}


void Foam::solvers::icoMulticomponentVoF::momentumTransportCorrector()
{
    momentumTransport.correct();
}


void Foam::solvers::icoMulticomponentVoF::
thermophysicalTransportCorrector()
{}


// ************************************************************************* //
