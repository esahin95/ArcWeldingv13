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

#include "multiphaseChangeVoF.H"
#include "geometricZeroField.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(multiphaseChangeVoF, 0);
    addToRunTimeSelectionTable(solver, multiphaseChangeVoF, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::multiphaseChangeVoF::multiphaseChangeVoF
(
    fvMesh& mesh
)
:
    multiphaseVoFSolver
    (
        mesh,
        autoPtr<multiphaseVoFMixture>
        (
            new multiphaseChangeVoFMixture(mesh)
        )
    ),

    mixture
    (
        refCast<multiphaseChangeVoFMixture>
        (
            multiphaseVoFSolver::mixture
        )
    ),

    phases(mixture.phases()),

    p(mixture.p()),

    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict(),
        false
    ),

    K("K", 0.5*magSqr(U)),

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

    momentumTransport(momentumTransport_())
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::multiphaseChangeVoF::~multiphaseChangeVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::multiphaseChangeVoF::prePredictor()
{
    multiphaseVoFSolver::prePredictor();

    contErr = fvc::ddt(rho)()() + fvc::div(rhoPhi)()();

    forAll(mixture.phases(), phasei)
    {
        const volScalarField& rho = phases[phasei].thermo().rho();
        contErr.ref() -= (fvModels().source(phases[phasei], rho)&rho)();
    }
}


void Foam::solvers::multiphaseChangeVoF::momentumTransportPredictor()
{
    momentumTransport.predict();
}


void Foam::solvers::multiphaseChangeVoF::
thermophysicalTransportPredictor()
{}


void Foam::solvers::multiphaseChangeVoF::momentumTransportCorrector()
{
    momentumTransport.correct();
}


void Foam::solvers::multiphaseChangeVoF::
thermophysicalTransportCorrector()
{}


// ************************************************************************* //
