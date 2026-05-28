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
//#include "geometricZeroField.H"
#include "fvcMeshPhi.H"
#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
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
    twoPhaseVoFSolver
    (
        mesh,
        autoPtr<twoPhaseVoFMixture>(new compressibleTwoPhaseVoFMixture(mesh))
    ),

    mixture_
    (
        refCast<compressibleTwoPhaseVoFMixture>(twoPhaseVoFSolver::mixture)
    ),

    p(mixture_.p()),

    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict(),
        false
    ),

    alphaRhoPhi1
    (
        IOobject::groupName("alphaRhoPhi", alpha1.group()),
        fvc::interpolate(mixture_.thermo1().rho())*alphaPhi1
    ),

    alphaRhoPhi2
    (
        IOobject::groupName("alphaRhoPhi", alpha2.group()),
        fvc::interpolate(mixture_.thermo2().rho())*alphaPhi2
    ),

    rhoCp_
    (
        "rhoCp",
        (
            mixture_.thermo1().rho() * alpha1 * mixture_.thermo1().Cp()
          + mixture_.thermo2().rho() * alpha2 * mixture_.thermo2().Cp()
        )
    ),

    Cp_
    (
        IOobject
        (
            "Cp",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rhoCp_ / rho
    ),

    rhoPhiCp_
    (
        "rhoPhiCp",
        rhoPhi * fvc::interpolate(Cp_)
    ),

    momentumTransport
    (
        rho,
        U,
        phi,
        rhoPhi,
        alphaPhi1,
        alphaPhi2,
        alphaRhoPhi1,
        alphaRhoPhi2,
        mixture_
    ),

    thermophysicalTransport(momentumTransport)
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

    if (!runTime.restart() || !divergent())
    {
        correctUphiBCs(U_, phi_, true);

        fv::correctPhi
        (
            phi_,
            U,
            p_rgh,
            rAU,
            autoPtr<volScalarField>(),
            pressureReference(),
            pimple
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::icoPhaseChangeVoF::~icoPhaseChangeVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::icoPhaseChangeVoF::prePredictor()
{
    twoPhaseVoFSolver::prePredictor();

    const volScalarField& rho1 = mixture_.thermo1().rho();
    const volScalarField& rho2 = mixture_.thermo2().rho();

    const volScalarField& Cp1 = mixture_.thermo1().Cp();
    const volScalarField& Cp2 = mixture_.thermo2().Cp();

    // Mass fluxes
    alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;
    alphaRhoPhi2 = fvc::interpolate(rho2)*alphaPhi2;

    rhoPhi = alphaRhoPhi1 + alphaRhoPhi2;

    // Heat capacity
    rhoCp_ = alpha1*rho1*Cp1 + alpha2*rho2*Cp2;

    Cp_ = rhoCp_ / rho;

    rhoPhiCp_ = alphaRhoPhi1*fvc::interpolate(Cp1)
              + alphaRhoPhi2*fvc::interpolate(Cp2);
}


void Foam::solvers::icoPhaseChangeVoF::momentumTransportPredictor()
{
    momentumTransport.predict();
}


void Foam::solvers::icoPhaseChangeVoF::thermophysicalTransportPredictor()
{
    thermophysicalTransport.predict();
}


void Foam::solvers::icoPhaseChangeVoF::pressureCorrector()
{
    incompressiblePressureCorrector(p);
}


void Foam::solvers::icoPhaseChangeVoF::momentumPredictor()
{
    twoPhaseVoFSolver::momentumPredictor();
}


void Foam::solvers::icoPhaseChangeVoF::thermophysicalPredictor()
{
    //const volScalarField& rho1(mixture_.rho1());
    //const volScalarField& rho2(mixture_.rho2());

    volScalarField& T = mixture_.T();

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp_,T)
        + fvm::div(rhoPhiCp_, T)
        - fvm::Sp(fvc::ddt(rhoCp_) + fvc::div(rhoPhiCp_), T)
        - fvm::laplacian(thermophysicalTransport.kappaEff(), T)
    );

    TEqn.relax();

    fvConstraints().constrain(TEqn);

    solve(TEqn);

    fvConstraints().constrain(T);

    mixture_.correctThermo();
    mixture_.correct();
}


void Foam::solvers::icoPhaseChangeVoF::momentumTransportCorrector()
{
    momentumTransport.correct();
}


void Foam::solvers::icoPhaseChangeVoF::thermophysicalTransportCorrector()
{
    thermophysicalTransport.correct();
}


// ************************************************************************* //
