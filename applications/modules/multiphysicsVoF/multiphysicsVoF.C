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

#include "multiphysicsVoF.H"
#include "localEulerDdtScheme.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcAverage.H"
#include "interfaceCompression.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(multiphysicsVoF, 0);
    addToRunTimeSelectionTable(solver, multiphysicsVoF, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::multiphysicsVoF::correctCoNum()
{
    twoPhaseSolver::correctCoNum();

    const scalarField sumPhi
    (
        interface.nearInterface()().primitiveField()
       *fvc::surfaceSum(mag(phi))().primitiveField()
    );

    alphaCoNum =
        0.5*gMax(sumPhi/mesh.V().primitiveField())*runTime.deltaTValue();

    const scalar meanAlphaCoNum =
        0.5
       *(gSum(sumPhi)/gSum(mesh.V().primitiveField()))
       *runTime.deltaTValue();

    Info<< "Interface Courant Number mean: " << meanAlphaCoNum
        << " max: " << alphaCoNum << endl;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::solvers::multiphysicsVoF::read()
{
    twoPhaseSolver::read();

    const dictionary& alphaControls = mesh.solution().solverDict(alpha1.name());

    vDotResidualAlpha =
        alphaControls.lookupOrDefault("vDotResidualAlpha", 1e-4);

    return true;
}

void Foam::solvers::multiphysicsVoF::correctInterface()
{
    interface.correct();
}

Foam::tmp<Foam::surfaceScalarField>
Foam::solvers::multiphysicsVoF::surfaceTensionForce() const
{    
    return interface.surfaceTensionForce() * corrf;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::multiphysicsVoF::multiphysicsVoF(fvMesh& mesh)
:
    twoPhaseSolver
    (
        mesh,
        autoPtr<twoPhaseVoFMixture>(new compressibleTwoPhaseVoFMixture(mesh))
    ),

    mixture_
    (
        refCast<compressibleTwoPhaseVoFMixture>(twoPhaseSolver::mixture)
    ),

    interface(mixture_, alpha1, alpha2, U),

    sigmaPtr
    (
        surfaceTensionModel::New
        (
            IOdictionary
            (
                IOobject
                (
                    "phaseProperties",
                    runTime.constant(),
                    mesh,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            ), mesh
        )
    ),

    p(mixture_.p()),

    vDot
    (
        IOobject
        (
            "vDot",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        alpha1()*fvc::div(phi)()()
    ),

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

    K("K", 0.5*magSqr(U)),

    corrf
    (
        "corrf",
        fvc::interpolate(scalar(2) * mixture_.rho() / (mixture_.rho1() + mixture_.rho2()))
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

    thermophysicalTransport(momentumTransport),

    mixture(mixture_)
{
    const word alphaScheme(mesh.schemes().div(divAlphaName)[1].wordToken());

    if (!compressionSchemes.found(alphaScheme))
    {
        WarningInFunction
            << "Scheme " << alphaScheme << " for " << divAlphaName
            << " is not an interface compression scheme:"
            << compressionSchemes.toc() << endl;
    }

    const dictionary& alphaControls = mesh.solution().solverDict(alpha1.name());

    if (alphaControls.found("cAlpha"))
    {
        FatalErrorInFunction
            << "Deprecated and unused cAlpha entry specified in "
            << alphaControls.name() << nl
            << "Please update the case to use one of the run-time "
               "selectable interface compression schemes:"
            << compressionSchemes.toc() << exit(FatalError);
    }

    if (transient())
    {
        correctCoNum();
    }
    
    read();

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

Foam::solvers::multiphysicsVoF::~multiphysicsVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::multiphysicsVoF::prePredictor()
{
    twoPhaseSolver::prePredictor();

    const volScalarField& rho1 = mixture_.thermo1().rho();
    const volScalarField& rho2 = mixture_.thermo2().rho();

    alphaRhoPhi1 = fvc::interpolate(rho1)*alphaPhi1;
    alphaRhoPhi2 = fvc::interpolate(rho2)*alphaPhi2;

    rhoPhi = alphaRhoPhi1 + alphaRhoPhi2;

    corrf = fvc::interpolate(scalar(2) * mixture_.rho() / (rho1 + rho2));

    contErr1 =
    (
        fvc::ddt(alpha1, rho1)()() + fvc::div(alphaRhoPhi1)()()
      - (fvModels().source(alpha1, rho1)&rho1)()
    );

    contErr2 =
    (
        fvc::ddt(alpha2, rho2)()() + fvc::div(alphaRhoPhi2)()()
      - (fvModels().source(alpha2, rho2)&rho2)()
    );
}


void Foam::solvers::multiphysicsVoF::momentumTransportPredictor()
{
    momentumTransport.predict();
}


void Foam::solvers::multiphysicsVoF::thermophysicalTransportPredictor()
{
    thermophysicalTransport.predict();
}


void Foam::solvers::multiphysicsVoF::momentumTransportCorrector()
{
    momentumTransport.correct();
}


void Foam::solvers::multiphysicsVoF::thermophysicalTransportCorrector()
{
    thermophysicalTransport.correct();
}


// ************************************************************************* //
