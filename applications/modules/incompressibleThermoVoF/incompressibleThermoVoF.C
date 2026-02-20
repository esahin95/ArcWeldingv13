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

#include "incompressibleThermoVoF.H"
#include "localEulerDdtScheme.H"
#include "fvCorrectPhi.H"
#include "geometricZeroField.H"
#include "addToRunTimeSelectionTable.H"

#include "fvcAverage.H"
#include "interfaceCompression.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(incompressibleThermoVoF, 0);
    addToRunTimeSelectionTable(solver, incompressibleThermoVoF, fvMesh);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::solvers::incompressibleThermoVoF::correctCoNum()
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

void Foam::solvers::incompressibleThermoVoF::correctInterface()
{
    interface.correct();
}


Foam::tmp<Foam::surfaceScalarField>
Foam::solvers::incompressibleThermoVoF::surfaceTensionForce() const
{
    return interface.surfaceTensionForce() * corrf;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::incompressibleThermoVoF::incompressibleThermoVoF(fvMesh& mesh)
:
    twoPhaseSolver
    (
        mesh,
        autoPtr<twoPhaseVoFMixture>(new incompressibleTwoPhaseVoFMixture(mesh))
    ),

    mixture
    (
        refCast<incompressibleTwoPhaseVoFMixture>(twoPhaseSolver::mixture)
    ),

    interface(mixture, alpha1, alpha2, U),

    corrf
    (
        IOobject
        (
            "corrf",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        scalar(2.0) * fvc::interpolate(mixture.rho()) / (mixture.rho1() + mixture.rho2())
    ),

    p
    (
        IOobject
        (
            "p",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rho*buoyancy.gh
    ),

    T
    (
        IOobject
        (
            "T",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),

    alphaSolid
    (
        IOobject
        (
            "alphaSolid",
            runTime.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),

    alphaSolidT_(),

    L_("L", dimensionSet(0, 2, -2, 0, 0), mixture.lookupOrDefault("L", 0.0)),

    Cu_("Cu", dimensionSet(0, 0, -1, 0, 0), mixture.lookupOrDefault("Cu", 0.0)),

    pressureReference_
    (
        p,
        p_rgh,
        pimple.dict()
    ),

    momentumTransport
    (
        U,
        phi,
        alphaPhi1,
        alphaPhi2,
        mixture
    ),

    rhoCp
    (
        IOobject
        (
            "rhoCp",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mixture.rho() * mixture.Cp() 
    ),

    rhoPhiCp
    (
        IOobject
        (
            "rhoPhiCp",
            runTime.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        phi * fvc::interpolate(rhoCp)
    )
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

    alphaSolidT_.reset
    (
        Function1<scalar>::New
        (
            "alphaSolidT",
            dimTemperature,
            unitFraction,
            refCast<dictionary>(mixture)
        ).ptr()
    );

    alphaSolid.primitiveFieldRef() = alpha1 * alphaSolidT_->value(T.primitiveField());
    alphaSolid.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::incompressibleThermoVoF::~incompressibleThermoVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleThermoVoF::prePredictor()
{    
    // Correct solid fraction field
    scalar relax(0.9);
    alphaSolid.storePrevIter();
    alphaSolid.primitiveFieldRef() = min
    (
        relax * alpha1 * alphaSolidT_->value(T.primitiveField()) + (1-relax) * alphaSolid.primitiveField(),
        alpha1
    );
    alphaSolid.correctBoundaryConditions();
    Info << "Maximum change in solid fraction: " << gMax(mag(alphaSolid.prevIter().primitiveField()-alphaSolid.primitiveField())) << endl;
    
    twoPhaseSolver::prePredictor();

    const dimensionedScalar& rho1 = mixture.rho1();
    const dimensionedScalar& rho2 = mixture.rho2();

    // Calculate the mass-flux
    rhoPhi = alphaPhi1*rho1 + alphaPhi2*rho2;

    // Correction factor
    corrf = fvc::interpolate(rho) * (scalar(2) / (rho1 + rho2));
}


void Foam::solvers::incompressibleThermoVoF::momentumTransportPredictor()
{
    momentumTransport.predict();
}


void Foam::solvers::incompressibleThermoVoF::thermophysicalTransportPredictor()
{
    //volScalarField Cp(mixture.Cp());
    
    //rhoCp == rho * Cp;

    //rhoPhiCp = rhoPhi * fvc::interpolate(Cp);
}

void Foam::solvers::incompressibleThermoVoF::momentumTransportCorrector()
{
    momentumTransport.correct();
}


void Foam::solvers::incompressibleThermoVoF::thermophysicalTransportCorrector()
{}


// ************************************************************************* //
