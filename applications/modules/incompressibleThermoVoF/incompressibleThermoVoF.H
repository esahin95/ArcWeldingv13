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

Class
    Foam::solvers::incompressibleThermoVoF

Description
    Solver module for 2 incompressible, isothermal immiscible fluids using a
    VOF (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    Either mixture or two-phase transport modelling may be selected.  In the
    mixture approach a single laminar, RAS or LES model is selected to model the
    momentum stress.  In the Euler-Euler two-phase approach separate laminar,
    RAS or LES selected models are selected for each of the phases.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, Lagrangian
    particles, surface film etc. and constraining or limiting the solution.

SourceFiles
    incompressibleThermoVoF.C

See also
    Foam::solvers::VoFSolver
    Foam::solvers::twoPhaseVoFSolver
    Foam::solvers::compressibleVoF

\*---------------------------------------------------------------------------*/

#ifndef incompressibleThermoVoF_H
#define incompressibleThermoVoF_H

#include "twoPhaseSolver.H"
#include "interfaceThermoProperties.H"

#include "incompressibleTwoPhaseVoFMixture.H"
#include "incompressibleInterPhaseTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class incompressibleThermoVoF Declaration
\*---------------------------------------------------------------------------*/

class incompressibleThermoVoF
:
    public twoPhaseSolver
{

protected:

    // Phase properties

        //- The compressible two-phase mixture
        incompressibleTwoPhaseVoFMixture& mixture;

    // Interface properties

        interfaceThermoProperties interface;

        surfaceScalarField corrf;

    // Thermophysical properties

        //- Static pressure field
        volScalarField p;

        //- Temperature field
        volScalarField T;

        //- Solid fraction field
        volScalarField alphaSolid;

        //- Solid fraction function
        autoPtr<Function1<scalar>> alphaSolidT_;

        //- Latent heat of fusion
        dimensionedScalar L_;

        //- Porosity
        dimensionedScalar Cu_;


    // Pressure reference

        //- Pressure reference
        Foam::pressureReference pressureReference_;


    // Momentum transport

        //- Momentum transport model
        incompressibleInterPhaseTransportModel momentumTransport;

    // Thermophysical transport

        //- Time coefficient
        volScalarField rhoCp;

        //- Convective transport
        surfaceScalarField rhoPhiCp;

private:

    // Private Member functions

        //- Correct the cached Courant numbers
        virtual void correctCoNum();

protected:

    // Protected Member Functions

        //- Adjust the rDeltaT in the vicinity of the interface
        virtual void setInterfaceRDeltaT(volScalarField& rDeltaT);

        //- Correct the interface properties following mesh-change
        //  and phase-fraction update
        virtual void correctInterface();

        //- Return the interface surface tension force for the momentum equation
        virtual tmp<surfaceScalarField> surfaceTensionForce() const;

        //- Return the pressure reference
        virtual const Foam::pressureReference& pressureReference() const
        {
            return pressureReference_;
        }

        //- The flow is incompressible
        virtual bool incompressible() const
        {
            return true;
        }

        //- Is the flow divergent?
        //  i.e. includes phase-fraction sources
        virtual bool divergent() const
        {
            return
                fvModels().addsSupToField(alpha1.name())
             || fvModels().addsSupToField(alpha2.name());
        }

        //- Return the mixture compressibility/density
        //  Not required for incompressible fluids
        virtual tmp<volScalarField> psiByRho() const
        {
            return tmp<volScalarField>(nullptr);
        }

        //- Calculate the alpha equation sources
        virtual void alphaSuSp
        (
            tmp<volScalarField::Internal>& Su,
            tmp<volScalarField::Internal>& Sp
        );

        //- Return the momentum equation stress term
        virtual tmp<fvVectorMatrix> divDevTau(volVectorField& U)
        {
            return momentumTransport.divDevTau(rho, U);
        }


public:

    //- Runtime type information
    TypeName("incompressibleThermoVoF");


    // Constructors

        //- Construct from region mesh
        incompressibleThermoVoF(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        incompressibleThermoVoF(const incompressibleThermoVoF&) = delete;


    //- Destructor
    virtual ~incompressibleThermoVoF();


    // Member Functions

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Predict the momentum transport
        virtual void momentumTransportPredictor();

        //- Predict thermophysical transport
        virtual void thermophysicalTransportPredictor();

        //- Construct and solve the energy equation,
        //  convert to temperature
        //  and update thermophysical and transport properties
        virtual void thermophysicalPredictor();

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the pressure equation in the PISO loop
        virtual void pressureCorrector();

        //- Correct the momentum transport
        virtual void momentumTransportCorrector();

        //- Correct the thermophysical transport
        virtual void thermophysicalTransportCorrector();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const incompressibleThermoVoF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
