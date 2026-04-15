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
    Foam::solvers::compressibleReactingVoF

Description
    Solver module for 2 compressible, non-isothermal immiscible fluids
    using a VOF (volume of fluid) phase-fraction based interface capturing
    approach, with optional mesh motion and mesh topology changes including
    adaptive re-meshing.

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
    compressibleReactingVoF.C

See also
    Foam::solvers::VoFSolver
    Foam::solvers::twoPhaseVoFSolver
    Foam::solvers::incompressibleReactingVoF

\*---------------------------------------------------------------------------*/

#ifndef compressibleReactingVoF_H
#define compressibleReactingVoF_H

#include "compressibleVoF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class compressibleReactingVoF Declaration
\*---------------------------------------------------------------------------*/

class compressibleReactingVoF
:
    public compressibleVoF
{

protected:

    // Phase change properties

        //- Properties
        const IOdictionary dict_;

        //- Active phase name
        const word phase_;

        //- Active volume fraction field
        const volScalarField& alphaVoF_;

        //- Solid fraction
        volScalarField solidFraction_;

        //- Solid phase fraction
        volScalarField alphaSolid_;

        //- Latent heat of fusion [J/kg]
        const dimensionedScalar L_;

        //- Solidus temperature
        const dimensionedScalar Tsol_;

        //- Liquidus temperature
        const dimensionedScalar Tliq_;

        //- Phase fraction under-relaxation coefficient
        const scalar relax_;

        //- Mushy region momentum sink coefficient [1/s]
        const dimensionedScalar Cu_;

        //- Coefficient used in porosity calc
        const scalar q_;


    // Protected Member Functions


public:

    // Access to the physical state


public:

    //- Runtime type information
    TypeName("compressibleReactingVoF");


    // Constructors

        //- Construct from region mesh
        compressibleReactingVoF(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        compressibleReactingVoF(const compressibleReactingVoF&) = delete;


    //- Destructor
    virtual ~compressibleReactingVoF();


    // Member Functions

        //- Called at the start of the PIMPLE loop
        virtual void prePredictor();

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the energy equation
        virtual void thermophysicalPredictor();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const compressibleReactingVoF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
