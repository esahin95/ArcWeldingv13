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
    Foam::solvers::icoPhaseChangeVoF

Description
    Solver module for 2 incompressible, isothermal immiscible fluids using a
    VOF (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing.

    The momentum and other fluid properties are of the "mixture_" and a single
    momentum equation is solved.

    Either mixture_ or two-phase transport modelling may be selected.  In the
    mixture_ approach a single laminar, RAS or LES model is selected to model the
    momentum stress.  In the Euler-Euler two-phase approach separate laminar,
    RAS or LES selected models are selected for each of the phases.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

    Optional fvModels and fvConstraints are provided to enhance the simulation
    in many ways including adding various sources, Lagrangian
    particles, surface film etc. and constraining or limiting the solution.

SourceFiles
    icoPhaseChangeVoF.C

See also
    Foam::solvers::VoFSolver
    Foam::solvers::twoPhaseVoFSolver
    Foam::solvers::compressibleVoF

\*---------------------------------------------------------------------------*/

#ifndef icoPhaseChangeVoF_H
#define icoPhaseChangeVoF_H

#include "icoThermoVoF.H"
#include "solidificationModel.H"
#include "evaporationModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class icoPhaseChangeVoF Declaration
\*---------------------------------------------------------------------------*/

class icoPhaseChangeVoF
:
    public icoThermoVoF
{

    // Phase change models

        //- Solidification model
        autoPtr<solidificationModel> solidificationModel_;

        //- Evaporation model
        autoPtr<evaporationModel> evaporationModel_;


    // Numerical Controls

        // Thermo corrections
        const label nThermoCorr_;

        // Thermo tolerance
        const scalar thermoTol_;


    // Private member functions

        // Add source terms
        template<class Type>
        void addSup(fvMatrix<Type>& eqn) const;

        // Correct phase change models
        scalar correctPhaseChange();


public:

    //- Runtime type information
    TypeName("icoPhaseChangeVoF");


    // Constructors

        //- Construct from region mesh
        icoPhaseChangeVoF(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        icoPhaseChangeVoF(const icoPhaseChangeVoF&) = delete;


    //- Destructor
    virtual ~icoPhaseChangeVoF();


    // Member Functions

        //- Construct and optionally solve the momentum equation
        virtual void momentumPredictor();

        //- Construct and solve the energy equation
        virtual void thermophysicalPredictor();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const icoPhaseChangeVoF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "icoPhaseChangeVoFTemplates.C"
#endif


#endif

// ************************************************************************* //
