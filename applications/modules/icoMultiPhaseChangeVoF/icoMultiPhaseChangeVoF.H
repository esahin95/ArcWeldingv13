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

Class
    Foam::solvers::icoMultiPhaseChangeVoF

Description
    Solver module for the solution of multiple compressible, isothermal
    immiscible fluids using a VOF (volume of fluid) phase-fraction based
    interface capturing approach, with optional mesh motion and mesh topology
    changes including adaptive re-meshing.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    A mixture approach for momentum transport is provided in which a single
    laminar, RAS or LES model is selected to model the momentum stress.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient and steady simulations.

SourceFiles
    icoMultiPhaseChangeVoF.C

See also
    Foam::solvers::VoFSolver
    Foam::solvers::multiphaseVoFSolver

\*---------------------------------------------------------------------------*/

#ifndef icoMultiPhaseChangeVoF_H
#define icoMultiPhaseChangeVoF_H

#include "icoMulticomponentVoF.H"

#include "solidificationModel.H"
#include "evaporationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class icoMultiPhaseChangeVoF Declaration
\*---------------------------------------------------------------------------*/

class icoMultiPhaseChangeVoF
:
    public icoMulticomponentVoF
{
protected:

    // Phase change models

        // Solidification model pointers
        PtrList<solidificationModel> solModels_;

        // Evaporation model pointers
        PtrList<evaporationModel> evaModels_;

public:

    //- Runtime type information
    TypeName("icoMultiPhaseChangeVoF");


    // Constructors

        //- Construct from region mesh
        icoMultiPhaseChangeVoF(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        icoMultiPhaseChangeVoF
        (
            const icoMultiPhaseChangeVoF&
        ) = delete;


    //- Destructor
    virtual ~icoMultiPhaseChangeVoF();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const icoMultiPhaseChangeVoF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
