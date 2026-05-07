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
    Foam::solvers::incompressibleReactingMultiphaseVoF

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
    incompressibleReactingMultiphaseVoF.C

See also
    Foam::solvers::VoFSolver
    Foam::solvers::multiphaseVoFSolver

\*---------------------------------------------------------------------------*/

#ifndef incompressibleReactingMultiphaseVoF_H
#define incompressibleReactingMultiphaseVoF_H

#include "compressibleMultiphaseVoF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class incompressibleReactingMultiphaseVoF Declaration
\*---------------------------------------------------------------------------*/

class incompressibleReactingMultiphaseVoF
:
    public compressibleMultiphaseVoF
{

protected:

    // Protected Data Members


private:

    // Private Member Functions



public:

    //- Runtime type information
    TypeName("incompressibleReactingMultiphaseVoF");


    // Constructors

        //- Construct from region mesh
        incompressibleReactingMultiphaseVoF(fvMesh& mesh);

        //- Disallow default bitwise copy construction
        incompressibleReactingMultiphaseVoF
        (
            const incompressibleReactingMultiphaseVoF&
        ) = delete;


    //- Destructor
    virtual ~incompressibleReactingMultiphaseVoF();


    // Member Functions




    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const incompressibleReactingMultiphaseVoF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
