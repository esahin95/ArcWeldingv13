/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2023 OpenFOAM Foundation
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
    Foam::fv::laserDTRM

Description
    Discrete Tray Radiation Method for collimated radiation flux, derived
    from the implementation given in OpenFOAM-v2512.

    The ray is considered planar and developed to work with a solver
    derived from compressibleVoF.

Usage
    Example usage:
    \verbatim
    laserDTRM
    {
        type    laserDTRM;
        phase   water;
        ***Add more stuff here***
    }
    \endverbatim

SourceFiles
    laserDTRM.C

\*---------------------------------------------------------------------------*/

#ifndef laserDTRM_H
#define laserDTRM_H

#include "fvModel.H"
#include "fluidThermo.H"

#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class laserDTRM Declaration
\*---------------------------------------------------------------------------*/

class laserDTRM
:
    public fvModel
{
    // Private Data

        //- The name of the Lagrangian phase
        const word phaseName_;

        //- Reference to the carrier phase thermo
        const fluidThermo& thermo_;

        //- Reference to the carrier phase field
        const volScalarField& alpha_;

        //- Total laser power
        const scalar Qtot_;

        //- Number of rays
        const label nRays_;

        //- Laser position
        const autoPtr<Function1<vector>> pos_;

        //- Laser radius
        const scalar rad_;

        //- Laser normal direction
        const vector normal_;

        //- Derived radial directions
        const vector radial1_;
        const vector radial2_;

        //- Radial power distribution
        const autoPtr<Function1<scalar>> powerDist_;

        //- Relaxation factor
        const scalar relax_;

        //- Absorption rate
        const scalar a_;

        //- Laser power deposition field
        mutable volScalarField Q_;

        //- Current time index (used for updating)
        mutable label curTimeIndex_;


public:

    //- Runtime type information
    TypeName("laserDTRM");


    // Constructors

        //- Construct from explicit source name and mesh
        laserDTRM
        (
            const word& sourceName,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        laserDTRM
        (
            const laserDTRM&
        ) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Correct

            //- Ray tracing and updating power deposition fields
            virtual void correct();


        // Add explicit and implicit contributions to compressible equation

            //- Add explicit contribution to phase energy equation
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& he,
                fvMatrix<scalar>& eqn
            ) const;


        // Mesh changes

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const laserDTRM&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
