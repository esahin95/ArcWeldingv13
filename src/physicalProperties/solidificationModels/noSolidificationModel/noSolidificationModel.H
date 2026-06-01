/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026 OpenFOAM Foundation
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
    Foam::noSolidificationModel

Description

SourceFiles
    noSolidificationModel.C

\*---------------------------------------------------------------------------*/

#ifndef noSolidificationModel_H
#define noSolidificationModel_H

#include "solidificationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace solidificationModels
{


/*---------------------------------------------------------------------------*\
                          Class noSolidificationModel Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public solidificationModel
{

public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from components
        none(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        none(const none&) = delete;


    //- Destructor
    virtual ~none()
    {}

    // Member functions

        //- Correct solid fraction
        virtual void correct(const volScalarField& T);

        //- Source term for energy equation
        virtual tmp<fvScalarMatrix> hSource
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const surfaceScalarField& alphaRhoPhi,
            const volScalarField& T
        ) const;

        //- Source term for momentum equation
        virtual tmp<fvVectorMatrix> USource
        (
            const volScalarField& alpha,
            const volScalarField& rho,
            const volVectorField& U
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const none&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidificationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
