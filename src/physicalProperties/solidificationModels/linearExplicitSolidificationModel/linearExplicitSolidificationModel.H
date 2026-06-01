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
    Foam::linearExplicitSolidificationModel

Description

SourceFiles
    linearExplicitSolidificationModel.C

\*---------------------------------------------------------------------------*/

#ifndef linearExplicitSolidificationModel_H
#define linearExplicitSolidificationModel_H

#include "solidificationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace solidificationModels
{


/*---------------------------------------------------------------------------*\
                          Class linearExplicitSolidificationModel Declaration
\*---------------------------------------------------------------------------*/

class linearExplicit
:
    public solidificationModel
{
protected:

    // Solidification properties

        dimensionedScalar Lm_;
        dimensionedScalar Tliq_;
        dimensionedScalar Tsol_;
        dimensionedScalar Cu_;
        scalar q_;
        scalar relax_;


public:

    //- Runtime type information
    TypeName("linearExplicit");

    // Constructors

        //- Construct from components
        linearExplicit(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        linearExplicit(const linearExplicit&) = delete;


    //- Destructor
    virtual ~linearExplicit()
    {}


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const linearExplicit&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidificationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
