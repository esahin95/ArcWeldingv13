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
    Foam::fresnelReflection

Description

SourceFiles
    fresnelReflectionI.H
    fresnelReflection.C
    fresnelReflectionIO.C

\*---------------------------------------------------------------------------*/

#ifndef fresnelReflection_H
#define fresnelReflection_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class fresnelReflection;
Istream& operator>>(Istream&, fresnelReflection&);
Ostream& operator<<(Ostream&, const fresnelReflection&);


/*---------------------------------------------------------------------------*\
                          Class fresnelReflection Declaration
\*---------------------------------------------------------------------------*/

class fresnelReflection
:
    public baseClassName
{
    // Private Data

        //- Description of data_
        dataType data_;


    // Private Member Functions


public:

    // Static Data Members

        //- Static data staticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        fresnelReflection();

        //- Construct from components
        fresnelReflection(const dataType& data);

        //- Construct from Istream
        fresnelReflection(Istream&);

        //- Copy constructor
        fresnelReflection(const fresnelReflection&);

        //- Disallow default bitwise copy construction
        fresnelReflection(const fresnelReflection&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<fresnelReflection> New();


    //- Destructor
    ~fresnelReflection();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const fresnelReflection&);

        //- Disallow default bitwise assignment
        void operator=(const fresnelReflection&) = delete;


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, fresnelReflection&);
        friend Ostream& operator<<(Ostream&, const fresnelReflection&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fresnelReflectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
