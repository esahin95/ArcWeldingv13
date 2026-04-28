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
    Foam::noReflection

Description

SourceFiles
    noReflectionI.H
    noReflection.C
    noReflectionIO.C

\*---------------------------------------------------------------------------*/

#ifndef noReflection_H
#define noReflection_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class noReflection;
Istream& operator>>(Istream&, noReflection&);
Ostream& operator<<(Ostream&, const noReflection&);


/*---------------------------------------------------------------------------*\
                          Class noReflection Declaration
\*---------------------------------------------------------------------------*/

class noReflection
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
        noReflection();

        //- Construct from components
        noReflection(const dataType& data);

        //- Construct from Istream
        noReflection(Istream&);

        //- Copy constructor
        noReflection(const noReflection&);

        //- Disallow default bitwise copy construction
        noReflection(const noReflection&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<noReflection> New();


    //- Destructor
    ~noReflection();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const noReflection&);

        //- Disallow default bitwise assignment
        void operator=(const noReflection&) = delete;


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, noReflection&);
        friend Ostream& operator<<(Ostream&, const noReflection&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "noReflectionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
