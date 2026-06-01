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
    Foam::solidificationModel

Description

SourceFiles
    solidificationModelI.H
    solidificationModel.C
    solidificationModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef solidificationModel_H
#define solidificationModel_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class solidificationModel;
Istream& operator>>(Istream&, solidificationModel&);
Ostream& operator<<(Ostream&, const solidificationModel&);


/*---------------------------------------------------------------------------*\
                          Class solidificationModel Declaration
\*---------------------------------------------------------------------------*/

class solidificationModel
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
        solidificationModel();

        //- Construct from components
        solidificationModel(const dataType& data);

        //- Construct from Istream
        solidificationModel(Istream&);

        //- Copy constructor
        solidificationModel(const solidificationModel&);

        //- Disallow default bitwise copy construction
        solidificationModel(const solidificationModel&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<solidificationModel> New();


    //- Destructor
    ~solidificationModel();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const solidificationModel&);

        //- Disallow default bitwise assignment
        void operator=(const solidificationModel&) = delete;


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, solidificationModel&);
        friend Ostream& operator<<(Ostream&, const solidificationModel&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidificationModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
