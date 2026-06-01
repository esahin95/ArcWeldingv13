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
    noSolidificationModelI.H
    noSolidificationModel.C
    noSolidificationModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef noSolidificationModel_H
#define noSolidificationModel_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class noSolidificationModel;
Istream& operator>>(Istream&, noSolidificationModel&);
Ostream& operator<<(Ostream&, const noSolidificationModel&);


/*---------------------------------------------------------------------------*\
                          Class noSolidificationModel Declaration
\*---------------------------------------------------------------------------*/

class noSolidificationModel
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
        noSolidificationModel();

        //- Construct from components
        noSolidificationModel(const dataType& data);

        //- Construct from Istream
        noSolidificationModel(Istream&);

        //- Copy constructor
        noSolidificationModel(const noSolidificationModel&);

        //- Disallow default bitwise copy construction
        noSolidificationModel(const noSolidificationModel&) = delete;


    // Selectors

        //- Select null constructed
        static autoPtr<noSolidificationModel> New();


    //- Destructor
    ~noSolidificationModel();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const noSolidificationModel&);

        //- Disallow default bitwise assignment
        void operator=(const noSolidificationModel&) = delete;


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, noSolidificationModel&);
        friend Ostream& operator<<(Ostream&, const noSolidificationModel&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "noSolidificationModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
