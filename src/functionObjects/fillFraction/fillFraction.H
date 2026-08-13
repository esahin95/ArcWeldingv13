/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2021 OpenFOAM Foundation
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
    Foam::functionObjects::fillFraction

Description
    This function object reports the height of the interface above a set of
    locations. For each location, it writes the vertical distance of the
    interface above both the location and the lowest boundary. It also writes
    the point on the interface from which these heights are computed. It uses
    an integral approach, so if there are multiple interfaces above or below a
    location then this method will generate average values.

    Example of function object specification:
    \verbatim
    fillFraction1
    {
        type           fillFraction;
        libs           ("libfieldFunctionObjects.so");
        alpha          alpha.water;
        locations      ((0 0 0) (10 0 0) (20 0 0));
    }
    \endverbatim

Usage
    \table
        Property     | Description               | Required | Default value
        type         | type name                 | yes      |
        alpha        | name of the alpha field   | no       | alpha
        locations    | list of locations to report the height at | yes |
        liquid       | is the alpha field that of the liquid | no | true
    \endtable

SourceFiles
    fillFraction.C

\*---------------------------------------------------------------------------*/

#ifndef fillFraction_H
#define fillFraction_H

#include "traceSurface.H"

//#include "Cloud.H"
//#include "tracerParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class fillFraction Declaration
\*---------------------------------------------------------------------------*/

class fillFraction
:
    public traceSurface
{

protected:

    // Protected Member Functions

        //- Output quantity
        virtual scalar output(lagrangian::Cloud<tracerParticle>& cloud);


public:

    //- Runtime type information
    TypeName("fillFraction");


    // Constructors

        //- Construct from Time and dictionary
        fillFraction
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~fillFraction();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
