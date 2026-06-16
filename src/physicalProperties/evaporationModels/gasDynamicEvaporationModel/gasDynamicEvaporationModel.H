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
    Foam::gasDynamicEvaporationModel

Description

SourceFiles
    gasDynamicEvaporationModel.C

\*---------------------------------------------------------------------------*/

#ifndef gasDynamicEvaporationModel_H
#define gasDynamicEvaporationModel_H

#include "evaporationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace evaporationModels
{


/*---------------------------------------------------------------------------*\
                          Class gasDynamicEvaporationModel Declaration
\*---------------------------------------------------------------------------*/

class gasDynamic
:
    public evaporationModel
{
protected:

    // Evaporation properties

        //- Evaporation properties dictionary
        const dictionary& dict_;

        //- Latent heat of evaporation
        const dimensionedScalar Lv_;

        //- Atmospheric properties
        const scalar p0_;
        const scalar T0_;
        const scalar R0_;
        const scalar g0_;

        //- Vapor properties
        const scalar Rv_;
        const scalar gv_;
        const scalar Tv_;

    // Numeric parameters

        // Relaxation
        const scalar relax_;

        // Stabilization term
        const scalar q_;


    // Precomputed factors

        // Mass flux coefficient
        scalar coeffA_;

        // Recoil pressure coefficient
        scalar coeffB_;


    //- Recoil pressure force
    volVectorField pRec_;


public:

    //- Runtime type information
    TypeName("gasDynamic");

    // Constructors

        //- Construct from components
        gasDynamic(const fvMesh& mesh, const word& group);

        //- Disallow default bitwise copy construction
        gasDynamic(const gasDynamic&) = delete;


    //- Destructor
    virtual ~gasDynamic()
    {}


    // Member functions

        //- Correct solid fraction
        virtual scalar correct(const bool relax=true);

        //- Source term for energy equation
        virtual void hSource(fvScalarMatrix& TEqn) const;

        //- Source term for momentum equation
        virtual void USource(fvVectorMatrix& UEqn) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const gasDynamic&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace evaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
