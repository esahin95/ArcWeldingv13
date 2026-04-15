/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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
    Foam::DTRMParticle

Description
    Laser particle class with coupling to the phase field.

SourceFiles
    DTRMParticleI.H
    DTRMParticle.C
    DTRMParticleIO.C

\*---------------------------------------------------------------------------*/

#ifndef DTRMParticle_H
#define DTRMParticle_H

#include "particle.H"
#include "Cloud.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class DTRMParticle;

Ostream& operator<<(Ostream&, const DTRMParticle&);


/*---------------------------------------------------------------------------*\
                        Class DTRMParticle Declaration
\*---------------------------------------------------------------------------*/

class DTRMParticle
:
    public particle
{
    // Private Data

        //- Size in bytes of the fields
        static const std::size_t sizeofFields_;

        /*
        //- Initial power
        const scalar q0_;

        //- Track index
        const label trackIndex_;

        //- Absorption rate at reflection
        const scalar a_;
        */

        //- Current power
        scalar q_;

        //- Direction
        vector d_;


public:

    friend class lagrangian::Cloud<DTRMParticle>;


    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Interpolators for continuous phase fields

            const interpolationCellPoint<scalar>& alphaInterp_;

        // Reference to heat source

            volScalarField& Q_;

    public:

        // Constructors

            inline trackingData
            (
                const lagrangian::Cloud<DTRMParticle>& cloud,
                const interpolationCellPoint<scalar>& alphaInterp_,
                volScalarField& Q
            );


        // Member Functions

            inline const interpolationCellPoint<scalar>& alphaInterp() const;
    };


    //- Runtime type information
    TypeName("DTRMParticle");


    // Constructors

        //- Construct from components
        /*
        DTRMParticle
        (
            const polyMesh& mesh,
            const vector& position,
            const vector& direction,
            const scalar power,
            const scalar absorption,
            const label trackIndex
        );
        */

        //- Construct from Istream
        DTRMParticle(Istream& is, bool readFields = true);

        //- Construct and return a clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new DTRMParticle(*this));
        }

        //- Construct from Istream and return
        static autoPtr<DTRMParticle> New(Istream& is)
        {
            return autoPtr<DTRMParticle>(new DTRMParticle(is));
        }


    // Member Functions

        // Access

            //- Return power
            inline scalar q() const;

            //- Return direction
            inline const vector& d() const;


        // Tracking

            //- Move
            bool move(lagrangian::Cloud<DTRMParticle>&, trackingData&);


        // Patch interactions

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(lagrangian::Cloud<DTRMParticle>& cloud, trackingData& td);


    // I-O

        //static void readFields(lagrangian::Cloud<DTRMParticle>& c);

        //static void writeFields(const lagrangian::Cloud<DTRMParticle>& c);


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const DTRMParticle&);
};


template<>
inline bool contiguous<DTRMParticle>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DTRMParticleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
