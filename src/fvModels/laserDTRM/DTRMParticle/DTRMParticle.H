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
    Simple solid spherical particle class with one-way coupling with the
    continuous phase.

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

class DTRMParticleCloud;


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

        //- Diameter
        scalar d_;

        //- Velocity of parcel
        vector U_;


public:

    friend class lagrangian::Cloud<DTRMParticle>;


    //- Class used to pass tracking data to the trackToFace function
    class trackingData
    :
        public particle::trackingData
    {
        // Interpolators for continuous phase fields

            const interpolationCellPoint<scalar>& rhoInterp_;
            const interpolationCellPoint<vector>& UInterp_;
            const interpolationCellPoint<scalar>& nuInterp_;

        //- Local gravitational or other body-force acceleration
        const vector& g_;


    public:

        // Constructors

            inline trackingData
            (
                const DTRMParticleCloud& spc,
                const interpolationCellPoint<scalar>& rhoInterp,
                const interpolationCellPoint<vector>& UInterp,
                const interpolationCellPoint<scalar>& nuInterp,
                const vector& g
            );


        // Member Functions

            inline const interpolationCellPoint<scalar>& rhoInterp() const;

            inline const interpolationCellPoint<vector>& UInterp() const;

            inline const interpolationCellPoint<scalar>& nuInterp() const;

            inline const vector& g() const;
    };


    //- Runtime type information
    TypeName("DTRMParticle");


    // Constructors

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

            //- Return diameter
            inline scalar d() const;

            //- Return velocity
            inline const vector& U() const;


        // Tracking

            //- Move
            bool move(DTRMParticleCloud&, trackingData&);


        // Patch interactions

            //- Overridable function to handle the particle hitting a wallPatch
            void hitWallPatch(DTRMParticleCloud& cloud, trackingData& td);

            //- Transform the physical properties of the particle
            //  according to the given transformation
            virtual void transformProperties(const transformer&);


    // I-O

        static void readFields(lagrangian::Cloud<DTRMParticle>& c);

        static void writeFields(const lagrangian::Cloud<DTRMParticle>& c);


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
