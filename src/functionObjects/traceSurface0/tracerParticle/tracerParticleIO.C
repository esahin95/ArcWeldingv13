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

\*---------------------------------------------------------------------------*/

#include "tracerParticle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::tracerParticle::sizeofFields_
(
    sizeof(tracerParticle) - sizeof(particle)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tracerParticle::tracerParticle(Istream& is, bool readFields)
:
    particle(is, readFields)
{
    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            q0_ = readScalar(is);
            is >> trackIndex_ >> q_ >> d_ >> transmissive_;
        }
        else
        {
            is.read(reinterpret_cast<char*>(&q0_), sizeofFields_);
        }
    }

    // Check state of Istream
    is.check("tracerParticle::tracerParticle(Istream&)");
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const tracerParticle& p)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.q0_
            << token::SPACE << p.trackIndex_
            << token::SPACE << p.q_
            << token::SPACE << p.d_
            << token::SPACE << p.transmissive_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.q0_),
            tracerParticle::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const tracerParticle&)");

    return os;
}


// ************************************************************************* //
