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

Application
    inheritance

Description
    Example plain (non-CFD) application, which can be compiled and tested with
    the following commands, which should print "42" followed by "true".

        inheritance 42 guess && echo true
        inheritance -multiply 4 3 number

\*---------------------------------------------------------------------------*/

#include "argList.H"

#include "derived.H"

#include <limits>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Minimal default options and do not print header when running
    #include "removeCaseOptions.H"
    writeInfoHeader = false;

    argList args(argc, argv);

    derived d(0.0);

    d.method();

    Info<<small<<endl;
    Info<<great<<endl;
    Info<<(1.0/great)<<endl;

    double x = std::numeric_limits<double>::epsilon();
    Info<<x<<endl;

    return 0;
}

// ************************************************************************* //
