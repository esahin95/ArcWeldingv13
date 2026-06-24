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
    useGraphs

Description
    Example plain (non-CFD) application, which can be compiled and tested with
    the following commands, which should print "42" followed by "true".

        useGraphs 42 guess && echo true
        useGraphs -multiply 4 3 number

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "volFields.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Minimal default options and do not print header when running
    #include "removeCaseOptions.H"
    writeInfoHeader = false;

    // Character width of the options listed by "-help"
    argList::usageMin = 24;

    /*
    // Example "-multiply <factor>" command line option
    argList::addOption
    (
        "multiply",
        "factor",
        "multiply the scalar by the factor"
    );
    */

    // Example mandatory arguments
    argList::validArgs.append("x");
    argList::validArgs.append("y");

    /*
    // Example description for "-help"
    argList::addNote
    (
        "Prints '<name> = <scalar>*<factor>', where <factor> is optional\n"
        "Returns true if the word == 'guess'"
    );
    */

    argList args(argc, argv);

    /*
    const scalar m = args.optionLookupOrDefault<scalar>("multiply", 1.0);
    */
    const scalar x = args.argRead<scalar>(1);
    const scalar y = args.argRead<scalar>(2);

    Info<< x << " , " << y << endl;

    scalarField yvals(4, 0.0);
    yvals[0] = 1.0;
    yvals[1] = 0.7;
    yvals[2] = 0.3;
    yvals[3] = 0.0;

    label n = yvals.size();

    scalarField xvals(4, 0.0);
    xvals[0] = 1410;
    xvals[3] = 1620;
    xvals[1] = (xvals[n-1]-xvals[0])*(1. - yvals[1])+xvals[0];
    xvals[2] = (xvals[n-1]-xvals[0])*(1. - yvals[2])+xvals[0];
    Info<<xvals<<endl;

    const scalar slope = 1e10;

    const dimensionedScalar Tsol("TSol", dimTemperature, xvals[0]);
    const dimensionedScalar Tliq("Tliq", dimTemperature, xvals[n-1]);

    const scalar thermoTol = 1e-6;

    scalar T0;
    scalar dFdT;

    if ((x < Tliq.value()) && (x > Tsol.value()))
    {
        label lo = 0;
        for (lo=0; lo<n && xvals[lo]>x; ++lo)
        {}
        Info<<lo<<endl;

        label low = lo;
        if (low < n)
        {
            for (label i=low; i<n; ++i)
            {
                if (xvals[i] > xvals[lo] && xvals[i] <= x)
                {
                    lo = i;
                }
            }
        }
        Info<<lo<<endl;

        label hi = 0;
        for (hi=0; hi<n && xvals[hi]<x; ++hi)
        {}
        Info<<hi<<endl;

        label high = hi;
        if (high < n)
        {
            for (label i=high; i<n; ++i)
            {
                if (xvals[i] < xvals[hi] && xvals[i] >= x)
                {
                    hi = i;
                }
            }
        }
        Info<<hi<<endl;

        dFdT = (yvals[hi] - yvals[lo]) / (xvals[hi] - xvals[lo]);
        Info<< (-1.0/(Tliq-Tsol).value()) << endl;

        T0 = xvals[lo] + (y - yvals[lo])/dFdT;
        Info<< (Tsol.value() + (1.0 - y)*(Tliq-Tsol).value()) << endl;
    }
    else if ((x >= Tliq.value()) && (y > thermoTol))
    {
        dFdT = slope;
        T0 = Tliq.value() - thermoTol;
    }
    else if ((x <= Tsol.value()) && (y < 1 - thermoTol))
    {
        dFdT = slope;
        T0 = Tsol.value() + thermoTol;
    }
    else
    {
        dFdT = 0.0;
        T0 = x;
    }

    Info<< "T0 = " << T0 << " , dFdT = " << dFdT << endl;

    return 0.0;
}

// ************************************************************************* //
