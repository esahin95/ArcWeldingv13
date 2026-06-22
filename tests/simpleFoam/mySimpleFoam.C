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
    mySimpleFoam

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "volFields.H"
#include "fvmLaplacian.H"
using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;


    runTime++;

    const dimensionedScalar kappa
    (
        "kappa",
        dimPower/dimLength/dimTemperature,
        1.0
    );

    for(label i=0; i<3; i++)
    {
        fvScalarMatrix TEqn
        (
            fvm::laplacian(kappa, T)
        );

        TEqn.solve();
    }

    T.write();

    volVectorField gradT
    (
        IOobject
        (
            "gradT",
            runTime.name(),
            mesh
        ),
        fvc::grad(T)
    );

    gradT.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
