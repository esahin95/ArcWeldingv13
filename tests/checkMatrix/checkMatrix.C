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
    checkMatrix

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "volFields.H"
#include "fvcFlux.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "simpleMatrix.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // Set root case
    argList args(argc, argv);
    if (!args.checkRootCase())
    {
        FatalError.exit();
    }

    // Create time
    Info<< "Create time" << nl << endl;
    Time runTime(Time::controlDictName, args);

    // Create mesh
    Info<< "Create mesh for time = " << runTime.name() << nl << endl;
    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.name(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // Create fields
    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    T.correctBoundaryConditions();
    const volScalarField::Boundary& TBoundary = T.boundaryField();
    forAll(TBoundary, patchI)
    {
        if (TBoundary[patchI].patch().name() == "left")
        {
            const fvPatchField<scalar>& TPatch = TBoundary[patchI];
            Info<<TPatch.patch().name()<<endl;
            forAll(TPatch, faceI)
            {
                Info<<TPatch[faceI]<<endl;
            }
        }
    }

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.name(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    forAll(U, cellI)
    {
        const vector& C = mesh.C()[cellI];
        U[cellI] = vector(C.x(),-C.y(),0);
    }

    volVectorField::Boundary& UBoundary = U.boundaryFieldRef();
    forAll(UBoundary, patchI)
    {
        fvPatchField<vector>& UPatch = UBoundary[patchI];
        const Field<vector>& Cf = UPatch.patch().Cf();
        forAll(UPatch, faceI)
        {
            const vector& C = Cf[faceI];
            UPatch[faceI] = vector(C.x(), -C.y(), 0);
        }
    }

    // Setup linear equation system
    const surfaceScalarField phi
    (
        "phi",
        fvc::flux(U)
    );

    const dimensionedScalar D
    (
        "D",
        sqr(dimLength)/dimTime,
        0.01
    );


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    {
        runTime++;

        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        const dimensionedScalar Q(dimless/dimTime, 100.0);
        fvScalarMatrix TEqn
        (
            fvm::div(phi, T) - fvm::laplacian(D, T) == Q*0.0
        );


        {
            scalarField& src = TEqn.source();
            scalarField& diag = TEqn.diag();
            forAll(src, cellI)
            {
                const scalar Vc = mesh.V()[cellI];

                src[cellI] += Vc*Q.value();
                diag[cellI] += Vc*Q.value();
            }
        }

        volScalarField su(Q*(0.0*T + 1.0));
        TEqn -= fvm::Sp(Q, T) - su;
        TEqn += fvm::Sp(su.v(), T);

        volScalarField::Internal TI(T);
        fvScalarMatrix XEqn
        (
            fvm::ddt(su, T)
        );

        // Adding a source term on the right of TEqn is the same as adding it
        // to source or subtracting from TEqn


        #include "writeMatrixCoeff.H"

        TEqn.solve();

        T.write();
        U.write();
    }


    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
