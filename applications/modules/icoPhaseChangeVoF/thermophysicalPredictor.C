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

\*---------------------------------------------------------------------------*/

#include "icoPhaseChangeVoF.H"
#include "fvcMeshPhi.H"
#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


void Foam::solvers::icoPhaseChangeVoF::thermophysicalPredictor()
{
    volScalarField& T = mixture_.T();

    scalar resT, resSf;
    label iter = 0;
    do
    {
        T.storePrevIter();

        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp,T) + fvm::div(rhoPhiCp, T)
            - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCp), T)
            - fvm::laplacian(thermophysicalTransport.kappaEff(), T)
            ==
            fvModels().source(rhoCp, T)
        );
        solidificationModel_->hSource(TEqn);

        TEqn.relax();

        fvConstraints().constrain(TEqn);

        solve(TEqn);

        fvConstraints().constrain(T);

        resSf = solidificationModel_->correct();

        resT = gMax(mag(T.prevIter() - T)().primitiveField());
        Info<< resT << " , " << resSf << endl;
    }
    while
    (
           ++iter<50 && resSf > 0.001
    );

    mixture_.correctThermo();
    mixture_.correct();
}


// ************************************************************************* //
