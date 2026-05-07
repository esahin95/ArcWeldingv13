/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "incompressibleReactingVoF.H"
#include "fvcMeshPhi.H"
#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvmLaplacian.H"
#include "fvcVolumeIntegrate.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::incompressibleReactingVoF::thermophysicalPredictor()
{
    const volScalarField& rho1(mixture.rho1());
    const volScalarField& rho2(mixture.rho2());

    volScalarField& T = mixture_.T();

    const volScalarField& e1(mixture.thermo1().he());
    const volScalarField& e2(mixture.thermo2().he());

    const fvScalarMatrix e1Source(fvModels().source(alpha1, rho1, e1));
    const fvScalarMatrix e2Source(fvModels().source(alpha2, rho2, e2));

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp_,T)
        + fvm::div(rhoPhiCp_, T)
        - fvm::Sp(fvc::ddt(rhoCp_) + fvc::div(rhoPhiCp_), T)
        - fvm::laplacian(thermophysicalTransport.kappaEff(), T)
        ==
          (e1Source&e1) + (e2Source&e2)
    );

    TEqn.relax();

    fvConstraints().constrain(TEqn);

    label iter = 0;
    scalar res;
    do
    {
        solidFraction_.storePrevIter();

        solve
        (
            TEqn
          ==
            L_ *
            (
                fvc::ddt(alpha1, rho1, solidFraction_)
              + fvc::div(alphaRhoPhi1, solidFraction_)
            )
        );

        fvConstraints().constrain(T);

        solidFraction_ =
            max
            (
                min
                (
                    solidFraction_
                    + (relax_ / L_) * mixture_.thermo1().Cp()
                    * (Tliq_ - T - solidFraction_ * (Tliq_ - Tsol_)),
                    1.0
                ),
                0.0
            );

        res =
        gMax
        (
            mag
            (
                  solidFraction_.primitiveField()
                - solidFraction_.prevIter().primitiveField()
            )
        );

        Info<< gMin(T.primitiveField()) << " , "
            << gMax(T.primitiveField()) << " , "
            << res << endl;
    }
    while
    (
        ++iter < 50 && res > 1e-3
    );

    Info<< "Iter = " << iter << " , res = " << res << endl;

    alphaSolid_ = solidFraction_ * alphaVoF_;

    mixture_.correctThermo();
    mixture_.correct();
}


// ************************************************************************* //
