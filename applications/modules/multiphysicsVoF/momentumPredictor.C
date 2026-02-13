/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "multiphysicsVoF.H"
#include "fvmSup.H"
#include "fvmDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcReconstruct.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::solvers::multiphysicsVoF::divDevTau
(
    volVectorField& U
)
{
    return
        momentumTransport.divDevTau(U)
      - fvm::Sp(contErr1() + contErr2(), U);
}


void Foam::solvers::multiphysicsVoF::momentumPredictor()
{    
    volVectorField& U = U_;

    volVectorField gradSigma(fvc::grad(sigmaPtr->sigma()));
    volVectorField nHatv(interface.n());
    volScalarField corr((scalar(2) * rho) / (mixture.rho1() + mixture.rho2()));

    tUEqn =
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + MRF.DDt(rho, U)
      + divDevTau(U)
     ==
        fvModels().source(rho, U)
      + (gradSigma - (gradSigma & nHatv) * nHatv) * mag(fvc::grad(alpha1)) * corr
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvConstraints().constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    surfaceTensionForce()
                  - buoyancy.ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvConstraints().constrain(U);

        K = 0.5*magSqr(U);
    }
}


// ************************************************************************* //
