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

#include "KnightEvaporationModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace evaporationModels
{
    defineTypeNameAndDebug(Knight, 0);

    addToRunTimeSelectionTable
    (
        evaporationModel,
        Knight,
        dictionary
    );
}
}

using namespace Foam::constant;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::evaporationModels::Knight::Knight
(
    const fvMesh& mesh,
    const word& group
)
:
    gasDynamic(mesh, group),

    Th_(dict_.lookup<scalar>("Th"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::evaporationModels::Knight::correct(const bool relax)
{
    /*
    mDot_.storePrevIter();
    const volScalarField& mDot0 = mDot_.prevIter();

    volVectorField gradAlpha = fvc::grad(alpha_);

    const scalar a = sqrt(gv_*Rv_/g0_/R0_);
    const scalar b = (g0_ + 1.0)/4.0;

    scalar tmp;
    scalar res = 0;
    scalar maxMDot = 0;
    forAll(mDot_, cellI)
    {
        const scalar T = T_[cellI];
        const scalar Ma = max(min((T-Tv_)/(Th_ - Tv_), 1.0), 0.0);
        const scalar m = sqrt(0.5*gv_)*Ma;

        tmp = 0.5*m*(gv_ - 1.0)/(gv_ + 1.0);
        const scalar sqrtTByTs =
            sqrt(1.0 + mathematical::pi*sqr(tmp)) - sqrt(mathematical::pi)*tmp;

        tmp = a*sqrt(T/T0_)*sqrtTByTs*Ma;
        const scalar pByP1 = 1 + g0_*tmp*(b*tmp + sqrt(1.0 + sqr(b*tmp)));

        mDot_[cellI] = sqrt(2.0/Rv_/T)*m*pByP1/sqrtTByTs*mag(gradAlpha[cellI]);
        pRec_[cellI] = ((1.0 + 2.0*sqr(m))*pByP1 - 1.0)*gradAlpha[cellI];

        res = max(res, mag(mDot_[cellI] - mDot0[cellI]));
        maxMDot = max(maxMDot, mDot_[cellI]);
    }
    mDot_.correctBoundaryConditions();
    pRec_.correctBoundaryConditions();

    res /= (maxMDot + q_);
    return returnReduce(res, maxOp<scalar>());
    */
    return gasDynamic::correct(relax);
}


// ************************************************************************* //
