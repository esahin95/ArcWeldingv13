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

    Th_
    (
        "Th",
        dimTemperature,
        dict_.lookup<scalar>("Th")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::evaporationModels::Knight::correct(const bool relax)
{
    mDot_.storePrevIter();

    const scalar Th = Th_.value();
    const scalar Tv = Tv_.value();
    const scalar T0 = T0_.value();
    const scalar Rv = Rv_.value();
    const scalar R0 = R0_.value();

    const scalar a = sqrt(gv_*Rv/g0_/R0);
    const scalar b = (g0_ + 1.0)/4.0;

    scalar tmp;
    forAll(mDot_, cellI)
    {
        const scalar T = T_[cellI];
        const scalar Ma = max(min((T-Tv)/(Th-Tv), 1.0), 0.0);
        const scalar m = sqrt(0.5*gv_)*Ma;

        tmp = 0.5*m*(gv_-1.0)/(gv_+1.0);
        const scalar sqrtTByTs =
            sqrt(1.0 + mathematical::pi*sqr(tmp)) - sqrt(mathematical::pi)*tmp;

        /*
        const scalar rSqrtTByTs = 1.0/sqrtTByTs;
        const scalar rhoByRhos =
            (
                (
                    (sqr(m)+0.5)*exp(sqr(m))*erfc(m)
                - m/sqrt(mathematical::pi)
                )*rSqrtTByTs
            + (
                    1.0 - sqrt(mathematical::pi)*m*exp(sqr(m))*erfc(m)
                )*sqr(rSqrtTByTs)/2.0
            );
        */

        //const scalar pByPs = rhoByRhos*sqr(sqrtTByTs);

        tmp = a*sqrt(T/T0)*sqrtTByTs*Ma;
        const scalar pByP1 = 1 + g0_*tmp*(b*tmp + sqrt(1.0 + sqr(b*tmp)));

        mDot_[cellI] = sqrt(2.0/Rv/T)*m*pByP1/sqrtTByTs;
        pRec_[cellI] = (1.0 + 2.0*sqr(m))*pByP1 - 1.0;
    }
    mDot_.correctBoundaryConditions();
    pRec_.correctBoundaryConditions();

    return
        gMax(mag(mDot_.prevIter() - mDot_)().primitiveField())
      /(gMax(mDot_.primitiveField()) + q_);
}


// ************************************************************************* //
