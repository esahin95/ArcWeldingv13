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

#include "gasDynamicEvaporationModel.H"
#include "addToRunTimeSelectionTable.H"

#include "fvcGrad.H"
#include "mathematicalConstants.H"
#include "physicoChemicalConstants.H"
#include "fvcVolumeIntegrate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace evaporationModels
{
    defineTypeNameAndDebug(gasDynamic, 0);

    addToRunTimeSelectionTable
    (
        evaporationModel,
        gasDynamic,
        dictionary
    );
}
}

using namespace Foam::constant;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::evaporationModels::gasDynamic::gasDynamic
(
    const fvMesh& mesh,
    const word& group
)
:
    evaporationModel(mesh, group),

    dict_(subDict(evaporationModel::typeName)),

    Lv_
    (
        "Lv",
        dimEnergy/dimMass,
        dict_.lookup<scalar>("Lv")
    ),

    p0_
    (
        "p0",
        dimPressure,
        dict_.lookup<scalar>("p0")
    ),

    T0_
    (
        "T0",
        dimTemperature,
        dict_.lookup<scalar>("T0")
    ),

    R0_
    (
        "R0",
        physicoChemical::R/
        dimensionedScalar
        (
            "M",
            dimMass/dimMoles,
            dict_.lookup<scalar>("M0")*1e-3 // expects g/mol
        )
    ),

    g0_(dict_.lookup<scalar>("g0")),

    Rv_
    (
        "Rv",
        physicoChemical::R/
        dimensionedScalar
        (
            "M",
            dimMass/dimMoles,
            dict_.lookup<scalar>("Mv")*1e-3 // expects g/mol
        )
    ),

    gv_(dict_.lookup<scalar>("gv")),

    Tv_
    (
        "Tv",
        dimTemperature,
        dict_.lookup<scalar>("Tv")
    ),

    coeffA_(0.0),

    coeffB_(0.0),

    pRec_
    (
        IOobject
        (
            IOobject::groupName("pRec", group),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPressure, 0.0)
    )
{
    // Jump conditions
    scalar m = sqrt(gv_/2.0);

    scalar tmp = 0.5*m*(gv_ - 1.0)/(gv_ + 1.0);

    scalar sqrtTByTs =
        sqrt(1 + mathematical::pi*sqr(tmp)) - sqrt(mathematical::pi)*tmp;

    scalar rhoByRhos =
        (
            (
                (sqr(m)+0.5)*exp(sqr(m))*erfc(m)
              - m/sqrt(mathematical::pi)
            )/sqrtTByTs
          + (
                1.0 - sqrt(mathematical::pi)*m*exp(sqr(m))*erfc(m)
            )/sqr(sqrtTByTs)/2.0
        );

    scalar pByPs = rhoByRhos*sqr(sqrtTByTs);

    // Precompute coefficients
    coeffA_ = 2*sqrt(mathematical::pi)*m*rhoByRhos*sqrtTByTs;
    coeffB_ = pByPs + sqr(coeffA_)/mathematical::twoPi/rhoByRhos;
    DebugInfo<< "Evaporation coefficient coeffA = " << coeffA_ <<endl;
    DebugInfo<< "Recoil pressure coefficient coeffB = " << coeffB_ <<endl;

    // Update mass transfer rate and recoil pressure
    correct(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::evaporationModels::gasDynamic::correct(const bool relax)
{
    mDot_.storePrevIter();

    // Saturated vapor pressure
    //DebugInfo<< "psat" <<endl;
    volScalarField pSat = p0_*exp(Lv_/Rv_/Tv_*(1.0 - Tv_/T_));

    // Mass transfer rate
    //DebugInfo<< "mDot" <<endl;
    //DebugInfo<< "Min. T = "
    //         << gMin(T_.primitiveField()) << endl;
    mDot_ = coeffA_*pSat/sqrt(mathematical::twoPi*Rv_*T_)*mag(fvc::grad(alpha_));

    // Recoil pressure
    //DebugInfo<< "prec" <<endl;
    pRec_ = coeffB_*pSat;

    return gMax(mag(mDot_.prevIter() - mDot_)().primitiveField());
}


void Foam::evaporationModels::gasDynamic::hSource(fvScalarMatrix& TEqn) const
{
    if (debug)
    {
        const dimensionedScalar hv = fvc::domainIntegrate(Lv_*mDot_);
        Info<< "Total evaporative enthalphy: " << hv << endl;
    }

    // Add latent heat of evaporation
    TEqn += Lv_*mDot_;
}


void Foam::evaporationModels::gasDynamic::USource(fvVectorMatrix& UEqn) const
{
    // Add recoil pressure term
    UEqn -= pRec_*fvc::grad(alpha_);
}

// ************************************************************************* //
