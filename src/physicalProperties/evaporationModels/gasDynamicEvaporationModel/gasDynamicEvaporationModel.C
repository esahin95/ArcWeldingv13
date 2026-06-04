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

    R_
    (
        "twoPiR",
        physicoChemical::R/
        dimensionedScalar
        (
            "M",
            dimMass/dimMoles,
            dict_.lookup<scalar>("M")*1e-3 // expects g/mol
        )
    ),

    pSat_
    (
        IOobject
        (
            IOobject::groupName("pSat", group),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPressure, 0.0)
    )
{
    correct(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::evaporationModels::gasDynamic::correct(const bool relax)
{
    // Saturated vapor pressure
    Info<<"Update pSat"<<endl;
    pSat_ = p0_*exp((Lv_/R_/T0_)*(1.0 - (T0_/T_)));
    forAll(pSat_, cellI)
    {
        if (pSat_[cellI] <= p0_.value())
        {
            pSat_[cellI] *= 0.0;
        }
    }

    // Mass transfer rate
    Info<<"Update mDot"<<endl;
    mDot_ =
        0.816*pSat_/sqrt((mathematical::twoPi*R_)*T_)*mag(fvc::grad(alpha_));

    return 0.0;
}


void Foam::evaporationModels::gasDynamic::hSource(fvScalarMatrix& TEqn) const
{
    // Add latent heat of evaporation
    TEqn += Lv_*mDot_;
}


void Foam::evaporationModels::gasDynamic::USource(fvVectorMatrix& UEqn) const
{
    // Add recoil pressure term
    UEqn -= 0.55*pSat_*fvc::grad(alpha_);
}

// ************************************************************************* //
