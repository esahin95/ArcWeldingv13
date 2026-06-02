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

#include "linearExplicitSolidificationModel.H"
#include "addToRunTimeSelectionTable.H"

#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidificationModels
{
    defineTypeNameAndDebug(linearExplicit, 0);

    addToRunTimeSelectionTable
    (
        solidificationModel,
        linearExplicit,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidificationModels::linearExplicit::linearExplicit
(
    const fvMesh& mesh,
    const word& group
)
:
    solidificationModel(mesh, group),

    Lm_
    (
        "Lm",
        dimEnergy/dimMass,
        dict_.lookup<scalar>("Lm")
    ),

    Tliq_
    (
        "Tliq",
        dimTemperature,
        dict_.lookup<scalar>("Tliq")
    ),

    Tsol_
    (
        "Tsol",
        dimTemperature,
        dict_.lookup<scalar>("Tsol")
    ),

    Cu_(dict_.lookupOrDefault<scalar>("Cu", 1e5)),

    q_(dict_.lookupOrDefault<scalar>("q", 0.001)),

    relax_(dict_.lookupOrDefault<scalar>("relax", 1.0))
{
    correct(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::solidificationModels::linearExplicit::correct(const bool relax)
{
    sf_.storePrevIter();

    if (relax)
    {
        const volScalarField& Cp = thermo_.Cp();

        sf_ += (relax_/Lm_)*Cp*(Tliq_ - T_ - sf_*(Tliq_ - Tsol_));
        sf_ = max(min(sf_, 1.0), 0.0);
    }
    else
    {
        sf_ = max(min((Tliq_ - T_)/(Tliq_ - Tsol_), 1.0), 0.0);
    }

    //latentHeat_ = alpha_*thermo_.rho()*Lm_*(1-sf_)/rho_;
    return gMax(mag(sf_.prevIter() - sf_)().primitiveField());
}


void Foam::solidificationModels::linearExplicit::hSource
(
    fvScalarMatrix& TEqn
) const
{
    TEqn -= Lm_*
        (
            fvc::ddt(alpha_, thermo_.rho(), sf_)
          + fvc::div(alphaRhoPhi_, sf_)
        );
}


void Foam::solidificationModels::linearExplicit::USource
(
    fvVectorMatrix& UEqn
) const
{
    const scalarField& V = mesh_.V();

    scalarField& Sp = UEqn.diag();
    forAll(Sp, cellI)
    {
        const scalar Vc = V[cellI];
        const scalar sf = alpha_[cellI]*sf_[cellI];
        Sp[cellI] += Vc*Cu_*sqr(sf)/(pow3(1.0 - sf) + q_);
    }
}

// ************************************************************************* //
