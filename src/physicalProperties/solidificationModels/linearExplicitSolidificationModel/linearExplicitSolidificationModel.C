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

    dict_(subDict(dictName_)),

    Lm_
    (
        "L",
        dimEnergy/dimMass,
        dict_.lookup<scalar>("L")
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

    Cu_
    (
        //"Cu",
        //dimless,
        dict_.lookupOrDefault<scalar>("Cu", 1e5)
    ),

    q_
    (
        //"q",
        //dimless,
        dict_.lookupOrDefault<scalar>("q", 0.001)
    ),

    relax_
    (
        //"relax",
        //dimless,
        dict_.lookupOrDefault<scalar>("relax", 1.0)
    ),

    alphaSolid_
    (
        IOobject
        (
            IOobject::groupName("alpha.solid", group),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0.0)
    )
{
    correct(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar
Foam::solidificationModels::linearExplicit::correct(const bool relax)
{
    tmp<volScalarField> tsfNew
    (
        volScalarField::New("sfNew", mesh_, dimless)
    );
    volScalarField& sfNew = tsfNew.ref();

    // Relaxation
    if (relax)
    {
        sfNew = max
            (
                min
                (
                    sf_ + relax_*thermo_.Cp()/Lm_*
                    (
                        Tliq_ - T_ - sf_*(Tliq_ - Tsol_)
                    ),
                    1.0
                ),
                0.0
            );
    }
    else
    {
        sfNew = max(min((Tliq_ - T_)/(Tliq_ - Tsol_), 1.0), 0.0);
    }

    // Residual
    const scalar res =
        gMax(mag(sf_.primitiveField() - sfNew.primitiveField()));

    // Update fields
    sf_ = sfNew;
    alphaSolid_ = alpha_*sf_;

    return res;
}


void Foam::solidificationModels::linearExplicit::addSup
(
    fvMatrix<scalar>& eqn
) const
{
    eqn -= Lm_*
        (
            fvc::ddt(alpha_, thermo_.rho(), sf_)
          + fvc::div(alphaRhoPhi_, sf_)
        );
}


void Foam::solidificationModels::linearExplicit::addSup
(
    fvMatrix<vector>& eqn
) const
{
    const scalarField& V = mesh_.V();
    scalarField& Sp = eqn.diag();
    forAll(Sp, cellI)
    {
        const scalar sf = alpha_[cellI]*sf_[cellI];
        Sp[cellI] += V[cellI]*Cu_*sqr(sf)/(pow3(1.0 - sf) + q_);
    }
}

// ************************************************************************* //
