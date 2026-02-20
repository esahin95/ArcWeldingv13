/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2017-2025 OpenFOAM Foundation
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

#include "VoFEnthalpyPorosity.H"
#include "fvcDdt.H"
#include "incompressibleTwoPhaseVoFMixture.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(VoFEnthalpyPorosity, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            VoFEnthalpyPorosity,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::VoFEnthalpyPorosity::readCoeffs(const dictionary& dict)
{
    alphaSolidT_.reset
    (
        Function1<scalar>::New
        (
            "alphaSolidT",
            dimTemperature,
            unitFraction,
            dict
        ).ptr()
    );
    L_ = dimensionedScalar("L", dimEnergy/dimMass, dict);
    relax_ = dict.lookupOrDefault<scalar>("relax", dimless, 0.9);
    Cu_ = dict.lookupOrDefault<scalar>("Cu", dimless/dimTime, 100000);
    q_ = dict.lookupOrDefault<scalar>("q", dimless, 0.001);
}


Foam::word Foam::fv::VoFEnthalpyPorosity::alphaSolidName() const
{    
    return IOobject::groupName(mixture_.alpha1().name(), "solid");
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::VoFEnthalpyPorosity::VoFEnthalpyPorosity
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),

    alphaSolidT_(),

    L_("L", dimEnergy/dimMass, NaN),

    relax_(NaN),

    Cu_(NaN),

    q_(NaN),

    mixture_
    (
        mesh().lookupObject<incompressibleTwoPhaseVoFMixture>
        (
            "phaseProperties"
        )
    ),

    alphaSolid_
    (
        IOobject
        (
            alphaSolidName(),
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    readCoeffs(coeffs(dict));

    const volScalarField& TVoF = mesh().lookupObject<volScalarField>("T");
    const volScalarField& alphaVoF = mixture_.alpha1();
    alphaSolid_.primitiveFieldRef() = alphaVoF * alphaSolidT_->value(TVoF.primitiveField());
    alphaSolid_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::VoFEnthalpyPorosity::addSupFields() const
{
    return wordList({"U", "T"});
}


void Foam::fv::VoFEnthalpyPorosity::addSup
(
    const volScalarField& rho,
    const volScalarField& Cp,
    const volScalarField& T,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }
    Info<< type() << ": applying source to " << eqn.psi().name() << endl;

    eqn += L_ * (fvc::ddt(rho, alphaSolid_));
}


void Foam::fv::VoFEnthalpyPorosity::addSup
(
    const volScalarField& rho,
    const volVectorField& U,
    fvMatrix<vector>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }
    Info<< type() << ": applying source to " << eqn.psi().name() << endl;

    eqn.diag() -= mesh().V() * rho * Cu_ * sqr(alphaSolid_) / (pow3(1 - alphaSolid_) + q_);
}


void Foam::fv::VoFEnthalpyPorosity::correct()
{
    if (debug)
    {
        Info<< type() << ": " << name()
            << " - updating solid phase fraction" << endl;
    }

    Info<< type() << ": " << name() << " - updating solid phase fraction" << endl;

    const volScalarField& alphaSolidOld = alphaSolid_.oldTime();

    const volScalarField& TVoF = mesh().lookupObject<volScalarField>("T");
    const volScalarField& alphaVoF = mixture_.alpha1();

    alphaSolid_.primitiveFieldRef() = relax_ * alphaVoF * alphaSolidT_->value(TVoF.primitiveField()) + (1-relax_) * alphaSolidOld.primitiveField() ;
    alphaSolid_.correctBoundaryConditions();

    Info<< "maximum change in solid fraction: " << gMax(mag(alphaSolid_.primitiveField() - alphaSolidOld.primitiveField())) << endl;
}


void Foam::fv::VoFEnthalpyPorosity::topoChange
(
    const polyTopoChangeMap& map
)
{}


void Foam::fv::VoFEnthalpyPorosity::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::VoFEnthalpyPorosity::distribute
(
    const polyDistributionMap& map
)
{}


bool Foam::fv::VoFEnthalpyPorosity::movePoints()
{
    return true;
}


bool Foam::fv::VoFEnthalpyPorosity::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs(coeffs(dict));
        return true;
    }
    else
    {
        return false;
    }

    return false;
}


// ************************************************************************* //
