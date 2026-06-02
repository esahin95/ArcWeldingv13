/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2024 OpenFOAM Foundation
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

#include "greyBodyRadiation.H"
#include "addToRunTimeSelectionTable.H"

#include "physicoChemicalConstants.H"
#include "fvcGrad.H"
#include "fvmSup.H"
#include "fvcVolumeIntegrate.H"


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(greyBodyRadiation, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            greyBodyRadiation,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::greyBodyRadiation::greyBodyRadiation
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),

    alpha_
    (
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", dict.lookup<word>("phase"))
        )
    ),

    eps_
    (
        dict.lookup<scalar>("eps") * constant::physicoChemical::sigma
    ),

    T0_
    (
        "T0",
        dimTemperature,
        dict.lookup<scalar>("T0")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::greyBodyRadiation::addSupFields() const
{
    return wordList(1, "T");
}


void Foam::fv::greyBodyRadiation::correct()
{}


void Foam::fv::greyBodyRadiation::addSup
(
    const volScalarField& rho,
    const volScalarField& T,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    fvScalarMatrix radEqn
    (
        fvm::Sp(4*eps_*pow3(T)*mag(fvc::grad(alpha_)), T)
        - eps_*(pow4(T0_) + 3*pow4(T))*mag(fvc::grad(alpha_))
    );

    if (debug)
    {
        const dimensionedScalar Qtot = fvc::domainIntegrate(radEqn&T);
        Info<< "Radiation in fvModel to source: " << Qtot << endl;
    }

    eqn += radEqn;
}


void Foam::fv::greyBodyRadiation::topoChange(const polyTopoChangeMap& map)
{}


void Foam::fv::greyBodyRadiation::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::greyBodyRadiation::distribute(const polyDistributionMap& map)
{}


bool Foam::fv::greyBodyRadiation::movePoints()
{
    return true;
}


// ************************************************************************* //
