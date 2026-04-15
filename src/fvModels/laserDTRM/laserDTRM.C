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

#include "laserDTRM.H"
#include "compressibleTwoPhaseVoFMixture.H"
#include "fvmSup.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(laserDTRM, 0);

        addToRunTimeSelectionTable
        (
            fvModel,
            laserDTRM,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::laserDTRM::laserDTRM
(
    const word& sourceName,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(sourceName, modelType, mesh, dict),

    phaseName_(dict.lookup("phase")),

    thermo_
    (
        mesh.lookupObject<fluidThermo>
        (
            IOobject::groupName(physicalProperties::typeName, phaseName_)
        )
    ),

    Qtot_(dict.lookup<scalar>("Q")),

    nRays_(dict.lookup<scalar>("nRays")),

    pos_
    (
        Function1<vector>::New
        (
            "position",
            dimless,
            dimless,
            dict
        )
    ),

    rad_(dict.lookup<scalar>("radius")),

    normal_(normalised(dict.lookup<vector>("normal"))),

    radial1_(normalised(perpendicular(normal_))),

    radial2_(normalised(normal_ ^ radial1_)),

    powerDist_
    (
        Function1<scalar>::New
        (
            "powerDist",
            dimless,
            dimless,
            dict
        )
    ),

    relax_(dict.lookupOrDefault<scalar>("relax", 1.0)),

    a_(dict.lookupOrDefault<scalar>("absorption", 1.0)),

    Q_
    (
        IOobject
        (
            "Q",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimPower, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),

    curTimeIndex_(-1)
{
    Q_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::laserDTRM::addSupFields() const
{
    return wordList({thermo_.he().name()});
}


void Foam::fv::laserDTRM::correct()
{
    if (curTimeIndex_ == mesh().time().timeIndex())
    {
        return;
    }

    Q_ = Zero;

    // Populate the cloud
    //IDLList

    // Tracking data


    // Ray tracing


    // Finalize computation

    curTimeIndex_ = mesh().time().timeIndex();
}


void Foam::fv::laserDTRM::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& he,
    fvMatrix<scalar>& eqn
) const
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name() << endl;
    }

    if (&he == &thermo_.he())
    {
        Info<< "This part is called with: &he == &thermo_.he()" << endl;
    }
    else
    {
        FatalErrorInFunction
            << "Support for field " << he.name() << " is not implemented"
            << exit(FatalError);
    }
}


void Foam::fv::laserDTRM::topoChange(const polyTopoChangeMap& map)
{}


void Foam::fv::laserDTRM::mapMesh(const polyMeshMap& map)
{}


void Foam::fv::laserDTRM::distribute(const polyDistributionMap& map)
{}


bool Foam::fv::laserDTRM::movePoints()
{
    return true;
}


// ************************************************************************* //
