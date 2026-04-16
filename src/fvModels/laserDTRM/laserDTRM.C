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
#include "Cloud.H"
#include "DTRMParticle.H"
#include "fvcVolumeIntegrate.H"

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

    alpha_
    (
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", phaseName_)
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
        dimensionedScalar(dimPower/dimVolume, 0.0),
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
    const meshSearch& searchEngine = meshSearch::New(mesh());

    lagrangian::Cloud<DTRMParticle> cloud
    (
        mesh(),
        "DTRMCloud",
        IDLList<DTRMParticle>()
    );

    List<vector> positions(4, pos_->value(0.0));
    positions[1].x() += 0.01;
    positions[2].x() -= 0.01;
    positions[3].x() += 0.02;

    label nLocateBoundaryHits = 0;
    forAll(positions, trackIndex)
    {
        const vector& position = positions[trackIndex];

        const label cellI = searchEngine.findCell(position);

        label candidate = -1;
        if (cellI != -1)
        {
            candidate = Pstream::myProcNo();
        }
        label owner = returnReduce(candidate, maxOp<label>());
        DebugInfo<< owner <<endl;

        if (owner == Pstream::myProcNo())
        {
            cloud.addParticle(new DTRMParticle
                (
                    searchEngine,
                    position,
                    cellI,
                    nLocateBoundaryHits,
                    normal_,
                    Qtot_,
                    a_,
                    trackIndex
                )
            );
        }

        if (owner == -1)
        {
            DebugInfo
                <<"Cannot find owner cell for position = "
                <<position<<endl;
        }
    }

    /*
    forAll(positions, trackIndex)
    {
        const vector& position = positions[trackIndex];

        const label cellI = searchEngine.findCell(position);

        if (cellI!=-1)
        {
            if
            (
                true//returnReduce(myProcNo, minOp<label>()) == Pstream::myProcNo()
            )
            {
                cloud.addParticle(new DTRMParticle
                    (
                        searchEngine,
                        position,
                        cellI,
                        nLocateBoundaryHits,
                        normal_,
                        Qtot_,
                        a_,
                        1
                    )
                );
            }
        }

        if (returnReduce(cellI, maxOp<label>()) == -1)
        {
            DebugInfo
                <<"Cannot find owner cell for position = "
                <<position<<endl;
        }

        label inProc = cellI == -1? 0 : 1;
        DebugInfo
            << "Particle added "
            << returnReduce(inProc, sumOp<label>())
            << " times"
            <<endl;
    }
    */

    // Tracking data
    interpolationCellPoint<scalar> alphaInterp(alpha_);
    DTRMParticle::trackingData td
    (
        cloud,
        alphaInterp,
        Q_
    );

    DebugInfo
        << "Cloud size at start: "
        << returnReduce(cloud.size(), sumOp<label>())
        << endl;

    // Ray tracing
    cloud.move(cloud, td);

    DebugInfo
        << "Cloud size at end: "
        << returnReduce(cloud.size(), sumOp<label>())
        << endl;
    forAllConstIter(lagrangian::Cloud<DTRMParticle>, cloud, iter)
    {
        const DTRMParticle& p = iter();

        DebugInfo<<p.position(mesh())<<endl;
    }

    // Finalize computation
    const scalarField& V = mesh().V();
    forAll(Q_, cellI)
    {
        Q_[cellI] /= V[cellI];
    }

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

        if (debug)
        {
            const dimensionedScalar Qtot = fvc::domainIntegrate(Q_);
            Info<< "Adding laser deposition to source: " << Qtot << endl;
        }

        eqn -= Q_;
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
