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
#include "writeFile.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
void gatherAndFlatten(DynamicField<Type>& field)
{
    List<List<Type>> gatheredField(Pstream::nProcs());
    gatheredField[Pstream::myProcNo()] = field;
    Pstream::gatherList(gatheredField);

    field =
        ListListOps::combine<List<Type>>
        (
            gatheredField,
            accessOp<List<Type>>()
        );
}

}


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

        randomGenerator laserDTRM::rndGen_(261782, true);
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

    a_(dict.lookupOrDefault<scalar>("absorption", 1e+6)),

    reflectionModelPtr_
    (
        reflectionModel::New
        (
            dict.subDict("reflectionModel"),
            mesh
        )
    ),

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

    curTimeIndex_(-1),

    formatterPtr_(setWriter::New(dict.lookup("setFormat"), dict)),

    outputPath_
    (
        mesh().time().globalPath()
        /functionObjects::writeFile::outputPrefix
        /name()
    ),

    allPositions_(),
    allTracks_(),
    allPowers_(),
    writeIndex_(-1)
{
    Q_.oldTime();

    mkDir(outputPath_);
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

    Q_.storePrevIter();
    Q_ = Zero;

    const meshSearch& searchEngine = meshSearch::New(mesh());

    lagrangian::Cloud<DTRMParticle> cloud
    (
        mesh(),
        "DTRMCloud",
        IDLList<DTRMParticle>()
    );
    DTRMParticle::nParticles = 0; // Maybe rename to nTracks

    // Compute particle positions
    List<vector> positions
    (
        nRays_, pos_->value(mesh().time().value())
    );

    forAll(positions, trackIndex)
    {
        // Rejection sampling from arbitrary distribution
        vector delta;
        scalar deltaMag;
        do
        {
            delta =
                  radial1_ * rndGen_.scalarAB(-rad_, rad_)
                + radial2_ * rndGen_.scalarAB(-rad_, rad_);

            deltaMag = mag(delta);
        }
        while
        (
            deltaMag > rad_ ||
            rndGen_.scalar01() > powerDist_->value(deltaMag)
        );

        positions[trackIndex] += delta;
    }

    // Populate cloud
    label nLocateBoundaryHits = 0;
    forAll(positions, trackIndex)
    {
        const vector& position = positions[trackIndex];

        const label cellI = searchEngine.findCell(position);

        // Generate particle for a single processor only
        label candidate = cellI!=-1? Pstream::myProcNo() : -1;
        label owner = returnReduce(candidate, maxOp<label>());
        if (owner == Pstream::myProcNo())
        {
            cloud.addParticle(new DTRMParticle
                (
                    searchEngine,
                    position,
                    cellI,
                    nLocateBoundaryHits,
                    normal_,
                    Qtot_ / nRays_
                )
            );
        }
    }

    // Tracking data fields
    const volVectorField nHat = fvc::grad(alpha_);
    const volScalarField absorp = a_ * (1 - alpha_);

    // Construct tracking data
    interpolationCellPoint<scalar> alphaInterp(alpha_);
    interpolationCellPoint<scalar> absorpInterp(absorp);
    interpolationCellPoint<vector> nHatInterp(nHat);
    allPositions_.clear();
    allTracks_.clear();
    allPowers_.clear();
    DTRMParticle::trackingData td
    (
        cloud,
        alphaInterp,
        absorpInterp,
        nHatInterp,
        Q_,
        allPositions_,
        allTracks_,
        allPowers_,
        searchEngine,
        reflectionModelPtr_()
    );

    // Ray tracing
    cloud.move(cloud, td);

    // Finalize computation
    Q_.primitiveFieldRef() /= mesh().V().primitiveField();

    Q_.relax(relax_);

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

        eqn += Q_;
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

bool Foam::fv::laserDTRM::write(const bool write) const
 {
    if (Pstream::parRun())
    {
        gatherAndFlatten(allPositions_);
        gatherAndFlatten(allTracks_);
        gatherAndFlatten(allPowers_);
    }

    if (Pstream::master() && allPositions_.size())
    {
        DebugInfo<< "Writing out rays for time "
                 << mesh().time().name()
                 << " in directory "
                 << outputPath_
                 <<endl;

        const word writeIndex = Time::timeName(++writeIndex_);
        formatterPtr_->write
        (
            outputPath_,
            IOobject::groupName("traced", writeIndex),
            coordSet(allTracks_, word::null, allPositions_),
            "Power",
            allPowers_
        );
    }

    return true;
 }


// ************************************************************************* //
