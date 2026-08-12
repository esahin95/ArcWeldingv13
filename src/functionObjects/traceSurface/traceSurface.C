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

#include "traceSurface.H"
#include "fvMesh.H"
#include "interpolation.H"
#include "IOmanip.H"
#include "lineCellFace.H"
#include "Time.H"
#include "uniformDimensionedFields.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(traceSurface, 0);
    addToRunTimeSelectionTable(functionObject, traceSurface, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::traceSurface::writePositions()
{
    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);

    autoPtr<interpolation<scalar>>
        interpolator
        (
            interpolation<scalar>::New(interpolationScheme_, alpha)
        );

    const vector d0 = corners_[0];
    const vector d1 = corners_[1] - corners_[0];
    const vector d2 = corners_[2] - corners_[1];
    const scalar dx = mag(d1) / scalar(nx_);
    const scalar dy = mag(d2) / scalar(ny_);
    const scalar dhx = dx / scalar(ns_);
    const scalar dhy = dy / scalar(ns_);

    scalar x = 0.0;
    for (label i=0; i<nx_; i++)
    {
        scalar y = 0.0;
        for (label j=0; j<ny_; j++)
        {



            Info << x+0.5*dx << " " << y+0.5*dy << endl;
            y += dy;
        }

        x += dx;
    }


    if (Pstream::master())
    {
        writeTime(file());
    }



    if (Pstream::master())
    {
        file().endl();
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::traceSurface::writeFileHeader(const label i)
{
    writeCommented(file(), "Location");
    file().endl();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::traceSurface::traceSurface
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(mesh_, name),
    alphaName_(dict.lookup("alpha")),
    corners_(dict.lookup("corners")),
    nx_(dict.lookup<label>("nx")),
    ny_(dict.lookup<label>("ny")),
    ns_(dict.lookup<label>("ns")),
    interpolationScheme_(
        dict.lookupOrDefault<word>("interpolationScheme", "cellPoint")
    )
{
    read(dict);

    resetName(typeName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::traceSurface::~traceSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::wordList Foam::functionObjects::traceSurface::fields() const
{
    return wordList(alphaName_);
}


bool Foam::functionObjects::traceSurface::execute()
{
    return true;
}


bool Foam::functionObjects::traceSurface::end()
{
    return true;
}


bool Foam::functionObjects::traceSurface::write()
{
    logFiles::write();

    writePositions();

    return true;
}


// ************************************************************************* //
