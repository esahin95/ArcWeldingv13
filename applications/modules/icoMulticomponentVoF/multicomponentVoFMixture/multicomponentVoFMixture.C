/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2025 OpenFOAM Foundation
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

#include "multicomponentVoFMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multicomponentVoFMixture, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multicomponentVoFMixture::multicomponentVoFMixture
(
    const fvMesh& mesh
)
:
    compressibleMultiphaseVoFMixture(mesh),

    metals_(phases().size(), false),

    rhoCp_
    (
        IOobject
        (
            "rhoCp",
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar("rhoCp", dimEnergy/dimVolume/dimTemperature, 0)
    )
{
    wordList metals(lookup("metals"));

    forAll(metals, phasei)
    {
        if (phases().found(metals[phasei]))
        {
            metals_[phasei] = true;
        }
    }

    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::multicomponentVoFMixture::correctThermo()
{
    compressibleMultiphaseVoFMixture::correctThermo();
}


void Foam::multicomponentVoFMixture::correct()
{
    compressibleMultiphaseVoFMixture::correct();

    rhoCp_ *= 0.0;
    forAll(phases(), phasei)
    {
        rhoCp_ += phases()[phasei]
                 *phases()[phasei].thermo().rho()
                 *phases()[phasei].thermo().Cp();
    }
}


Foam::tmp<Foam::volScalarField> Foam::multicomponentVoFMixture::kappaEff
(
    const volScalarField& nut
)
{
    tmp<volScalarField> tkappaEff
    (
        phases()[0]
       *(
            phases()[0].thermo().kappa()
          + phases()[0].thermo().rho()*phases()[0].thermo().Cp()*nut
        )
    );

    for (label phasei=1; phasei<phases().size(); phasei++)
    {
        tkappaEff.ref() +=
            phases()[phasei]
           *(
               phases()[phasei].thermo().kappa()
             + phases()[phasei].thermo().rho()*phases()[phasei].thermo().Cp()*nut
            );
    }

    return tkappaEff;
}


// ************************************************************************* //
