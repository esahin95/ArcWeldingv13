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

#include "linearImplicitSolidificationModel.H"
#include "addToRunTimeSelectionTable.H"

#include "fvcDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidificationModels
{
    defineTypeNameAndDebug(linearImplicit, 0);

    addToRunTimeSelectionTable
    (
        solidificationModel,
        linearImplicit,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidificationModels::linearImplicit::linearImplicit
(
    const fvMesh& mesh,
    const word& group
)
:
    linearExplicit(mesh, group),

    rhoCpApp_
    (
        IOobject::groupName("rhoCp", group),
        alpha_*thermo_.rho()*thermo_.Cp()
    ),

    rhoPhiCpApp_
    (
        IOobject::groupName("rhoPhiCp", group),
        alphaRhoPhi_*fvc::interpolate(thermo_.Cp())
    )
{
    rhoCpApp_.oldTime();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::scalar
Foam::solidificationModels::linearImplicit::correct(const bool relax)
{
    // Update heat capacities
    rhoCpApp_ =
        Lm_/(Tliq_-Tsol_)*alpha_*thermo_.rho()*pos(Tliq_-T_)*pos(T_-Tsol_);
    rhoPhiCpApp_ = alphaRhoPhi_*fvc::interpolate(rhoCpApp_/rho_);

    // Update solid fraction
    return linearExplicit::correct(relax);
}


void Foam::solidificationModels::linearImplicit::addSup
(
    fvMatrix<scalar>& eqn
) const
{
    eqn +=
    (
        fvm::ddt(rhoCpApp_, T_) + fvm::div(rhoPhiCpApp_, T_)
      - fvm::Sp(fvc::ddt(rhoCpApp_) + fvc::div(rhoPhiCpApp_), T_)
    );
}


// ************************************************************************* //
