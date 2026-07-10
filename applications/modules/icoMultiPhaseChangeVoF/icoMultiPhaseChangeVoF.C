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

#include "icoMultiPhaseChangeVoF.H"

#include "addToRunTimeSelectionTable.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(icoMultiPhaseChangeVoF, 0);
    addToRunTimeSelectionTable(solver, icoMultiPhaseChangeVoF, fvMesh);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::icoMultiPhaseChangeVoF::icoMultiPhaseChangeVoF
(
    fvMesh& mesh
)
:
    icoMulticomponentVoF(mesh),
    solModels_(phases.size()),
    evaModels_(phases.size())
{
    forAll(phases, phasei)
    {
        solModels_.set
        (
            phasei,
            solidificationModel::New(mesh, phases[phasei].name())
        );

        evaModels_.set
        (
            phasei,
            evaporationModel::New(mesh, phases[phasei].name())
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::icoMultiPhaseChangeVoF::~icoMultiPhaseChangeVoF()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //



// ************************************************************************* //
