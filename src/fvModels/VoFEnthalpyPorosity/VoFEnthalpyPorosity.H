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

Class
    Foam::fv::VoFEnthalpyPorosity

Description
    Solidification and melting model for VoF simulations.

    The presence of the solid phase in the flow field is incorporated into the
    model as a momentum porosity contribution; the energy associated with the
    phase change is added as an enthalpy contribution.  The solid fraction as a
    function of temperature \c alphaSolidT is specified as a Foam::Function1.

    The model writes the field \c alpha[01].solid which can be visualised to to
    show the solid distribution.

Usage
    Example usage:
    \verbatim
    VoFEnthalpyPorosity1
    {
        type            VoFEnthalpyPorosity;

        alphaSolidT     table
        (
            (330 1)
            (335 0)
        );

        L               334000;
    }
    \endverbatim

    Where:
    \table
        Property     | Description             | Required    | Default value
        alphaSolidT  | Solid fraction as function of temperature | yes |
        L            | Latent heat of fusion [J/kg] | yes    |
        relax        | Relaxation coefficient [0-1] | no     | 0.9
        Cu           | Model coefficient       | no          | 100000
        q            | Model coefficient       | no          | 0.001
    \endtable

See also
    Foam::fv::solidificationMeltingSource
    Foam::Function1

SourceFiles
    VoFEnthalpyPorosity.C
    VoFEnthalpyPorosityIO.C

\*---------------------------------------------------------------------------*/

#ifndef VoFEnthalpyPorosity_H
#define VoFEnthalpyPorosity_H

#include "fvModel.H"
#include "fvMesh.H"
#include "volFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class incompressibleTwoPhaseVoFMixture;

namespace fv
{

/*---------------------------------------------------------------------------*\
               Class VoFEnthalpyPorosity Declaration
\*---------------------------------------------------------------------------*/

class VoFEnthalpyPorosity
:
    public fvModel
{
    // Private Data

        //- Solid fraction as a function of temperature
        autoPtr<Function1<scalar>> alphaSolidT_;

        //- Latent heat of fusion [J/kg]
        dimensionedScalar L_;

        //- Phase fraction under-relaxation coefficient
        scalar relax_;

        //- Mushy region momentum sink coefficient [1/s]; default = 10^5
        scalar Cu_;

        //- Coefficient used in porosity calc - default = 0.001
        scalar q_;

        //- VoF mixture
        const incompressibleTwoPhaseVoFMixture& mixture_;

        //- Solid phase fraction
        mutable volScalarField alphaSolid_;


    // Private Member Functions

        //- Non-virtual read
        void readCoeffs(const dictionary& dict);

        //- Return the name of the solid phase fraction
        word alphaSolidName() const;

        //- Helper function to apply to the energy equation
        template<class RhoFieldType>
        void apply(const RhoFieldType& rho, fvMatrix<scalar>& eqn) const;


public:

    //- Runtime type information
    TypeName("VoFEnthalpyPorosity");


    // Constructors

        //- Construct from explicit source name and mesh
        VoFEnthalpyPorosity
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        VoFEnthalpyPorosity
        (
            const VoFEnthalpyPorosity&
        ) = delete;


    // Member Functions

        // Checks

            //- Return the list of fields for which the option adds source term
            //  to the transport equation
            virtual wordList addSupFields() const;


        // Add explicit and implicit contributions to compressible equation

            //- Add explicit contribution to energy equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& Cp,
                const volScalarField& T,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add implicit contribution to mixture momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                const volVectorField& U,
                fvMatrix<vector>& eqn
            ) const;


        // Mesh motion

            //- Update topology using the given map
            virtual void topoChange(const polyTopoChangeMap&);

            //- Update from another mesh using the given map
            virtual void mapMesh(const polyMeshMap&);

            //- Redistribute or update using the given distribution map
            virtual void distribute(const polyDistributionMap&);

            //- Update for mesh motion
            virtual bool movePoints();


        //- Correct the model
        virtual void correct();


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const VoFEnthalpyPorosity&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
