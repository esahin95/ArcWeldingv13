#include "constantThermoModel.H"
#include "addToRunTimeSelectionTable.H"


namespace Foam
{
namespace thermoModels
{
    defineTypeNameAndDebug(constant, 0);

    addToRunTimeSelectionTable(thermoModel, constant, dictionary);

    addNamedToRunTimeSelectionTable
    (
        thermoModel,
        constant,
        dictionary,
        Simple
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermoModels::constant::constant
(
    const fvMesh& mesh, 
    const word& group
)
:
    thermoModel(mesh, group),

    Cp0_("Cp", dimensionSet(0, 2, -2, -1, 0), *this),

    kappa0_("kappa", dimensionSet(1, 1, -3, -1, 0), *this),

    beta0_("beta", dimensionSet(0, 0, 0, -1, 0), *this),

    Cp_
    (
        IOobject
        (
            IOobject::groupName("Cp", group),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        Cp0_
    ),

    kappa_
    (
        IOobject
        (
            IOobject::groupName("kappa", group),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        kappa0_
    ),
    
    beta_
    (
        IOobject
        (
            IOobject::groupName("beta", group),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        beta0_
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::thermoModels::constant::read()
{
    if (thermoModel::read())
    {
        Cp0_.read(*this);
        Cp_ = Cp0_;

        kappa0_.read(*this);
        kappa_ = kappa0_;

        beta0_.read(*this);
        beta_ = beta0_;

        return true;
    }
    else 
    {
        return false;
    }
}