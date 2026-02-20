#include "thermoModel.H"
#include "fvMesh.H"

namespace Foam
{
    defineTypeNameAndDebug(thermoModel, 0);
    defineRunTimeSelectionTable(thermoModel, dictionary);
}

Foam::thermoModel::thermoModel(const fvMesh& mesh, const word& group)
:
    physicalProperties(mesh, group)
{}

bool Foam::thermoModel::read()
{
    return physicalProperties::read();
}