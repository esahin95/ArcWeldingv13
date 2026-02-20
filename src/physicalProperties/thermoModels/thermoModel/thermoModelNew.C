#include "thermoModel.H"
#include "fvMesh.H"

Foam::autoPtr<Foam::thermoModel> Foam::thermoModel::New
(
    const fvMesh& mesh,
    const word& group
)
{
    const word modelType
    (
        IOdictionary
        (
            thermoModel::findModelDict(mesh, group)
        ).lookupBackwardsCompatible
        (
            {
                "thermoModel"
            }
        )
    );

    Info<<"Selecting thermo model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter = 
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown thermo model " << modelType << nl << nl
            << "Valid thermo models are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<thermoModel>(cstrIter()(mesh, group));
}