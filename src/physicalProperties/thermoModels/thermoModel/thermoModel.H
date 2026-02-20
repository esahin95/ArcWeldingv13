#ifndef thermoModel_H
#define thermoModel_H

#include "physicalProperties.H"

#include "primitiveFieldsFwd.H"
#include "volFieldsFwd.H"

namespace Foam
{

class fvMesh;

class thermoModel
:
    public physicalProperties
{
    
public:
    
    TypeName("thermoModel");

    declareRunTimeSelectionTable
    (
        autoPtr,
        thermoModel,
        dictionary,
        (
            const fvMesh& mesh,
            const word& group
        ),
        (mesh, group)
    );

    thermoModel(const fvMesh& mesh, const word& group);

    thermoModel(const thermoModel&) = delete;

    static autoPtr<thermoModel> New
    (
        const fvMesh& mesh,
        const word& group = word::null
    );

    virtual ~thermoModel()
    {}

    virtual bool read() = 0;

    virtual tmp<volScalarField> Cp() const = 0;

    virtual tmp<scalarField> Cp(const label patchi) const = 0;

    virtual tmp<volScalarField> kappa() const = 0;

    virtual tmp<scalarField> kappa(const label patchi) const = 0;

    virtual tmp<volScalarField> beta() const = 0;

    virtual tmp<scalarField> beta(const label patchi) const = 0;

    virtual void correct() = 0;

    void operator=(const thermoModel&) = delete;
};


}

#endif