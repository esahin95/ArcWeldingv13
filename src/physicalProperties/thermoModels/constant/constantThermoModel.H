#ifndef constantThermoModel_H
#define constantThermoModel_H

#include "thermoModel.H"
#include "volFields.H"

namespace Foam
{
namespace thermoModels
{


class constant
:
    public thermoModel
{    
    dimensionedScalar Cp0_;
    dimensionedScalar kappa0_;
    dimensionedScalar beta0_;

    volScalarField Cp_;
    volScalarField kappa_;
    volScalarField beta_;

public:

    TypeName("constant");

    constant(const fvMesh& mesh, const word& group);

    virtual ~constant()
    {}

    virtual tmp<volScalarField> Cp() const
    {
        return Cp_;
    }

    virtual tmp<scalarField> Cp(const label patchi) const
    {
        return Cp_.boundaryField()[patchi];
    }

    virtual tmp<volScalarField> kappa() const
    {
        return kappa_;
    }

    virtual tmp<scalarField> kappa(const label patchi) const
    {
        return kappa_.boundaryField()[patchi];
    }

    virtual tmp<volScalarField> beta() const
    {
        return beta_;
    }

    virtual tmp<scalarField> beta(const label patchi) const
    {
        return beta_.boundaryField()[patchi];
    }

    virtual void correct()
    {}

    virtual bool read();
};

}
}

#endif