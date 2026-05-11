#include "ErrorFunction1.H"
#include "mathematicalConstants.H"

namespace Foam
{
    namespace Function1s
    {
        addFunction1(ErrorFunction, scalar);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::ErrorFunction<Type>::ErrorFunction
(
    const word& name,
    const unitConversions& units,
    const dictionary& dict
)
:
    FieldFunction1<Type, ErrorFunction<Type>>(name),
    Tliq_(dict.lookup<scalar>("Tliq")),
    Tsol_(dict.lookup<scalar>("Tsol")),
    Tmid_(0.5 * (Tsol_ + Tliq_)),
    a_(4.0 / (max(1e-6, Tliq_ - Tsol_)))
{}

template<class Type>
Foam::Function1s::ErrorFunction<Type>::ErrorFunction(const ErrorFunction& errf)
:
    FieldFunction1<Type, ErrorFunction<Type>>(errf),
    Tliq_(errf.Tliq_),
    Tsol_(errf.Tsol_),
    Tmid_(errf.Tmid_),
    a_(errf.a_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1s::ErrorFunction<Type>::~ErrorFunction()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1s::ErrorFunction<Type>::value(const scalar x) const 
{
    Type y(0.5 * (1.0 - Foam::erf(a_ * (x - Tmid_))));

    return y;
}

template<class Type>
Type Foam::Function1s::ErrorFunction<Type>::integral
(
    const scalar x1, 
    const scalar x2
) const 
{
    // Not implemented
    Type y(Zero);

    return y;
}

template<class Type>
void Foam::Function1s::ErrorFunction<Type>::write
(
    Ostream& os, 
    const unitConversions& units
) const 
{
    writeEntry(os, "Tsol", Tsol_);
    writeEntry(os, "Tliq", Tliq_);
}