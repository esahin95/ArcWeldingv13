#ifndef ErrorFunction1_H
#define ErrorFunction1_H

#include "Function1.H"

namespace Foam
{
namespace Function1s
{

/*---------------------------------------------------------------------------*\
                         Class ErrorFunction Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ErrorFunction
:
    public FieldFunction1<Type, ErrorFunction<Type>>
{
    
    // Private data
    
        scalar Tliq_;
        scalar Tsol_;
        scalar Tmid_;
        scalar a_;

public:

    //- Runtime type information
    TypeName("errorFunction");

    // Constructors

        //- Construct from name and dictionary
        ErrorFunction
        (
            const word& name, 
            const unitConversions& units,
            const dictionary& dict
        );

        //- Copy constructor
        ErrorFunction(const ErrorFunction& errf);

    // Destructor
    
        virtual ~ErrorFunction();

    
    // Member functions

        //- Return error function value as a function of scalar x
        virtual Type value(const scalar x) const;

        //- Integrate between two scalar fields
        virtual Type integral(const scalar x1, const scalar x2) const;
    
        //- Write in dictionary format
        virtual void write(Ostream& os, const unitConversions& units) const;

    // Member operators

        //- Disallow default bitwise assignment
        void operator=(const ErrorFunction<Type>&) = delete;

};

} // End namespace Function1s
} // End namespace Foam

#endif