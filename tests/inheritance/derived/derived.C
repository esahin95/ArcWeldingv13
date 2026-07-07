#include "derived.H"

#include "IOstreams.H"

Foam::derived::derived(Foam::label value)
:
    base(value)
{
    Info<< "Running Derived constructor" << endl;
}

void Foam::derived::method()
{
    base::method();
    Info<< "Running Derived method " << second() << endl;
}