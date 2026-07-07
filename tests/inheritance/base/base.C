#include "base.H"

#include "IOstreams.H"

Foam::base::base(Foam::label value)
:
    value(value)
{
    Info<< "Running Base constructor" << endl;
}

void Foam::base::method()
{
    Info<< "Running Base method " << second() << endl;
}