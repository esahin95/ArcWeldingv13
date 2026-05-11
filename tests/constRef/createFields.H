Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField newp
(
    IOobject
    (
        "p",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(p.dimensions(), 1e+3)
);

const volScalarField& pRef = p;

forAll(p, cellI)
{
    p[cellI] = Zero;
}
p.correctBoundaryConditions();

newp = 5 * p;




p.write();
newp.write();

Info<<"p address at: "<<(&p)<<endl;
Info<<p<<endl;

Info<<"newp address at: "<<(&p)<<endl;
Info<<newp<<endl;

Info<<"pRef address at: "<<(&pRef)<<endl;
Info<<pRef<<endl;

std::cout<<(&p)<<std::endl;
std::cout<<(&newp)<<std::endl;
std::cout<<(&pRef)<<std::endl;