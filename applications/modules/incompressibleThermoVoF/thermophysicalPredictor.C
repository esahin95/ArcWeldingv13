
#include "incompressibleThermoVoF.H"
#include "fvmDiv.H"
#include "fvcDdt.H"
#include "fvmLaplacian.H"

void Foam::solvers::incompressibleThermoVoF::thermophysicalPredictor()
{
    rhoCp == rho * mixture.Cp();
    rhoPhiCp == rhoPhi * fvc::interpolate(mixture.Cp());
    
    //- update temperature field
    T.storePrevIter();
    while (pimple.correctNonOrthogonal())
    {                
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp, T)
            + fvm::div(rhoPhiCp, T)
            - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoPhiCp), T)
            - fvm::laplacian(mixture.kappa(), T)
            ==
            L_ * fvc::ddt(rho, alphaSolid)
            + fvModels().source(rho, rhoCp/rho, T)
        );

        TEqn.relax();
        
        fvConstraints().constrain(TEqn);
        
        TEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            T.relax();
            Info<< "maximum change in T: " << gMax(mag(T.prevIter().primitiveField()-T.primitiveField())) << endl;
        }
    }
}