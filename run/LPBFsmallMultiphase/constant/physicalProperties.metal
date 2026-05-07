/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  13
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "constant";
    object      physicalProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

thermoType
{
    type            heRhoThermo;
    mixture         pureMixture;
    transport       const;
    thermo          eConst;
    equationOfState rhoConst;//Boussinesq;
    specie          specie;
    energy          sensibleInternalEnergy;
}

mixture
{
    specie
    {
        molWeight   55.8; // Predominantly iron
    }

    equationOfState
    {
        rho0        8000.0;
        T0          300;
        beta        5e-6;
        rho         8000.0;
    }

    thermodynamics
    {
        Cv          520;
        hf          0; // L = 2.7e5
    }

    transport
    {
        mu          0.004; // nu = 5e-7
        kappa       10.0;
    }
}


// ************************************************************************* //
