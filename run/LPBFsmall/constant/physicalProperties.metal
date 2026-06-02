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
    thermo          hConst;
    equationOfState rhoConst;//Boussinesq;
    specie          specie;
    energy          sensibleEnthalpy;
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
        Cp          520;
        hf          0; // L = 2.7e5
    }

    transport
    {
        mu          0.004; // nu = 5e-7
        kappa       10.0;
    }
}

solidificationModel
{
    type    linearExplicit;
    Lm      2.7e5;
    Tsol    1658.0;
    Tliq    1723.0;
    Cu      1e10;
    q       1e-3;
    relax   0.8;
}


// ************************************************************************* //
