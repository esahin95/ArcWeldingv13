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
    equationOfState rhoConst;
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
        rho         8000.0;
    }

    thermodynamics
    {
        Cp          520;
        hf          0;
    }

    transport
    {
        mu          0.004; // nu = 5e-7
        kappa       10.0;
    }
}

solidification
{
    type    none;
}

evaporation
{
    type    Knight;
    L       7.45e6;
    Tv      3068.0;
    Th      4210.0;
    Mv      55.8;
    p0      1e5;
    T0      300;
    M0      40.0;
    g0      1.6667;
    gv      1.6667;
    relax   0.5;
    q       1e-3;
}


// ************************************************************************* //
