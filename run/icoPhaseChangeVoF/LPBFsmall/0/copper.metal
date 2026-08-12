/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  13
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      binary;
    class       volScalarField;
    location    "0";
    object      copper.metal;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [];

internalField   uniform 0;

boundaryField
{
    ymax
    {
        type            zeroGradient;
    }
    xmin
    {
        type            zeroGradient;
    }
    xmax
    {
        type            zeroGradient;
    }
    ymin
    {
        type            zeroGradient;
    }
    zmin
    {
        type            zeroGradient;
    }
    zmax
    {
        type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;
    }
}


// ************************************************************************* //
