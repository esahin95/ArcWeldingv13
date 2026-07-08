
// Mixture molecular mass
const volScalarField Wm(this->thermo().W());

// Sum of mole fraction by diffusion
volScalarField sumXbyD
(
    volScalarField::New
    (
        "sumXbyD",
        T.mesh(),
        dimless/dimKinematicViscosity/Wm.dimensions()
    )
);

forAll(Dm_, i)
{
    // Initialize to zero
    sumXbyD = Zero;

    forAll(Y, j)
    {
        if (j != i)
        {
            sumXbyD +=
                Y[j]
                /(
                    this->thermo().Wi(j)
                    *(
                        i < j
                        ? evaluate
                        (
                            DFuncs_[i][j],
                            dimKinematicViscosity,
                            p,
                            T
                        )
                        : evaluate
                        (
                            DFuncs_[j][i],
                            dimKinematicViscosity,
                            p,
                            T
                        )
                    )
                );
        }
    }

    Dm_.set
    (
        i,
        (1/Wm - Y[i]/this->thermo().Wi(i))
        /max(sumXbyD, dimensionedScalar(sumXbyD.dimensions(), small))
    );
}