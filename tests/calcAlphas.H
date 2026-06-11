Foam::multiphaseMixtureThermo::solveAlphas
(
    volScalarField *massdotterm
)
{
    const word alphaScheme("div(phi,alpha)");
    const word alpharScheme("div(phirb,alpha)");

    volScalarField condensate
    (
        IOobject
        (
            "condensate",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("condensate", dimensionSet(0, 0, 0, 0, 0), 1.0)
    );

    for (phaseModel& alpha : phases_)
    {
        if (alpha.isGaseous())
        {
            condensate -= alpha;
        }
    }

    int phasei = 0;
    PtrList<surfaceScalarField> alphaPhiCorrs(phases_.size());
    for (phaseModel& alpha : phases_)
    {
        alphaPhiCorrs.set
        (
            phasei,
            new surfaceScalarField
            (
                phi_.name() + alpha.name(),
                fvc::flux
                (
                    phi_,
                    alpha,
                    alphaScheme
                )
            )
        );

        surfaceScalarField& alphaPhiCorr = alphaPhiCorrs[phasei];
        for (phaseModel& alpha2 : phases_)
        {
            if (&alpha2 == &alpha) continue;

            surfaceScalarField phic
            (
                (mag(phi_))/mesh_.magSf()
            );

            surfaceScalarField phir(0.0*phic*nHatf(condensate, (1.0-condensate)));//JUST TO INITIALISE TO ZERO - must be a neater way

            surfaceScalarField phirD(0.0*mag(phi_/mesh_.magSf())*nHatf(condensate, (1.0-condensate)));//setto zero
            phirD*=0.0;

            scalarCoeffSymmCTable::const_iterator cAlpha
            (
                cAlphas_.find(interfacePair(alpha, alpha2))
            );

            if (cAlpha != cAlphas_.end())
            {
                surfaceScalarField phic2
                (
                    (mag(phi_) + mag(phir))/mesh_.magSf()
                );

                phir +=(min(cAlpha()*phic2, max(phic2))*(nHatf(alpha, alpha2)));
            }

            scalarCoeffSymmDTable::const_iterator dAlpha
            (
                dAlphas_.find(interfacePair(alpha, alpha2))
            );

            if (dAlpha != dAlphas_.end())
            {
                volScalarField coefffield
                (
                    IOobject
                    (
                        "coefffield",
                        mesh_.time().timeName(),
                        mesh_
                    ),
                    mesh_,
                    dimensionedScalar("valdiff",dimdiff_,dAlpha())
                );

                phirD-= fvc::interpolate(coefffield)*mesh_.magSf()*((fvc::interpolate(alpha2)*fvc::snGrad(alpha))-(fvc::interpolate(alpha)*fvc::snGrad(alpha2)));;
            }

            alphaPhiCorr += phirD + fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha,
                alpharScheme
            );
        }

        MULES::limit
        (
            1.0/mesh_.time().deltaT().value(),
            geometricOneField(),
            alpha,
            phi_,
            alphaPhiCorr,
            zeroField(),
            zeroField(),
            oneField(),
            zeroField(),
            true
        );

        ++phasei;
    }

    MULES::limitSum(alphaPhiCorrs);

    rhoPhi_ = dimensionedScalar(dimensionSet(1, 0, -1, 0, 0), Zero);

    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar(dimless, Zero)
    );

    volScalarField divU(fvc::div(fvc::absolute(phi_, U_)));

    phasei = 0;
    for (phaseModel& alpha : phases_)
    {
        surfaceScalarField& alphaPhi = alphaPhiCorrs[phasei];
        alphaPhi += upwind<scalar>(mesh_, phi_).flux(alpha);

        volScalarField::Internal Sp
        (
            IOobject
            (
                "Sp",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(alpha.dgdt().dimensions(), Zero)
        );

        volScalarField::Internal Su
        (
            IOobject
            (
                "Su",
                mesh_.time().timeName(),
                mesh_
            ),
            divU*min(alpha, 1.0) + Sus[phasei] - PCR*(min(alpha, 1.0))
        );

        {
            const volScalarField& rhoK = alpha.thermo().rho();
            const scalarField& rhoNew = rhoK.primitiveField();
            const scalarField& rhoOld = rhoK.oldTime().primitiveField();
            const scalar dt = mesh_.time().deltaTValue();

            forAll(Sp, celli)
            {
                const scalar drhodt =
                    (rhoNew[celli] - rhoOld[celli]) / dt;
                const scalar rhoSafe =
                    Foam::max(rhoNew[celli], scalar(1e-6));

                Sp[celli] -= drhodt / rhoSafe;
            }
        }

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha,
            alphaPhi,
            Sp,
            Su
        );

        rhoPhi_ += fvc::interpolate(alpha.thermo().rho())*alphaPhi;

        Info<< alpha.name() << " volume fraction, min, max = "
            << alpha.weightedAverage(mesh_.V()).value()
            << ' ' << min(alpha).value()
            << ' ' << max(alpha).value()
            << endl;

        sumAlpha += alpha;

        ++phasei;
    }

    Info<< "Phase-sum volume fraction, min, max = "
        << sumAlpha.weightedAverage(mesh_.V()).value()
        << ' ' << min(sumAlpha).value()
        << ' ' << max(sumAlpha).value()
        << endl;
}

/* Deleted finite volume fields
    - alphasum
    - evaprate
    - condrate
    - alphagen
    - rhogen
    - massgen
    - divU

volScalarField alphasum
(
    IOobject
    (
        "alphasum",
        mesh_.time().timeName(),
        mesh_
    ),
    mesh_,
    dimensionedScalar("alphasum", dimensionSet(0, 0, 0, 0, 0), 0.0)
);

volScalarField evaprate
(
    IOobject
    (
        "evaprate",
        mesh_.time().timeName(),
        mesh_
    ),
    mesh_,
    dimensionedScalar("evaprate", dimensionSet(0, 0, -1, 0, 0), 0)
);

volScalarField condrate
(
    IOobject
    (
        "condrate",
        mesh_.time().timeName(),
        mesh_
    ),
    mesh_,
    dimensionedScalar("condrate", dimensionSet(0, 0, -1, 0, 0), 0)
);

volScalarField alphagen
(
    IOobject
    (
        "alphagen",
        mesh_.time().timeName(),
        mesh_
    ),
    mesh_,
    dimensionedScalar("alphagen", dimensionSet(0, 0, -1, 0, 0), 0)
);

volScalarField rhogen
(
    IOobject
    (
        "rhogen",
        mesh_.time().timeName(),
        mesh_
    ),
    mesh_,
    dimensionedScalar("rhogen", dimensionSet(1, -3, 0, 0, 0), 0)
);

volScalarField massgen
(
    IOobject
    (
        "massgen",
        mesh_.time().timeName(),
        mesh_
    ),
    mesh_,
    dimensionedScalar("massgen", dimensionSet(1, -3, -1, 1, 0), 0)
);

volScalarField divU
(
    IOobject
    (
        "divU",
        mesh_.time().timeName(),
        mesh_
    ),
    mesh_,
    dimensionedScalar("divU", dimensionSet(0, 0, -1, 0, 0), 0)
);

*/