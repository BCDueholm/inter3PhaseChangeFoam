    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + turbulence->divDevTau(rho, U)
     ==
        fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
      + fvModels.source(rho, U)
    );

    UEqn.relax();

    fvConstraints.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    mixture->surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvConstraints.constrain(U);
    }
