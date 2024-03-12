surfaceScalarField phic(mag(phi/mesh.magSf()));
phic = min(mixture->cAlpha()*phic, max(phic));

if (nAlphaSubCycles > 1)
{    
    dimensionedScalar totalDeltaT = runTime.deltaT();
    
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), 0)
    );

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

//Update density after solution of alphas
volScalarField rhoNew(alpha1*rho1 + alpha2*rho2 + alpha3*rho3);
rho == rhoNew;
