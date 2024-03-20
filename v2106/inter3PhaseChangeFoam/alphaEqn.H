{
    word alphaScheme("div(phi,alpha)");

    surfaceScalarField phir
    (
        IOobject
        (
            "phir",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mixture->cAlpha()*mag(phi/mesh.magSf())*mixture->nHatf()
    );
    
    //construct the source terms from phaseChangeThreePhaseMixture
    Pair<tmp<volScalarField>> vDotAlphal = mixture->vDotAlphal();
    
    //Assign source terms to variables.
    const volScalarField& vDotcAlphal = vDotAlphal[0]();
    const volScalarField& vDotvAlphal = vDotAlphal[1]();
    const volScalarField vDotvmcAlphal(vDotvAlphal - vDotcAlphal);
      
    //create temporary variables
    tmp<volScalarField> rhoNew;
    
    tmp<surfaceScalarField> talphaPhi1;
    tmp<surfaceScalarField> talphaPhi3;

    //If MULEScorr is set to yes in fvOptions start by solving alpha1 and alpha3
    if (MULESCorr)
    {
        //create fvMatrix of alpha1
        fvScalarMatrix alpha1Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha1)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha1)
          - fvm::Sp(divU, alpha1)
          ==
            fvm::Sp(vDotvmcAlphal, alpha1)
          + vDotcAlphal*(1.0-alpha3)
        );
        
        //solve alpha1
        alpha1Eqn.solve();

        //output info about solution
        Info<< "Phase-1 volume fraction = "
            << alpha1.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
            << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
            << endl;

        //store flux of alpha1
        talphaPhi1 = alpha1Eqn.flux();
        
        //create fvMatrix of alpha3
        fvScalarMatrix alpha3Eqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha3)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha3)
        );
        
        //solve alpha3
        alpha3Eqn.solve();
        
        //output info about solution
        Info<< "Phase-3 volume fraction = "
        << alpha3.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha3.name() << ") = " << min(alpha3).value()
        << "  Max(" << alpha3.name() << ") = " << max(alpha3).value()
        << endl;

        //store flux of alpha3
        talphaPhi3 = alpha3Eqn.flux();
        
        //calculate alpha2 from alpha1 and alpha3
        alpha2 = 1.0 - alpha1 - alpha3;
        
         //correct from phaseChangeThreePhaseMixture.C
        mixture->correct(); 
    }

    volScalarField alpha10("alpha10", alpha1);
    volScalarField alpha30("alpha30", alpha3);
    
    surfaceScalarField alphaPhiCorr1("alphaPhiCorr1", phi);
    surfaceScalarField alphaPhiCorr3("alphaPhiCorr3", phi);
    
    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {     

        tmp<surfaceScalarField> talphaPhiCorr1
        (
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
        );

        tmp<surfaceScalarField> talphaPhiCorr3
        (
            fvc::flux
            (
                phi,
                alpha3,
                alphaScheme
            )
        );
        //Do flux corrections of alpha1 and alpha3 with MULES::correct.
        if (MULESCorr)
        {
            talphaPhiCorr1.ref() -= talphaPhi1();

            volScalarField alpha100("alpha100", alpha10);
            alpha10 = alpha1;

            MULES::correct
            (
                geometricOneField(),
                alpha1,
                talphaPhi1(),
                talphaPhiCorr1.ref(),
                vDotvmcAlphal,
                (
                  - vDotvmcAlphal*alpha10
                )(),
                oneField(),
                zeroField()
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                talphaPhi1.ref() += talphaPhiCorr1();
            }
            else
            {
                alpha1 = 0.5*alpha1 + 0.5*alpha10;
                talphaPhi1.ref() += 0.5*talphaPhiCorr1();
            }
            
            //######## phase 3 ########
            talphaPhiCorr3.ref() -= talphaPhi3();

            volScalarField alpha300("alpha300", alpha30);
            alpha30 = alpha3;

            MULES::correct
            (
                geometricOneField(),
                alpha3,
                talphaPhi3(),
                talphaPhiCorr3.ref(),
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                talphaPhi3.ref() += talphaPhiCorr3();
            }
            else
            {
                alpha3 = 0.5*alpha3 + 0.5*alpha30;
                talphaPhi3.ref() += 0.5*talphaPhiCorr3();
            }
            
        }
        //If mulesCorr is set to "no" explicitSolve alpha1 and 3 using MULES. The is done in the same manner as above, first solving alpha1 then alpha3 then alpha2 from the 2 others.
        else
        {
        

            surfaceScalarField alphaPhi1
            (
                fvc::flux
                (
                    phi,
                    alpha1,
                    alphaScheme
                )

            );

            surfaceScalarField alphaPhi3
            (
                fvc::flux
                (
                    phi,
                    alpha3,
                    alphaScheme
                )
            );
            
            
            volScalarField Sp1 = vDotvmcAlphal;
            volScalarField Su1(vDotcAlphal*(scalar(1)-alpha3) + divU*min(alpha1, scalar(1)));
            
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha1,
                phi,
                alphaPhi1,
                Sp1,
                Su1,
                oneField(),
                zeroField()
            );
            
            alphaPhiCorr1 = alphaPhi1;
            
           
            
            MULES::explicitSolve
            (
                geometricOneField(),
                alpha3,
                phi,
                alphaPhi3,
                zeroField(),
                zeroField(),
                oneField(),
                zeroField()
            );
            
            alphaPhiCorr3 = alphaPhi3;
        }
        
        alpha2 = 1.0 - alpha1 - alpha3;
        
        mixture->correct();
    }
    
    //Calculate rhoPhi for the momentum equation.
    if (MULESCorr)
    {
        rhoPhi = phi*rho2 + talphaPhi1*(rho1-rho2) + talphaPhi3*(rho3-rho2);
    }
    else
    {
        rhoPhi = phi*rho2 + alphaPhiCorr1*(rho1-rho2) + alphaPhiCorr3*(rho3-rho2);
    }
    
    //output final solution of alphas
    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;

    Info<< "Vapour phase volume fraction = "
        << alpha2.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha2.name() << ") = " << min(alpha2).value()
        << "  Max(" << alpha2.name() << ") = " << max(alpha2).value()
        << endl;
        
    Info<< "Air phase volume fraction = "
        << alpha3.weightedAverage(mesh.V()).value()
        << "  Min(" << alpha3.name() << ") = " << min(alpha3).value()
        << "  Max(" << alpha3.name() << ") = " << max(alpha3).value()
        << endl;
}
