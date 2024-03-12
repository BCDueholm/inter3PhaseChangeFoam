/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "SchnerrSauer3.H"
#include "mathematicalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeThreePhaseMixtures
{
    defineTypeNameAndDebug(SchnerrSauer3, 0);
    addToRunTimeSelectionTable
    (
        phaseChangeThreePhaseMixture,
        SchnerrSauer3,
        components
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeThreePhaseMixtures::SchnerrSauer3::SchnerrSauer3
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    phaseChangeThreePhaseMixture(typeName, U, phi),

    n_("n", dimless/dimVolume, phaseChangeThreePhaseMixtureCoeffs_),
    dNuc_("dNuc", dimLength, phaseChangeThreePhaseMixtureCoeffs_),
    Cc_("Cc", dimless, phaseChangeThreePhaseMixtureCoeffs_),
    Cv_("Cv", dimless, phaseChangeThreePhaseMixtureCoeffs_),

    p0_("0", pSat().dimensions(), 0.0)
{
    correct();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


//Bubble radius of nuclei in the liquid. 
Foam::tmp<Foam::volScalarField>
Foam::phaseChangeThreePhaseMixtures::SchnerrSauer3::Rb
(
    const volScalarField& limitedAlpha1
) const
{
    return pow
    (
        (3.0/(4.0*constant::mathematical::pi*n_))
        *max((1.0 + alphaNuc() - limitedAlpha1)
        /(max(limitedAlpha1,ROOTVSMALL))
        ,ROOTVSMALL),
        1.0/3.0
    );
}

//alpha of initial nuclei.
Foam::dimensionedScalar
Foam::phaseChangeThreePhaseMixtures::SchnerrSauer3::alphaNuc() const
{
    dimensionedScalar Vnuc = n_*constant::mathematical::pi*pow3(dNuc_)/6;
    return Vnuc/(1 + Vnuc);
}

// common part of numerator of both source terms.
Foam::tmp<Foam::volScalarField>
Foam::phaseChangeThreePhaseMixtures::SchnerrSauer3::pCoeff
(
    const volScalarField& p
) const
{
    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    volScalarField limitedAlpha2(min(max(alpha2_, scalar(0)), scalar(1)));
    volScalarField limitedAlpha3(min(max(alpha3_, scalar(0)), scalar(1)));
    
    volScalarField rho
    (
        limitedAlpha1*rho1() + limitedAlpha2*rho2() + limitedAlpha3*rho3()
    );

    return
        ((rho1()*rho2())/(rho + limitedAlpha3*(rho1()-rho3())))*
        sqrt(2.0/(3.0*rho1()))*sqrt(1.0/((mag(p - pSat()) + 0.001*pSat())));
}

//source terms used for alpha quations.
Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeThreePhaseMixtures::SchnerrSauer3::mDotAlphal() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeff(this->pCoeff(p));

    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    volScalarField limitedAlpha2(min(max(alpha2_, scalar(0)), scalar(1)));
    volScalarField limitedAlpha3(min(max(alpha3_, scalar(0)), scalar(1)));
    
    volScalarField rho
    (
        limitedAlpha1*rho1() + limitedAlpha2*rho2() + limitedAlpha3*rho3()
    );
    volScalarField densityfrac
    (
        (rho+limitedAlpha3*(rho2()-rho3()))/(rho+limitedAlpha3*(rho1()-rho3()))
    );

    return Pair<tmp<volScalarField>>
    (
        (-Cc_)
        *((3.0*pCoeff*max(p - pSat(), p0_))
        /(Rb(limitedAlpha1) + pow(Rb(limitedAlpha1),4.0)
        *(4.0/3.0)*n_*constant::mathematical::pi*densityfrac)),

        (-Cv_)
        *((4.0*n_*constant::mathematical::pi*pow(Rb(limitedAlpha1),2.0)
        *pCoeff*min(p - pSat(), p0_))
        /(1.0 + pow(Rb(limitedAlpha1),3.0)
        *(4.0/3.0)*n_*constant::mathematical::pi*densityfrac))
    );
}

//source terms used for p_rgh equation
Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeThreePhaseMixtures::SchnerrSauer3::mDotP() const
{
    const volScalarField& p = alpha1_.db().lookupObject<volScalarField>("p");
    volScalarField pCoeff(this->pCoeff(p));

    volScalarField limitedAlpha1(min(max(alpha1_, scalar(0)), scalar(1)));
    volScalarField limitedAlpha2(min(max(alpha2_, scalar(0)), scalar(1)));
    volScalarField limitedAlpha3(min(max(alpha3_, scalar(0)), scalar(1)));
    
    volScalarField rho
    (
        limitedAlpha1*rho1() + limitedAlpha2*rho2() + limitedAlpha3*rho3()
    );
    volScalarField densityfrac
    (
        (rho+limitedAlpha3*(rho2()-rho3()))/(rho+limitedAlpha3*(rho1()-rho3()))
    );
    
    volScalarField apCoeff(limitedAlpha1*pCoeff);

    return Pair<tmp<volScalarField>>
    (
        (-Cc_)*(1.0 - limitedAlpha1 - limitedAlpha3)
        *((3.0)
        /(Rb(limitedAlpha1) + pow(Rb(limitedAlpha1),4.0)
        *(4.0/3.0)*n_*constant::mathematical::pi*densityfrac))
        *pos0(p - pSat())*pCoeff,

        (Cv_)*(limitedAlpha1)
        *((4.0*n_*constant::mathematical::pi*pow(Rb(limitedAlpha1),2.0))
        /(1.0 + pow(Rb(limitedAlpha1),3.0)
        *(4.0/3.0)*n_*constant::mathematical::pi*densityfrac))
        *neg(p - pSat())*pCoeff
    );
}


void Foam::phaseChangeThreePhaseMixtures::SchnerrSauer3::correct()
{
    phaseChangeThreePhaseMixture::correct();
}


bool Foam::phaseChangeThreePhaseMixtures::SchnerrSauer3::read()
{
    if (phaseChangeThreePhaseMixture::read())
    {
        phaseChangeThreePhaseMixtureCoeffs_ = optionalSubDict(type() + "Coeffs");

        phaseChangeThreePhaseMixtureCoeffs_.lookup("n") >> n_;
        phaseChangeThreePhaseMixtureCoeffs_.lookup("dNuc") >> dNuc_;
        phaseChangeThreePhaseMixtureCoeffs_.lookup("Cc") >> Cc_;
        phaseChangeThreePhaseMixtureCoeffs_.lookup("Cv") >> Cv_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
