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

#include "phaseChangeThreePhaseMixture.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(phaseChangeThreePhaseMixture, 0);
    defineRunTimeSelectionTable(phaseChangeThreePhaseMixture, components);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseChangeThreePhaseMixture::phaseChangeThreePhaseMixture
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    incompressibleThreePhaseMixture(U, phi),
    threePhaseInterfaceProperties
    (
        static_cast<incompressibleThreePhaseMixture&>(*this)
    ),
    phaseChangeThreePhaseMixtureCoeffs_(optionalSubDict(type + "Coeffs")),
    pSat_("pSat", dimPressure, lookup("pSat"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


//source terms to be used in alpha equations.
Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeThreePhaseMixture::vDotAlphal() const
{
    volScalarField alphalCoeff(-1.0/rho1() + alpha1()*(1.0/rho1() - 1.0/rho2()));
    Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal();

    return Pair<tmp<volScalarField>>
    (
        alphalCoeff*mDotAlphal[0],
        alphalCoeff*mDotAlphal[1]
    );
}

//not used in current implementation
Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeThreePhaseMixture::vDotAlphav() const
{
    dimensionedScalar alphavCoeff(1.0/rho2());
    Pair<tmp<volScalarField>> mDotAlphal = this->mDotAlphal();

    return Pair<tmp<volScalarField>>
    (
        alphavCoeff*mDotAlphal[0],
        alphavCoeff*mDotAlphal[1]
    );
}

//source terms to be used in p_rgh equation.
Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::phaseChangeThreePhaseMixture::vDotP() const
{
    dimensionedScalar pCoeff(1.0/rho2() - 1.0/rho1());
    Pair<tmp<volScalarField>> mDotP = this->mDotP();

    return Pair<tmp<volScalarField>>(pCoeff*mDotP[0], pCoeff*mDotP[1]);
}


void Foam::phaseChangeThreePhaseMixture::correct()
{
    incompressibleThreePhaseMixture::correct();
    threePhaseInterfaceProperties::correct(); 
}


bool Foam::phaseChangeThreePhaseMixture::read()
{
    if (incompressibleThreePhaseMixture::read())
    {
        phaseChangeThreePhaseMixtureCoeffs_ = optionalSubDict(type() + "Coeffs");
        lookup("pSat") >> pSat_;

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
