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

#ifndef SchnerrSauer3_H
#define SchnerrSauer3_H

#include "phaseChangeThreePhaseMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseChangeThreePhaseMixtures
{

/*---------------------------------------------------------------------------*\
                              Class SchnerrSauer3
\*---------------------------------------------------------------------------*/

class SchnerrSauer3
:
    public phaseChangeThreePhaseMixture
{
    // Private Data

        //- Bubble number density
        dimensionedScalar n_;

        //- Nucleation site diameter
        dimensionedScalar dNuc_;

        //- Condensation rate coefficient
        dimensionedScalar Cc_;

        //- Vapourisation rate coefficient
        dimensionedScalar Cv_;

        dimensionedScalar p0_;

        //- Nucleation site volume-fraction
        dimensionedScalar alphaNuc() const;

        //bubble radius
        tmp<volScalarField>Rb(const volScalarField& limitedAlpha1) const;

        //- Part of the condensation and vapourisation rates
        tmp<volScalarField> pCoeff(const volScalarField& p) const;


public:

    //- Runtime type information
    TypeName("SchnerrSauer3");


    // Constructors

        //- Construct from components
        SchnerrSauer3
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~SchnerrSauer3()
    {}


    // Member Functions

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient use in alpha equations
        virtual Pair<tmp<volScalarField>> mDotAlphal() const;

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient use in p_rhg equations
        virtual Pair<tmp<volScalarField>> mDotP() const;

        //- Correct the SchnerrSauer3 phaseChange model
        virtual void correct();

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace phaseChangeThreePhaseMixtures
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
