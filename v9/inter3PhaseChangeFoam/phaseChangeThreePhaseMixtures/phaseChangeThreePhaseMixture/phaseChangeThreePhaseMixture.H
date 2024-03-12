/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

#ifndef phaseChangeThreePhaseMixture_H
#define phaseChangeThreePhaseMixture_H

#include "incompressibleThreePhaseMixture.H"
#include "threePhaseInterfaceProperties.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "autoPtr.H"
#include "Pair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phaseChangeThreePhaseMixture
:
    public incompressibleThreePhaseMixture,
    public threePhaseInterfaceProperties
{
protected:

    // Protected data

        dictionary phaseChangeThreePhaseMixtureCoeffs_;

        //- Saturation vapour pressure
        dimensionedScalar pSat_;


public:

    //- Runtime type information
    TypeName("phaseChangeThreePhaseMixture");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            phaseChangeThreePhaseMixture,
            components,
            (
                const volVectorField& U,
                const surfaceScalarField& phi
            ),
            (U, phi)
        );


    // Constructors

        //- Construct from components
        phaseChangeThreePhaseMixture
        (
            const word& type,
            const volVectorField& U,
            const surfaceScalarField& phi
        );

        //- Disallow default bitwise copy construction
        phaseChangeThreePhaseMixture(const phaseChangeThreePhaseMixture&);


    // Selectors

        //- Return a reference to the selected phaseChange model
        static autoPtr<phaseChangeThreePhaseMixture> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~phaseChangeThreePhaseMixture()
    {}


    // Member Functions

        //- Return const-access to the saturation vapour pressure
        const dimensionedScalar& pSat() const
        {
            return pSat_;
        }

        virtual Pair<tmp<volScalarField>> mDotAlphal() const = 0;

        virtual Pair<tmp<volScalarField>> mDotP() const = 0;

        Pair<tmp<volScalarField>> vDotAlphal() const;
        
        //not used in current implementation
        Pair<tmp<volScalarField>> vDotAlphav() const;

        Pair<tmp<volScalarField>> vDotP() const;

        //- Correct the phaseChange model
        virtual void correct() = 0;

        //- Read the transportProperties dictionary and update
        virtual bool read() = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const phaseChangeThreePhaseMixture&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
