/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::incompressible::RASModels::windVeeringABLkEpsilon

Description
    Modified standard k-epsilon turbulence model for incompressible flows
    that allows for wind veering effects in ABL, according to Apsley and
    Castro (1997), Feng and Gu (2019) and Feng and Gu (2020).

    The default model coefficients correspond to the following:
    @verbatim
        windVeeringABLkEpsilonCoeffs
        {
            Cmu         0.09;
            C1          1.44;
            C2          1.92;
            sigmaEps    1.3;
	    sigmaK	1.;
	    lmax	28;
        }
    @endverbatim

SourceFiles
    windVeeringABLkEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef windVeeringABLkEpsilon_H
#define windVeeringABLkEpsilon_H

#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                     Class windVeeringABLkEpsilon Declaration
\*---------------------------------------------------------------------------*/

class windVeeringABLkEpsilon
:
    public RASModel
{
    // Private data

        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C2_;
            dimensionedScalar sigmaEps_;
	    dimensionedScalar sigmaK_;
	    dimensionedScalar lmax_;


        // Fields

            volScalarField k_;
            volScalarField epsilon_;
            volScalarField nut_;

public:

    //- Runtime type information
    TypeName("windVeeringABLkEpsilon");

    // Constructors

        //- Construct from components
        windVeeringABLkEpsilon
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport
        );


    //- Destructor
    virtual ~windVeeringABLkEpsilon()
    {}


    // Member Functions

        //- Return the turbulence viscosity
        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for k
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DkEff", nut_/sigmaK_ + nu())
            );
        }

        //- Return the effective diffusivity for epsilon
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField("DepsilonEff", nut_/sigmaEps_ + nu())
            );
        }

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the Reynolds stress tensor
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        //- Read RASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
