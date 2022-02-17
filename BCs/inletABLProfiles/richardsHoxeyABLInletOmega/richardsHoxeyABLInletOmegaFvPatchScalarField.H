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

Description
    Inlet dissipation rate profile for homogeneous shear stress driven ABL
    simulation. BC is presented in:
    Richards, P.J., Hoxey, R.P. (1993)
    Appropriate boundary conditions for computational wind engineering 
    models using the k-epsilon turbulence model.
    Journal of Wind Engineering and Industrial Aerodynamics (JWEIA) 
    46 & 47, 145-153.
	
    inlet
    {
        type	richardsHoxeyABLInletOmega;
        uTau	0.3186;
        kappa	0.4327;
	Cmu	0.09;
        z0	0.004;
        z	(0 1 0);
    }

SourceFiles
    richardsHoxeyABLInletOmegaFvPatchVectorField.C

Author
    October, 2020
    Mihael Cindori, Faculty of Mechanical Engineering and Naval
    Architecture, University of Zagreb, Croatia.  
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef richardsHoxeyABLInletOmegaFvPatchScalarField_H
#define richardsHoxeyABLInletOmegaFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class richardsHoxeyABLInletOmegaFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class richardsHoxeyABLInletOmegaFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Friction velocity
        const scalar uTau_;

        //- Von Karman constant
        const scalar kappa_;

	//- Cmu turbulence constant
	const scalar Cmu_;

        //- Roughness length
        const scalar z0_;

        //- Direction of the z-coordinate
        vector z_;

public:

    //- Runtime type information
    TypeName("richardsHoxeyABLInletOmega");


    // Constructors

        //- Construct from patch and internal field
        richardsHoxeyABLInletOmegaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        richardsHoxeyABLInletOmegaFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given richardsHoxeyABLInletOmegaFvPatchScalarField
        //  onto a new patch
        richardsHoxeyABLInletOmegaFvPatchScalarField
        (
            const richardsHoxeyABLInletOmegaFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new richardsHoxeyABLInletOmegaFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        richardsHoxeyABLInletOmegaFvPatchScalarField
        (
            const richardsHoxeyABLInletOmegaFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new richardsHoxeyABLInletOmegaFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
