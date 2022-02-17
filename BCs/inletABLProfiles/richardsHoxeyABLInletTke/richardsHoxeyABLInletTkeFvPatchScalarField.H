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
    Inlet turb. kin. energy profile for homogeneous shear stress driven ABL
    simulation. BC is presented in:
    Richards, P.J., Hoxey, R.P. (1993)
    Appropriate boundary conditions for computational wind engineering 
    models using the k-epsilon turbulence model.
    Journal of Wind Engineering and Industrial Aerodynamics (JWEIA) 
    46 & 47, 145-153.
	
    inlet
    {
        type	richardsHoxeyABLInletTke;
        uTau	0.3186;
	Cmu	0.09;
	z	(0 1 0);
    }
		
SourceFiles
    richardsHoxeyABLInletVelocityFvPatchVectorField.C

Author
    October, 2019
    Mihael Cindori, Faculty of Mechanical Engineering and Naval
    Architecture, University of Zagreb, Croatia.  
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef richardsHoxeyABLInletTkeFvPatchScalarField_H
#define richardsHoxeyABLInletTkeFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class richardsHoxeyABLInletTkeFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class richardsHoxeyABLInletTkeFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Friction velocity
        const scalar uTau_;

        //- Cmu constant
        const scalar Cmu_;

        //- Direction of the z-coordinate
        vector z_;

public:

    //- Runtime type information
    TypeName("richardsHoxeyABLInletTke");


    // Constructors

        //- Construct from patch and internal field
        richardsHoxeyABLInletTkeFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        richardsHoxeyABLInletTkeFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given richardsHoxeyABLInletTkeFvPatchScalarField
        //  onto a new patch
        richardsHoxeyABLInletTkeFvPatchScalarField
        (
            const richardsHoxeyABLInletTkeFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new richardsHoxeyABLInletTkeFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        richardsHoxeyABLInletTkeFvPatchScalarField
        (
            const richardsHoxeyABLInletTkeFvPatchScalarField&,
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
                new richardsHoxeyABLInletTkeFvPatchScalarField(*this, iF)
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
