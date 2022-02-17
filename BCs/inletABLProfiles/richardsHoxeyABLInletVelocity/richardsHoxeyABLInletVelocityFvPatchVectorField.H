/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    richardsHoxeyABLInletVelocityFvPatchVectorField

Description
    Inlet velocity profile for homogeneous shear stress driven ABL
    simulation. BC is presented in:
    Richards, P.J., Hoxey, R.P. (1993)
    Appropriate boundary conditions for computational wind engineering 
    models using the k-epsilon turbulence model.
    Journal of Wind Engineering and Industrial Aerodynamics (JWEIA) 
    46 & 47, 145-153.
	
    inlet
    {
        type     richardsHoxeyABLInletVelocity;
        kappa    0.4327; //default value
        uTau	 0.3186; 
        z0       0.004;
        n        (1 0 0); // flow direction
        z        (0 1 0); // height direction
    }
		

SourceFiles
    richardsHoxeyABLInletVelocityFvPatchVectorField.C

Author
    October, 2019
    Mihael Cindori, Faculty of Mechanical Engineering and Naval
    Architecture, University of Zagreb, Croatia.  
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef richardsHoxeyABLInletVelocityFvPatchVectorField_H
#define richardsHoxeyABLInletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class richardsHoxeyABLInletVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class richardsHoxeyABLInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Kappa constant (von Karman constant)
        const scalar kappa_;

        //- Friction velocity
        const scalar uTau_;

        //- Roughness length
        const scalar z0_;

        //- Flow direction
        vector n_;

        //- Direction of the z-coordinate
        vector z_;

public:

    //- Runtime type information
    TypeName("richardsHoxeyABLInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        richardsHoxeyABLInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        richardsHoxeyABLInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given richardsHoxeyABLInletVelocityFvPatchVectorField
        //  onto a new patch
        richardsHoxeyABLInletVelocityFvPatchVectorField
        (
            const richardsHoxeyABLInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new richardsHoxeyABLInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        richardsHoxeyABLInletVelocityFvPatchVectorField
        (
            const richardsHoxeyABLInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new richardsHoxeyABLInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

	     //- Return flow direction
         const vector& n() const
         {
             return n_;
         }

         //- Return z direction
         const vector& z() const
         {
             return z_;
         }

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
