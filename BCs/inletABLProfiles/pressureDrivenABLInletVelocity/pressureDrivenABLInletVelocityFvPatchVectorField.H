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
    pressureDrivenABLInletVelocityFvPatchVectorField

Description
    Boundary condition specifies an extended log-law profile at inlet
    (fixed value). BC is presented in: 
    Richards, P.J., Norris S.E. (2015) 
    Appropriate boundary conditions for a pressure driven boundary layer. 
    Journal of Wind Engineering and Industrial Aerodynamics (JWEIA) 
    142, 43-52.
	
    inlet
    {
        type	pressureDrivenABLInletVelocity;
        kappa	0.4327; //default value
        Cu1	0.528;
        Cu2	0.385;
        Cu3	-1.09;  // default value is 0
        Cu4	0.243;  // default value is 0
        hd	400;    // domain height
	uTau	0.3186; 
        z0      0.004;
        n	(1 0 0); // flow direction
        z	(0 1 0); // height direction
    }
		

SourceFiles
    pressureDrivenABLInletVelocityFvPatchVectorField.C

Author
    October, 2019
    Mihael Cindori, Faculty of Mechanical Engineering and Naval
    Architecture, University of Zagreb, Croatia.  
    All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef pressureDrivenABLInletVelocityFvPatchVectorField_H
#define pressureDrivenABLInletVelocityFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class pressureDrivenABLInletVelocityFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class pressureDrivenABLInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Kappa constant (von Karman constant)
        const scalar kappa_;

        //- Velocity profile constant
        const scalar Cu1_;

		//- Velocity profile constant
        const scalar Cu2_;

		//- Velocity profile constant
        const scalar Cu3_;

		//- Velocity profile constant
        const scalar Cu4_;

        //- Domain height
        const scalar hd_;

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
    TypeName("pressureDrivenABLInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        pressureDrivenABLInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureDrivenABLInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given pressureDrivenABLInletVelocityFvPatchVectorField
        //  onto a new patch
        pressureDrivenABLInletVelocityFvPatchVectorField
        (
            const pressureDrivenABLInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new pressureDrivenABLInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pressureDrivenABLInletVelocityFvPatchVectorField
        (
            const pressureDrivenABLInletVelocityFvPatchVectorField&,
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
                new pressureDrivenABLInletVelocityFvPatchVectorField(*this, iF)
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
