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
    Boundary condition specifies a turb. kin. energy profile at inlet
    (fixed value) for simulation of homogeneous pressure driven ABL. 
    BC is presented in: 
    Richards, P.J., Norris S.E. (2015) 
    Appropriate boundary conditions for a pressure driven boundary layer. 
    Journal of Wind Engineering and Industrial Aerodynamics (JWEIA) 
    142, 43-52.
	
    inlet 
    {
        type	pressureDrivenABLInletTke;
        uTau	0.3186;
        hd	400;
        Ck1	0.921;
        Ck2	3.533;
        Ck3	-1.926;
        Ck4	0.805;
        z	(0 1 0);
    }

SourceFiles
    pressureDrivenABLInletVelocityFvPatchVectorField.C

Author
	October, 2019
    Mihael Cindori, Faculty of Mechanical Engineering and Naval
	Architecture, University of Zagreb, Croatia.  
	All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef pressureDrivenABLInletTkeFvPatchScalarField_H
#define pressureDrivenABLInletTkeFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class pressureDrivenABLInletTkeFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class pressureDrivenABLInletTkeFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Friction velocity
        const scalar uTau_;

        //- Domain height
        const scalar hd_;

        //- Tke profile constant
        const scalar Ck1_;

		//- Tke profile constant
        const scalar Ck2_;

		//- Tke profile constant
        const scalar Ck3_;

		//- Tke profile constant
        const scalar Ck4_;

        //- Direction of the z-coordinate
        vector z_;

public:

    //- Runtime type information
    TypeName("pressureDrivenABLInletTke");


    // Constructors

        //- Construct from patch and internal field
        pressureDrivenABLInletTkeFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureDrivenABLInletTkeFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given pressureDrivenABLInletTkeFvPatchScalarField
        //  onto a new patch
        pressureDrivenABLInletTkeFvPatchScalarField
        (
            const pressureDrivenABLInletTkeFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new pressureDrivenABLInletTkeFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pressureDrivenABLInletTkeFvPatchScalarField
        (
            const pressureDrivenABLInletTkeFvPatchScalarField&,
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
                new pressureDrivenABLInletTkeFvPatchScalarField(*this, iF)
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
