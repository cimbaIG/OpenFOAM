/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::incompressible::RASModels::nutABLWallFunctionFvPatchScalarField

Description
    Boundary condition for turbulent (kinematic) viscosity when using wall
    functions to simulate neutrally-stratified unidirectional atmospheric 
    boundary layer.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nutABLWallFunction;
	z0		uniform 0.001;
	value		0.001;
    }
    \endverbatim

See also
    Foam::nutWallFunctionFvPatchScalarField

SourceFiles
    nutABLWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutABLWallFunctionFvPatchScalarField_H
#define nutABLWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
             Class nutABLWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutABLWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    // Protected data
    
	//- Roughness height
	scalarField z0_;

    // Protected Member Functions

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> calcNut() const;


public:

    //- Runtime type information
    TypeName("nutABLWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutABLWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutABLWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutABLWallFunctionFvPatchScalarField
        //  onto a new patch
        nutABLWallFunctionFvPatchScalarField
        (
            const nutABLWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutABLWallFunctionFvPatchScalarField
        (
            const nutABLWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutABLWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutABLWallFunctionFvPatchScalarField
        (
            const nutABLWallFunctionFvPatchScalarField&,
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
                new nutABLWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Calculate and return the yPlus at the boundary
        virtual tmp<scalarField> yPlus() const;

	// Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );

        // I-O

            //- Write            
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
