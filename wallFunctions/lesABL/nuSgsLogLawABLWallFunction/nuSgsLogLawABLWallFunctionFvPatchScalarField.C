/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "nuSgsLogLawABLWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace LESModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nuSgsLogLawABLWallFunctionFvPatchScalarField::
nuSgsLogLawABLWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    nuName_("nu"),
    kappa_(0.41),
    z0_(0.001)
{}


nuSgsLogLawABLWallFunctionFvPatchScalarField::
nuSgsLogLawABLWallFunctionFvPatchScalarField
(
    const nuSgsLogLawABLWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    nuName_(ptf.nuName_),
    kappa_(ptf.kappa_),
    z0_(ptf.z0_)
{}


nuSgsLogLawABLWallFunctionFvPatchScalarField::
nuSgsLogLawABLWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu")),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    z0_(dict.lookupOrDefault<scalar>("z0", 0.001))
{}


nuSgsLogLawABLWallFunctionFvPatchScalarField::
nuSgsLogLawABLWallFunctionFvPatchScalarField
(
    const nuSgsLogLawABLWallFunctionFvPatchScalarField& nwfpsf
)
:
    fixedValueFvPatchScalarField(nwfpsf),
    UName_(nwfpsf.UName_),
    nuName_(nwfpsf.nuName_),
    kappa_(nwfpsf.kappa_),
    z0_(nwfpsf.z0_)
{}


nuSgsLogLawABLWallFunctionFvPatchScalarField::
nuSgsLogLawABLWallFunctionFvPatchScalarField
(
    const nuSgsLogLawABLWallFunctionFvPatchScalarField& nwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(nwfpsf, iF),
    UName_(nwfpsf.UName_),
    nuName_(nwfpsf.nuName_),
    kappa_(nwfpsf.kappa_),
    z0_(nwfpsf.z0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void nuSgsLogLawABLWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    const scalarField& ry = patch().deltaCoeffs();

    const fvPatchVectorField& U =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    scalarField magUp = mag(U.patchInternalField() - U);

    const scalarField& nuw =
        patch().lookupPatchField<volScalarField, scalar>(nuName_);

    scalarField& nuSgsw = *this;

    scalarField magFaceGradU = mag(U.snGrad());

    forAll(nuSgsw, facei)
    {
        scalar magUpara = magUp[facei];

        scalar utau = ( kappa_ * magUpara ) / ( log( (1/ry[facei]) * z0_ ) );
        
        nuSgsw[facei] =
            max(sqr(max(utau, 0))/magFaceGradU[facei] - nuw[facei], 0.0);
    }

    fixedValueFvPatchScalarField::evaluate();
}


void nuSgsLogLawABLWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("z0") << z0_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nuSgsLogLawABLWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LESModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
