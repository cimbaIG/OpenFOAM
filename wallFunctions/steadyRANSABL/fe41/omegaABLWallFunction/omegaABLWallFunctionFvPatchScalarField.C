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

\*---------------------------------------------------------------------------*/

#include "omegaABLWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void omegaABLWallFunctionFvPatchScalarField::checkType()
{
    if (!patch().isWall())
    {
        FatalErrorIn("omegaABLWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

omegaABLWallFunctionFvPatchScalarField::omegaABLWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF),
    UName_("U"),
    kName_("k"),
    GName_("RASModel::G"),
    nuName_("nu"),
    nutName_("nut"),
    Cmu_(0.09),
    E_(9.8),
    kappa_(0.41),
    z0_(p.size(), SMALL)
{
    checkType();
}


omegaABLWallFunctionFvPatchScalarField::omegaABLWallFunctionFvPatchScalarField
(
    const omegaABLWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedInternalValueFvPatchField<scalar>(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    kName_(ptf.kName_),
    GName_(ptf.GName_),
    nuName_(ptf.nuName_),
    nutName_(ptf.nutName_),
    Cmu_(ptf.Cmu_),
    E_(ptf.E_),
    kappa_(ptf.kappa_),
    z0_(ptf.z0_)
{
    checkType();
}


omegaABLWallFunctionFvPatchScalarField::omegaABLWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    GName_(dict.lookupOrDefault<word>("G", "RASModel::G")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu")),
    nutName_(dict.lookupOrDefault<word>("nut", "nut")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    z0_("z0", dict, p.size())
{
    checkType();
}


omegaABLWallFunctionFvPatchScalarField::omegaABLWallFunctionFvPatchScalarField
(
    const omegaABLWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedInternalValueFvPatchField<scalar>(owfpsf),
    UName_(owfpsf.UName_),
    kName_(owfpsf.kName_),
    GName_(owfpsf.GName_),
    nuName_(owfpsf.nuName_),
    nutName_(owfpsf.nutName_),
    Cmu_(owfpsf.Cmu_),
    E_(owfpsf.E_),
    kappa_(owfpsf.kappa_),
    z0_(owfpsf.z0_)
{
    checkType();
}


omegaABLWallFunctionFvPatchScalarField::omegaABLWallFunctionFvPatchScalarField
(
    const omegaABLWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(owfpsf, iF),
    UName_(owfpsf.UName_),
    kName_(owfpsf.kName_),
    GName_(owfpsf.GName_),
    nuName_(owfpsf.nuName_),
    nutName_(owfpsf.nutName_),
    Cmu_(owfpsf.Cmu_),
    E_(owfpsf.E_),
    kappa_(owfpsf.kappa_),
    z0_(owfpsf.z0_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void omegaABLWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // If G field is not present, execute zero gradient evaluation
    // HJ, 20/Mar/2011
    if (!db().foundObject<volScalarField>(GName_))
    {
        InfoIn("void omegaABLWallFunctionFvPatchScalarField::updateCoeffs()")
            << "Cannot access " << GName_ << " field for patch "
            << patch().name() << ".  Evaluating as zeroGradient"
            << endl;

        fvPatchScalarField::updateCoeffs();
        zeroGradientFvPatchScalarField::evaluate();

        return;
    }

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalarField& y = rasModel.y()[patch().index()];

    const scalar Cmu25 = pow(Cmu_, 0.25);

    volScalarField& G = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>(GName_));

    // Note: omega is now a refValue and set in fixedInternalValueFvPatchField
    // HJ, 3/Aug/2011
    scalarField& omega = refValue();

    const scalarField& k = db().lookupObject<volScalarField>(kName_);

    const scalarField& nuw =
        lookupPatchField<volScalarField, scalar>(nuName_);

    const scalarField& nutw =
        lookupPatchField<volScalarField, scalar>(nutName_);

    const fvPatchVectorField& Uw =
        lookupPatchField<volVectorField, vector>(UName_);

    vectorField n = patch().nf();

    const scalarField magGradUw = mag(Uw.snGrad());

    const labelList& faceCells = patch().faceCells();

    // Set omega and G
    forAll(nutw, faceI)
    {
        label faceCellI = faceCells[faceI];

        scalar omega_ = sqrt(k[faceCellI])/(Cmu25*kappa_*(y[faceI]+z0_[faceI]));

        omega[faceI] = omega_;

        G[faceCellI] = ( nutw[faceI] + nuw[faceI] ) * magGradUw[faceI]
		* Cmu25 * sqrt(k[faceCellI])
                / ( kappa_ * ( y[faceI] + z0_[faceI] ) );
    }

    // TODO: perform averaging for cells sharing more than one boundary face

    fixedInternalValueFvPatchField<scalar>::updateCoeffs();
}


void omegaABLWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedInternalValueFvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    writeEntryIfDifferent<word>(os, "G", "RASModel::G", GName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("z0") << z0_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    omegaABLWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
