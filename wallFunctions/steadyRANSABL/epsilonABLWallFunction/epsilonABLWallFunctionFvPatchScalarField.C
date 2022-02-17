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

#include "epsilonABLWallFunctionFvPatchScalarField.H"
#include "RASModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void epsilonABLWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorIn("epsilonABLWallFunctionFvPatchScalarField::checkType()")
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

epsilonABLWallFunctionFvPatchScalarField::epsilonABLWallFunctionFvPatchScalarField
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
    kappa_(0.41),
    z0_(p.size(), SMALL)
{
    checkType();
}


epsilonABLWallFunctionFvPatchScalarField::epsilonABLWallFunctionFvPatchScalarField
(
    const epsilonABLWallFunctionFvPatchScalarField& ptf,
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
    kappa_(ptf.kappa_),
    z0_(ptf.z0_, mapper)
{
    checkType();
}


epsilonABLWallFunctionFvPatchScalarField::epsilonABLWallFunctionFvPatchScalarField
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
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    z0_("z0", dict, p.size())
{
    checkType();
}


epsilonABLWallFunctionFvPatchScalarField::epsilonABLWallFunctionFvPatchScalarField
(
    const epsilonABLWallFunctionFvPatchScalarField& ewfpsf
)
:
    fixedInternalValueFvPatchField<scalar>(ewfpsf),
    UName_(ewfpsf.UName_),
    kName_(ewfpsf.kName_),
    GName_(ewfpsf.GName_),
    nuName_(ewfpsf.nuName_),
    nutName_(ewfpsf.nutName_),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    z0_(ewfpsf.z0_)
{
    checkType();
}


epsilonABLWallFunctionFvPatchScalarField::epsilonABLWallFunctionFvPatchScalarField
(
    const epsilonABLWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(ewfpsf, iF),
    UName_(ewfpsf.UName_),
    kName_(ewfpsf.kName_),
    GName_(ewfpsf.GName_),
    nuName_(ewfpsf.nuName_),
    nutName_(ewfpsf.nutName_),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    z0_(ewfpsf.z0_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void epsilonABLWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");
    const scalarField& y = rasModel.y()[patch().index()];

    const scalar Cmu25 = pow(Cmu_, 0.25);

    volScalarField& G = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>(GName_));

    volScalarField& epsilon = const_cast<volScalarField&>
        (db().lookupObject<volScalarField>(dimensionedInternalField().name()));

    const volScalarField& k = db().lookupObject<volScalarField>(kName_);

    const scalarField& nuw =
        patch().lookupPatchField<volScalarField, scalar>(nuName_);

    const scalarField& nutw =
        patch().lookupPatchField<volScalarField, scalar>(nutName_);

    const fvPatchVectorField& Uw =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const scalarField magGradUw = mag(Uw.snGrad());
    const scalar srCmu = sqrt(Cmu_);

    // Set epsilon and G
    forAll(nutw, faceI)
    {
        const label faceCellI = patch().faceCells()[faceI];

	const scalar U_tau = Cmu25 * sqrt(k[faceCellI]);
	const scalar z0 = z0_[faceI];
	const scalar dUdy = U_tau / (kappa_ * (y[faceI] + z0));

	//- set epsilon in the near-wall cell
        epsilon[faceCellI] = srCmu*k[faceCellI]*dUdy-nuw[faceI]*sqr(dUdy);

	//- set the production term in the near-wall cell
	const scalar nuEff = nuw[faceI] + nutw[faceI];
	G[faceCellI] = nutw[faceI]*magGradUw[faceI]*dUdy-nuw[faceI]*sqr(dUdy);
    }

    // TODO: perform averaging for cells sharing more than one boundary face

    fixedInternalValueFvPatchField<scalar>::updateCoeffs();
}


void epsilonABLWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedInternalValueFvPatchField<scalar>::evaluate(commsType);
}


void epsilonABLWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedInternalValueFvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    writeEntryIfDifferent<word>(os, "G", "RASModel::G", GName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    z0_.writeEntry("z0", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    epsilonABLWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
