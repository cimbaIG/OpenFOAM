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

#include "richardsHoxeyABLInletEpsFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

richardsHoxeyABLInletEpsFvPatchScalarField::richardsHoxeyABLInletEpsFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    uTau_(0),
    kappa_(0),
    z0_(0),
    z_(0, 1, 0)
{}


richardsHoxeyABLInletEpsFvPatchScalarField::richardsHoxeyABLInletEpsFvPatchScalarField
(
    const richardsHoxeyABLInletEpsFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    uTau_(ptf.uTau_),
    kappa_(ptf.kappa_),
    z0_(ptf.z0_),
    z_(ptf.z_)
{}


richardsHoxeyABLInletEpsFvPatchScalarField::richardsHoxeyABLInletEpsFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    uTau_(readScalar(dict.lookup("uTau"))),
    kappa_(readScalar(dict.lookup("kappa"))),
    z0_(readScalar(dict.lookup("z0"))),
    z_(dict.lookup("z"))
{
    if (mag(z_) < SMALL)
    {
        FatalErrorIn("richardsHoxeyABLInletEpsFvPatchScalarField(dict)")
            << "z given with zero size not correct"
            << abort(FatalError);
    }

    evaluate();
}


richardsHoxeyABLInletEpsFvPatchScalarField::richardsHoxeyABLInletEpsFvPatchScalarField
(
    const richardsHoxeyABLInletEpsFvPatchScalarField& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    uTau_(fcvpvf.uTau_),
    kappa_(fcvpvf.kappa_),
    z0_(fcvpvf.z0_),
    z_(fcvpvf.z_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void richardsHoxeyABLInletEpsFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    z_ /= mag(z_);

    const vectorField& cf = patch().Cf();

    const scalarField coord(cf & z_);

    scalarField Eps(coord.size());

    forAll(coord, i)
    {
	Eps[i] = Foam::pow(uTau_,3) / (kappa_ * (coord[i] + z0_));
    }

    scalarField::operator=(Eps);
}


// Write
void richardsHoxeyABLInletEpsFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("uTau")
        << uTau_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa")
        << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("z0")
        << z0_ << token::END_STATEMENT << nl;
    os.writeKeyword("z")
        << z_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, richardsHoxeyABLInletEpsFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
