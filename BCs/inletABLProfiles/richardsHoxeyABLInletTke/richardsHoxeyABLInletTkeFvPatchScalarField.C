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

#include "richardsHoxeyABLInletTkeFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

richardsHoxeyABLInletTkeFvPatchScalarField::richardsHoxeyABLInletTkeFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    uTau_(0),
    Cmu_(0), 
    z_(0, 1, 0)
{}


richardsHoxeyABLInletTkeFvPatchScalarField::richardsHoxeyABLInletTkeFvPatchScalarField
(
    const richardsHoxeyABLInletTkeFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    uTau_(ptf.uTau_),
    Cmu_(ptf.Cmu_),
    z_(ptf.z_)
{}


richardsHoxeyABLInletTkeFvPatchScalarField::richardsHoxeyABLInletTkeFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    uTau_(readScalar(dict.lookup("uTau"))),
    Cmu_(readScalar(dict.lookup("Cmu"))),
    z_(dict.lookup("z"))
{
    if (mag(z_) < SMALL)
    {
        FatalErrorIn("richardsHoxeyABLInletTkeFvPatchScalarField(dict)")
            << "z given with zero size not correct"
            << abort(FatalError);
    }

    evaluate();
}


richardsHoxeyABLInletTkeFvPatchScalarField::richardsHoxeyABLInletTkeFvPatchScalarField
(
    const richardsHoxeyABLInletTkeFvPatchScalarField& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    uTau_(fcvpvf.uTau_),
    Cmu_(fcvpvf.Cmu_),
    z_(fcvpvf.z_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void richardsHoxeyABLInletTkeFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    z_ /= mag(z_);

    const vectorField& cf = patch().Cf();

    const scalarField coord(cf & z_);

    scalarField tke(coord.size());

    forAll(coord, i)
    {
        tke[i] = sqr(uTau_) / sqrt(Cmu_);
    }

    scalarField::operator=(tke);
}


// Write
void richardsHoxeyABLInletTkeFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("uTau")
        << uTau_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cmu")
        << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("z")
        << z_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, richardsHoxeyABLInletTkeFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
