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

\*---------------------------------------------------------------------------*/

#include "richardsHoxeyABLInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

richardsHoxeyABLInletVelocityFvPatchVectorField::richardsHoxeyABLInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    kappa_(0.4327),
    uTau_(0),
    z0_(0),
    n_(pTraits<vector>::zero),
    z_(pTraits<vector>::zero)
{}


richardsHoxeyABLInletVelocityFvPatchVectorField::richardsHoxeyABLInletVelocityFvPatchVectorField
(
    const richardsHoxeyABLInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    kappa_(ptf.kappa_),
    uTau_(ptf.uTau_),
    z0_(ptf.z0_),
    n_(ptf.n_),
    z_(ptf.z_)
{}


richardsHoxeyABLInletVelocityFvPatchVectorField::richardsHoxeyABLInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.4327)),
    uTau_(readScalar(dict.lookup("uTau"))),
    z0_(readScalar(dict.lookup("z0"))),
    n_(dict.lookup("n")),
    z_(dict.lookup("z"))

{

    if (mag(n_) < SMALL || mag(z_) < SMALL)
    {
        FatalErrorIn
        (
            "richardsHoxeyABLInletVelocityFvPatchVectorField"
            "("
                "const fvPatch&, "
                "const DimensionedField<vector, volMesh>&, "
                "onst dictionary&"
            ")"
        )
            << "magnitude of n or z must be greater than zero"
            << abort(FatalError);
    }

    evaluate();
}


richardsHoxeyABLInletVelocityFvPatchVectorField::richardsHoxeyABLInletVelocityFvPatchVectorField
(
    const richardsHoxeyABLInletVelocityFvPatchVectorField& fcvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(fcvpvf, iF),
    kappa_(fcvpvf.kappa_),
    uTau_(fcvpvf.uTau_),
    z0_(fcvpvf.z0_),
    n_(fcvpvf.n_),
    z_(fcvpvf.z_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void richardsHoxeyABLInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    n_ /= mag(n_);
    z_ /= mag(z_);

    const vectorField& cf = patch().Cf();

    const scalarField coord(cf & z_);

    scalarField Un(coord.size());

    forAll(coord, i)
    {
        scalar ref = (coord[i] + z0_)/z0_;

        Un[i] = (uTau_/kappa_) * Foam::log(ref);
    }

    vectorField::operator=(Un*n_);
}


// Write
void richardsHoxeyABLInletVelocityFvPatchVectorField::write(Ostream& os) const
{

    fvPatchVectorField::write(os);
    os.writeKeyword("kappa")
        << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("uTau")
        << uTau_ << token::END_STATEMENT << nl;
    os.writeKeyword("z0")
        << z0_ << token::END_STATEMENT << nl;
    os.writeKeyword("n")
        << n_ << token::END_STATEMENT << nl;
    os.writeKeyword("z")
        << z_ << token::END_STATEMENT << nl;
    writeEntry("value", os);

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, richardsHoxeyABLInletVelocityFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
