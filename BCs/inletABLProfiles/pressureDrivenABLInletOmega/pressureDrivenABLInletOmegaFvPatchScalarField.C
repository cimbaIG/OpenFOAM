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

#include "pressureDrivenABLInletOmegaFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pressureDrivenABLInletOmegaFvPatchScalarField::pressureDrivenABLInletOmegaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    uTau_(0),
    kappa_(0),
    hd_(0),
	Ck1_(0),
    Ck2_(0),
    Ck3_(0),
    Ck4_(0),
	Cu1_(0),
    Cu2_(0),
    Cu3_(0), 
    z_(0, 1, 0)
{}


pressureDrivenABLInletOmegaFvPatchScalarField::pressureDrivenABLInletOmegaFvPatchScalarField
(
    const pressureDrivenABLInletOmegaFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    uTau_(ptf.uTau_),
    kappa_(ptf.kappa_),
	hd_(ptf.hd_),
	Ck1_(ptf.Ck1_),
	Ck2_(ptf.Ck2_),
	Ck3_(ptf.Ck3_),
	Ck4_(ptf.Ck4_),
	Cu1_(ptf.Cu1_),
	Cu2_(ptf.Cu2_),
	Cu3_(ptf.Cu3_),
    z_(ptf.z_)
{}


pressureDrivenABLInletOmegaFvPatchScalarField::pressureDrivenABLInletOmegaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    uTau_(readScalar(dict.lookup("uTau"))),
    kappa_(readScalar(dict.lookup("kappa"))),
    hd_(readScalar(dict.lookup("hd"))),
    Ck1_(readScalar(dict.lookup("Ck1"))),
    Ck2_(readScalar(dict.lookup("Ck2"))),
    Ck3_(readScalar(dict.lookup("Ck3"))),
    Ck4_(readScalar(dict.lookup("Ck4"))),
    Cu1_(readScalar(dict.lookup("Cu1"))),
    Cu2_(readScalar(dict.lookup("Cu2"))),
    Cu3_(readScalar(dict.lookup("Cu3"))),
    z_(dict.lookup("z"))
{
    if (mag(z_) < SMALL)
    {
        FatalErrorIn("pressureDrivenABLInletOmegaFvPatchScalarField(dict)")
            << "z given with zero size not correct"
            << abort(FatalError);
    }

    evaluate();
}


pressureDrivenABLInletOmegaFvPatchScalarField::pressureDrivenABLInletOmegaFvPatchScalarField
(
    const pressureDrivenABLInletOmegaFvPatchScalarField& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    uTau_(fcvpvf.uTau_),
    kappa_(fcvpvf.kappa_),
    hd_(fcvpvf.hd_),
    Ck1_(fcvpvf.Ck1_),
    Ck2_(fcvpvf.Ck2_),
    Ck3_(fcvpvf.Ck3_),
    Ck4_(fcvpvf.Ck4_),
    Cu1_(fcvpvf.Cu1_),
    Cu2_(fcvpvf.Cu2_),
    Cu3_(fcvpvf.Cu3_),
    z_(fcvpvf.z_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void pressureDrivenABLInletOmegaFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    z_ /= mag(z_);

    const vectorField& cf = patch().Cf();

    const scalarField coord(cf & z_);

    scalarField Tke(coord.size());

    scalarField Omega(coord.size());

    //Info << coord << endl;
    
    // Get patch range and calculate domain height
    //boundBox bb(patch().patch().localPoints(), true);
    //scalar hd = bb.maxDim() - bb.minDim();

    forAll(coord, i)
    {

		scalar zRef = 1 - (coord[i]/hd_);
		scalar zRef_ = 1 - zRef;

        Tke[i] = sqr(uTau_) * 
      ( Ck1_
      + Ck2_ * Foam::pow(zRef,2)
	  + Ck3_ * Foam::pow(zRef,4)
	  + Ck4_ * Foam::pow(zRef,6) );

		Omega[i] = ((Tke[i]) / (kappa_ * uTau_ * coord[i])) * 
	  ( 1 + (1 + Cu1_)*zRef_ + ( 1 + Cu1_ + 2*Cu2_ )*Foam::pow(zRef_,2) +
      ( 1 + Cu1_ + 2*Cu2_ + 3*Cu3_)*Foam::pow(zRef_,3) );

    }

    scalarField::operator=(Omega);
}


// Write
void pressureDrivenABLInletOmegaFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("uTau")
        << uTau_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa")
        << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("hd")
        << hd_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ck1")
        << Ck1_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ck2")
        << Ck2_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ck3")
        << Ck3_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ck4")
        << Ck4_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cu1")
        << Cu1_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cu2")
        << Cu2_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cu3")
        << Cu3_ << token::END_STATEMENT << nl;
    os.writeKeyword("z")
        << z_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, pressureDrivenABLInletOmegaFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
