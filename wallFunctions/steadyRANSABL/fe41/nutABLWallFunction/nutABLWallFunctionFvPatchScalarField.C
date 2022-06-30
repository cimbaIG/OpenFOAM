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

\*---------------------------------------------------------------------------*/

#include "nutABLWallFunctionFvPatchScalarField.H"
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

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> nutABLWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalarField& y = turbModel.y()[patchI];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const scalarField& nuw = turbModel.nu().boundaryField()[patchI];

    const scalar Cmu25 = pow025(Cmu_);

    tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
    scalarField& nutw = tnutw();

    // Get face cells
    const unallocLabelList& fc = patch().faceCells();

    forAll(nutw, faceI)
    {
        const label faceCellI = fc[faceI];

	const scalar uStar = Cmu25*sqrt(k[faceCellI]);
        const scalar yPlus = uStar*y[faceI]/nuw[faceI];
	const scalar y0 = z0_[faceI];

        nutw[faceI] = nuw[faceI]*(yPlus*kappa_/log((y[faceI] + y0)/y0) - 1.0);
    }

    return tnutw;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutABLWallFunctionFvPatchScalarField::nutABLWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    z0_(p.size(), 0.0)
{}


nutABLWallFunctionFvPatchScalarField::nutABLWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    z0_("z0", dict, p.size())
{}


nutABLWallFunctionFvPatchScalarField::nutABLWallFunctionFvPatchScalarField
(
    const nutABLWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    z0_(ptf.z0_, mapper)
{}


nutABLWallFunctionFvPatchScalarField::nutABLWallFunctionFvPatchScalarField
(
    const nutABLWallFunctionFvPatchScalarField& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf),
    z0_(wfpsf.z0_)
{}


nutABLWallFunctionFvPatchScalarField::nutABLWallFunctionFvPatchScalarField
(
    const nutABLWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF),
    z0_(wfpsf.z0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> nutABLWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchI = patch().index();

    const turbulenceModel& turbModel =
        db().lookupObject<turbulenceModel>("turbulenceModel");

    const scalarField& y = turbModel.y()[patchI];
    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    const scalarField kwc = k.boundaryField()[patchI].patchInternalField();
    const scalarField& nuw = turbModel.nu().boundaryField()[patchI];

    return pow025(Cmu_)*y*sqrt(kwc)/nuw;
}

void nutABLWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m 
)
{
    nutWallFunctionFvPatchScalarField::autoMap(m);
    z0_.autoMap(m);
}

void nutABLWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    nutWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const nutABLWallFunctionFvPatchScalarField& nrwfpsf =
	    refCast<const nutABLWallFunctionFvPatchScalarField>(ptf);

    z0_.rmap(nrwfpsf.z0_, addr);
}

void nutABLWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    z0_.writeEntry("z0", os);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, nutABLWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
