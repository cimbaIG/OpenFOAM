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

#include "BlockFvmDivSigma.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "solidPolyMesh.H"
#include "divSigmaScheme.H"
#include "blockLaplacianScheme.H"
#include "blockLaplacianTransposeScheme.H"
#include "blockLaplacianTraceScheme.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace BlockFvm
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<BlockLduMatrix<vector> >
divSigma
(
    const solidPolyMesh& solidMesh,
    const surfaceScalarField& muf,
    const surfaceScalarField& lambdaf,
    GeometricField<vector, fvPatchField, volMesh>& U,
    Field<vector>& blockB
)
{
    return fv::divSigmaScheme::New
    (
        U.mesh(),
        U.mesh().schemesDict().divScheme("fvmDiv(sigma)")
    )().fvmDivSigma(solidMesh, muf, lambdaf, U, blockB);
}


tmp<BlockLduMatrix<vector> >
laplacian
(
    const solidPolyMesh& solidMesh,
    const surfaceScalarField& muf,
    GeometricField<vector, fvPatchField, volMesh>& U,
    Field<vector>& blockB
)
{
    return fv::blockLaplacian::New
    (
        U.mesh(),
        U.mesh().schemesDict().laplacianScheme
        (
            "fvmBlockLaplacian(" + U.name() + ')'
        )
    )().fvmBlockLaplacian(solidMesh, muf, U, blockB);
}


tmp<BlockLduMatrix<vector> >
laplacianTranspose
(
    const solidPolyMesh& solidMesh,
    const surfaceScalarField& muf,
    GeometricField<vector, fvPatchField, volMesh>& U,
    Field<vector>& blockB
)
{
    return fv::blockLaplacianTranspose::New
    (
        U.mesh(),
        U.mesh().schemesDict().laplacianScheme
        (
            "fvmBlockLaplacian(" + U.name() + ')'
        )
    )().fvmBlockLaplacianTranspose(solidMesh, muf, U, blockB);
}


tmp<BlockLduMatrix<vector> >
laplacianTrace
(
    const solidPolyMesh& solidMesh,
    const surfaceScalarField& lambdaf,
    GeometricField<vector, fvPatchField, volMesh>& U,
    Field<vector>& blockB
)
{
    return fv::blockLaplacianTrace::New
    (
        U.mesh(),
        U.mesh().schemesDict().laplacianScheme
        (
            "fvmBlockLaplacian(" + U.name() + ')'
        )
    )().fvmBlockLaplacianTrace(solidMesh, lambdaf, U, blockB);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace BlockFvm

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
