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

Description
    Mesh motion solver for a polyMesh.  Based on solving the
    vertex-based totalLaplace motion equation.  The boundary motion is set as a
    boundary condition on the motion velocity variable motionU.

\*---------------------------------------------------------------------------*/

#include "totalLaplaceTetMotionSolver.H"
#include "motionDiff.H"
#include "addToRunTimeSelectionTable.H"
#include "tetFem.H"
#include "elementFields.H"
#include "tetFec.H"
#include "fixedValueTetPolyPatchFields.H"
#include "tetPolyPatchInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(totalLaplaceTetMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        tetMotionSolver,
        totalLaplaceTetMotionSolver,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::totalLaplaceTetMotionSolver::totalLaplaceTetMotionSolver
(
    const polyMesh& mesh,
    Istream&
)
:
    tetMotionSolver(mesh),
    fvMesh_(refCast<const fvMesh>(mesh)),
    diffusionPtr_(motionDiff::New(*this).ptr()),
    firstMotion_(true),
    solverPerf_()
{
    frozen_ = Switch(lookup("frozenDiffusion"));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::totalLaplaceTetMotionSolver::
~totalLaplaceTetMotionSolver()
{
    deleteDemandDrivenData(diffusionPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::totalLaplaceTetMotionSolver::solve()
{
    pointIOField refAllPoints
    (
        IOobject
        (
            "points",
            fvMesh_.time().constant(),
            polyMesh::meshSubDir,
            fvMesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    vectorField oldAllPoints = fvMesh_.allPoints();
  
    // Transform fixedValue bc from incremental to total
    forAll(motionU().boundaryField(), patchI)
    {
        if
        (
            isA<fixedValueTetPolyPatchVectorField>
            (
                motionU().boundaryField()[patchI]
            )
        )
        {
            fixedValueTetPolyPatchVectorField& patchMotionU =
                refCast<fixedValueTetPolyPatchVectorField>
                (
                    motionU().boundaryField()[patchI]
                );

            const labelList& meshPoints = 
                fvMesh_.boundaryMesh()[patchI].meshPoints();

            vectorField refPointU(meshPoints.size(), vector::zero);
            
            forAll (refPointU, pointI)
            {
                refPointU[pointI] =
                (
                    oldAllPoints[meshPoints[pointI]] 
                  - refAllPoints[meshPoints[pointI]]
                );
            }
            refPointU /= fvMesh_.time().value();

            tetPolyPatchInterpolation tppiAPatch
            (
                refCast<const faceTetPolyPatch>
                (
                    patchMotionU.patch()
                )
            );
                
            patchMotionU == patchMotionU
              + tppiAPatch.pointToPointInterpolate(refPointU);
        }
    }

    fvMesh& mesh = const_cast<fvMesh&>(fvMesh_);
    mesh.movePoints(refAllPoints);

    // Solve for mesh motion

    if (!frozen_ && !firstMotion_)
    {
        Info << "Correct mesh motion diffusion field." << endl;

        diffusionPtr_->correct();
    }

    tetFemVectorMatrix motionEqn
    (
        tetFem::laplacian
        (
            diffusion().motionGamma(),
            motionU()
        )
    );

    // Apply motion constraints
    applyConstraints(motionEqn);

    // Solve the motion equation
    if (firstMotion_)
    {
        firstMotion_ = false;

        // In the first solution, solve the motion twice to avoid relative
        // tolerance problem
        for (label i = 0; i < 2; i++)
        {
            solverPerf_ = motionEqn.solve();
        }
    }
    else
    {
        solverPerf_ = motionEqn.solve();
    }

    if (needTotDisplacement())
    {
        totDisplacement() += motionU()*tetMesh().time().deltaT();
    }

    mesh.movePoints(oldAllPoints);
}


void Foam::totalLaplaceTetMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    firstMotion_ = true;
    tetMotionSolver::updateMesh(mpm);
}


// ************************************************************************* //
