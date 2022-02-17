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

Application
    ABLFoam

Description
    Solver for generation of Atmospheric Boundary Layer.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"
#       include "initConvergenceCheck.H"

        p.storePrevIter();

        // Pressure-velocity SIMPLE corrector
        {
#           include "UEqn.H"
#           include "pEqn.H"
        }

        turbulence->correct();

    if(iterCounter == WStime && iterWScorrect <= nWScorrect)
    { 
        scalar Zref = 0.202;

        volTensorField nutGradU
	(
            IOobject
            (
		"nutGradU",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
            -turbulence->nut() * fvc::grad(U)
        );

        volScalarField diffR
	(
            IOobject
            (
		"diffR",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ, //bolje NO_READ
		IOobject::AUTO_WRITE
	    ),
            measurementsR.component(tensor::XY) - nutGradU.component(tensor::YX),
            zeroGradientFvPatchScalarField::typeName
        );
        diffR.write();

        volScalarField ratioR
	(
            IOobject
            (
		"ratioR",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ, //bolje NO_READ
		IOobject::AUTO_WRITE
	    ),
            measurementsR.component(tensor::XY) / nutGradU.component(tensor::YX),
            zeroGradientFvPatchScalarField::typeName
        );
        ratioR.write();

    volScalarField WScorrected
	(
            IOobject
            (
		"WScorrected",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
            WS.component(vector::X) * ( ratioR ),
            zeroGradientFvPatchScalarField::typeName
        );
        WScorrected.write();

        forAll(WS,cellI)
        {
 
            WS[cellI].x() = WScorrected[cellI];
            WS[cellI].y() = 0;
            WS[cellI].z() = 0;
          
        }
        WS.write();

        Info << "A new value of the wind source term has been calculated!" << nl;

        volTensorField gradU = fvc::grad(U);

        volVectorField gradNutGradU = fvc::grad(turbulence->nut() * gradU.component(tensor::YX));
        gradU.write();
        nutGradU.write();
        gradNutGradU.write();

        iterCounter = 0;
        ++iterWScorrect;
    }

    ++iterCounter;
    //Info << "iterCounter: " << iterCounter << nl;

    Rxy = -turbulence->nut() * fvc::grad(U);
    stressResidual = measurementsR.component(tensor::XY) - Rxy.component(tensor::YX);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

#       include "convergenceCheck.H"

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
