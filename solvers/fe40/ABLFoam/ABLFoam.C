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
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    simpleControl simple(mesh);

#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

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

        volScalarField Cd
	(
            IOobject
            (
		"Cd",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ, //bolje NO_READ
		IOobject::AUTO_WRITE
	    ),
            //2*WS.component(vector::X)*Lref / nutGradU.component(tensor::YX),
            //0.25*WS.component(vector::X)*Lref / nutGradU.component(tensor::YX),
            WS.component(vector::X) * Zref / nutGradU.component(tensor::YX),
            zeroGradientFvPatchScalarField::typeName
        );
        Cd.write();

        // correction due to existence of the pressure gradient 
        WS -= fvc::grad(p);

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

        volScalarField WScorrection
	(
            IOobject
            (
		"WScorrection",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
            (Cd/Zref) * ( diffR ),
            zeroGradientFvPatchScalarField::typeName
        );
        WScorrection.write();

        //const volVectorField& centers = mesh.C();

        forAll(WS,cellI)
        {
 
            WS[cellI].x() += WScorrection[cellI];
            WS[cellI].y() = 0;
            WS[cellI].z() = 0;
          
        }
        WS.write();
        
        // correction due to incorrect shear stress over height
        //WS += fvc::div(- turbulence->nuEff() * (fvc::grad(measurementsU) - fvc::grad(U)) );
        //WS.write();
        //volSymmTensorField Rstress = turbulence->R();

	//volScalarField Rxy = turbulence->R()().component(symmTensor::XY);

	   /*forAll(Rstress, cellI)
	   {
		Rstress[cellI].xx() = 0;
		Rstress[cellI].yy() = 0;
		Rstress[cellI].zz() = 0;
		Rstress[cellI].xz() = 0;
		Rstress[cellI].yz() = 0;
	   }

           forAll(Rstress.boundaryField(),patchI)
           {
		forAll(Rstress.boundaryField()[patchI], faceI)
                {
                    Rstress.boundaryField()[patchI][faceI].xx() = 0;
	            Rstress.boundaryField()[patchI][faceI].yy() = 0;
		    Rstress.boundaryField()[patchI][faceI].zz() = 0;
		    Rstress.boundaryField()[patchI][faceI].xz() = 0;
		    Rstress.boundaryField()[patchI][faceI].yz() = 0;
                }
           }

        Rstress.write();*/

        Info << "A new value of the wind source term has been calculated!" << nl;

        volTensorField gradU = fvc::grad(U);

        volVectorField gradNutGradU = fvc::grad(turbulence->nut() * gradU.component(tensor::YX));
        gradU.write();
        nutGradU.write();
        gradNutGradU.write();

        /*forAll(WS, cellI)
        {
            const cell& c = mesh.cells()[cellI];

            scalar minY(VGREAT), maxY(-VGREAT);

            forAll(c, fI)
            {
                const face& f = mesh.faces()[c[fI]];

                forAll(f, pI)
                {
                    minY = min(minY, mesh.points()[f[pI]].y());
		    maxY = max(maxY, mesh.points()[f[pI]].y());
                }
            }

            WScorrection[cellI].x() = (1.0 / (maxY - minY)) * (measurementsR[cellI].xy() - nutGradU[cellI].yx());
            
            diffR[cellI].x() = (measurementsR[cellI].xy() - nutGradU[cellI].yx());

            WS[cellI].x() += WScorrection[cellI].x();
            WS[cellI].y() = 0;
            WS[cellI].z() = 0;

        } 
        WS.write();*/

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

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
