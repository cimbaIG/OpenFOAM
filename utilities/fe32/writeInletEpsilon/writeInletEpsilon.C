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
    writeInletEpsilon

Description
    Utility for defining velocity, turbulent kinetic energy and
	dissipation rate inlet boundary conditions. Utility reads 
	converged profiles and set them at the inlet of the computational
	domain.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "fixedGradientFvPatchFields.H"
#include "cubicSpline.H"

int main(int argc, char *argv[])
{
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
	#include "createFields.H"

	Info<< "\nStarting applying inlet boundary conditions...\n" << endl;

	const volVectorField& centers = mesh.C();
	volScalarField& eps = turbulence->epsilon()();

	//Load converged epsilon field
	CubicSpline splineEps;
	std::ifstream fileEps("inletEpsilon.dat");
        std::string line;
	while (std::getline(fileEps, line))
	{
	    std::stringstream iss(line);
	    double height, Eps;
	    iss >> height >> Eps;
	    splineEps.addPoint(height, Eps);
	}
	fileEps.close();
    
	// Get index of inlet patch
	label inletPatchID = mesh.boundaryMesh().findPatchID("inlet");
	
	//Get reference to boundary value
	const fvPatchVectorField& faceCenters = mesh.C().boundaryField()[inletPatchID];
	fvPatchScalarField& inletEpsBC = eps.boundaryField()[inletPatchID];

	Info<< "\nWriting epsilon inlet boundary profile.\n" << endl;

	//Look over patch faces
	forAll(inletEpsBC, faceI)
	{
		//Get coordinates of face centers
		const vector& c = faceCenters[faceI];
		//Write values of converged epsilon field to inletEps
		scalar inletEps = splineEps.evaluate(c.y());
		//Write BC
		inletEpsBC[faceI] = inletEps;
	}

	eps.write();

	Info<< "End\n" << endl;

	return 0;
}

// ************************************************************************* //
