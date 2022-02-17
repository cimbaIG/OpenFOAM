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
    writeTopBCs

Description
    Utility for defining turbulent kinetic energy and
	dissipation rate top boundary conditions. Utility reads 
	converged profiles and set them at the top of the computational
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

	Info<< "\nStarting applying top boundary conditions...\n" << endl;

	volScalarField& k = turbulence->k()();
	volScalarField& eps = turbulence->epsilon()();

	Info<< "\nReading converged turbulent kinetic energy data.\n" << endl;

	//Load converged k field
	CubicSpline splineK;
	std::ifstream fileK("convergedK.dat");
	std::string line;
	while (std::getline(fileK, line))
	{
	    std::stringstream iss(line);
	    double height, K;
	    iss >> height >> K;
	    splineK.addPoint(height, K);
	}
	fileK.close();  

	Info<< "\nReading converged turbulent viscosity data.\n" << endl;

	//Load converged nut field
	CubicSpline splineNut;
	std::ifstream fileNut("convergedNut.dat");
	while (std::getline(fileNut, line))
	{
	    std::stringstream iss(line);
	    double height, Nut;
	    iss >> height >> Nut;
	    splineNut.addPoint(height, Nut);
	}
	fileNut.close();  

	Info<< "\nReading converged dissipation rate data.\n" << endl;

	//Load converged epsilon field
	CubicSpline splineEps;
	std::ifstream fileEps("convergedEps.dat");
	while (std::getline(fileEps, line))
	{
	    std::stringstream iss(line);
	    double height, Eps;
	    iss >> height >> Eps;
	    splineEps.addPoint(height, Eps);
	}
	fileEps.close();

	Info<< "\nReading converged wind source data.\n" << endl;

	//Load converged WS field
	CubicSpline splineWS;
	std::ifstream fileWS("convergedWS.dat");
	while (std::getline(fileWS, line))
	{
	    std::stringstream iss(line);
	    double height, WS;
	    iss >> height >> WS;
	    splineWS.addPoint(height, WS);
	}
	fileWS.close();  

	// Get index of top patch
	label topPatchID = mesh.boundaryMesh().findPatchID("top");
	
	//Get reference to boundary value
	const fvPatchVectorField& faceCenters = mesh.C().boundaryField()[topPatchID];
	fvPatchScalarField& topKBC = k.boundaryField()[topPatchID];	
        fvPatchScalarField& topNUTBC = turbulence->nut()().boundaryField()[topPatchID];
	fvPatchScalarField& topEpsBC = eps.boundaryField()[topPatchID];
	fvPatchVectorField& topWSBC = WS.boundaryField()[topPatchID];

	Info<< "\nWriting turbulent kinetic energy at the top boundary.\n" << endl;

	//Look over patch faces
	forAll(topKBC, faceI)
	{
		//Get coordinates of face centers
		const vector& c = faceCenters[faceI];
		//Write values of converged k field to inletTKE
		scalar topTKE = splineK.evaluate(c.y());
		//Write BC
		topKBC[faceI] = topTKE;
	}

	Info<< "\nWriting turbulent viscosity at the top boundary.\n" << endl;

	//Look over patch faces
	forAll(topNUTBC, faceI)
	{
		//Get coordinates of face centers
		const vector& c = faceCenters[faceI];
		//Write values of converged k field to inletTKE
		scalar topNUT = splineNut.evaluate(c.y());
		//Write BC
		topNUTBC[faceI] = topNUT;
	}

	Info<< "\nWriting dissipation rate at the top boundary.\n" << endl;

	//Look over patch faces
	forAll(topEpsBC, faceI)
	{
		//Get coordinates of face centers
		const vector& c = faceCenters[faceI];
		//Write values of converged epsilon field to inletEps
		scalar topEps = splineEps.evaluate(c.y());
		//Write BC
		topEpsBC[faceI] = topEps;
	}

	Info<< "\nWriting wind source at the top boundary.\n" << endl;

	//Look over patch faces
	forAll(topWSBC, faceI)
	{
		//Get coordinates of face centers
		const vector& c = faceCenters[faceI];
		//Write values of converged WS field to vector inletWS
		vector topWS(splineWS.evaluate(c.y()), 0, 0);
		//Write BC
		topWSBC[faceI] = topWS;
	}

	k.write();
	eps.write();
	WS.write();
        turbulence->nut()().write();

	Info<< "End\n" << endl;

	return 0;
}

// ************************************************************************* //
