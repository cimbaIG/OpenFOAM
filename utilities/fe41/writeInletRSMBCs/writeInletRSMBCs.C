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
    writeInletRSMBCs

Description
    Utility for defining velocity, turbulent kinetic energy,
	dissipation rate and Reynolds stress inlet boundary conditions. 
	Utility reads converged profiles and set them at the inlet of the 
	computational domain. This utility should be used with RSM
	turbulence models.

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

	volScalarField& k = turbulence->k()();
	volScalarField& eps = turbulence->epsilon()();
        volScalarField& nut = turbulence->nut()();
	volSymmTensorField& R = turbulence->R()();

	//- Load converged WS field
	CubicSpline splineWS;
	std::ifstream fileWS("convergedWS.dat");
	std::string line;
	while (std::getline(fileWS, line))
	{
	    std::stringstream iss(line);
	    double height, WS;
	    iss >> height >> WS;
	    splineWS.addPoint(height, WS);
	}
	fileWS.close();

	//Load converged U field
	CubicSpline splineU;
	std::ifstream fileU("convergedU.dat");
	while (std::getline(fileU, line))
	{
	    std::stringstream iss(line);
	    double height, velocity;
	    iss >> height >> velocity;
	    splineU.addPoint(height, velocity);
	}
	fileU.close();

	//Load converged k field
	CubicSpline splineK;
	std::ifstream fileK("convergedK.dat");
	while (std::getline(fileK, line))
	{
	    std::stringstream iss(line);
	    double height, K;
	    iss >> height >> K;
	    splineK.addPoint(height, K);
	}
	fileK.close();  

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

        //Load converged epsilon field
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

	//Load converged Rxx field
	CubicSpline splineRxx;
        std::ifstream fileRxx("convergedRxx.dat");
        while (std::getline(fileRxx, line))
        {
            std::stringstream iss(line);
            double height, Rxx;
            iss >> height >> Rxx;
            splineRxx.addPoint(height, Rxx);
        }
        fileRxx.close();
    
	//Load converged Rxy field
        CubicSpline splineRxy;
        std::ifstream fileRxy("convergedRxy.dat");
        while (std::getline(fileRxy, line))
        {
            std::stringstream iss(line);
            double height, Rxy;
            iss >> height >> Rxy;
            splineRxy.addPoint(height, Rxy);
        }
        fileRxy.close();

	//Load converged Ryy field
        CubicSpline splineRyy;
        std::ifstream fileRyy("convergedRyy.dat");
        while (std::getline(fileRyy, line))
        {
            std::stringstream iss(line);
            double height, Ryy;
            iss >> height >> Ryy;
            splineRyy.addPoint(height, Ryy);
        }
        fileRyy.close();

        //Load converged Rzz field
        CubicSpline splineRzz;
        std::ifstream fileRzz("convergedRzz.dat");
        while (std::getline(fileRzz, line))
        {
            std::stringstream iss(line);
            double height, Rzz;
            iss >> height >> Rzz;
            splineRzz.addPoint(height, Rzz);
        }
        fileRzz.close();

	// Get index of inlet patch
	label inletPatchID = mesh.boundaryMesh().findPatchID("inlet");
	
	//Get reference to boundary value
	const fvPatchVectorField& faceCenters = mesh.C().boundaryField()[inletPatchID];
	fvPatchVectorField& inletUBC = U.boundaryField()[inletPatchID];
	fvPatchScalarField& inletKBC = k.boundaryField()[inletPatchID];	
	fvPatchScalarField& inletEpsBC = eps.boundaryField()[inletPatchID];
    	fvPatchScalarField& inletNutBC = nut.boundaryField()[inletPatchID];
	fvPatchVectorField& inletWSBC = WS.boundaryField()[inletPatchID];
	fvPatchSymmTensorField& inletRBC = R.boundaryField()[inletPatchID];

	Info<< "\nWriting velocity inlet boundary profile.\n" << endl;

	//Look over patch faces
	forAll(inletUBC, faceI)
	{
		//Get coordinates of face centers
		const vector& c = faceCenters[faceI];
		//Write values of converged U field to vector inletVelocity
		//vector inletVelocity(splineU.evaluate(c.y()), 0, 0);
		vector inletVelocity(splineU.evaluate(c.y()), 0, 0);
        	//Write BC
		inletUBC[faceI] = inletVelocity;
	}

	Info<< "\nWriting turbulent kinetic energy inlet boundary profile.\n" << endl;

	//Look over patch faces
	forAll(inletKBC, faceI)
	{
		//Get coordinates of face centers
		const vector& c = faceCenters[faceI];
		//Write values of converged k field to inletTKE
		scalar inletTKE = splineK.evaluate(c.y());
		//Write BC
		inletKBC[faceI] = inletTKE;
	}

    /*Info<< "\nWriting specific dissipation rate inlet boundary profile.\n" << endl;
    //Look over patch faces
    forAll(inletOmegaBC, faceI)
    {
       //Get coordinates of face centers
       const vector& c = faceCenters[faceI];
       //Write values of converged omega field to inletOmega
       scalar inletOmega = splineOmega.evaluate(c.y());
       //Write BC
       inletOmegaBC[faceI] = inletOmega;
    }*/

	Info<< "\nWriting dissipation rate inlet boundary profile.\n" << endl;

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

    Info<< "\nWriting turb. viscosity inlet boundary profile.\n" << endl;

	//Look over patch faces
	forAll(inletNutBC, faceI)
	{
		//Get coordinates of face centers
		const vector& c = faceCenters[faceI];
		//Write values of converged epsilon field to inletEps
		scalar inletNut = splineNut.evaluate(c.y());
		//Write BC
		inletNutBC[faceI] = inletNut;
	}

	Info<< "\nWriting wind source inlet boundary profile.\n" << endl;

	//Look over patch faces
	forAll(inletWSBC, faceI)
	{
		//Get coordinates of face centers
		const vector& c = faceCenters[faceI];
		//Write values of converged WS field to vector inletWS
		vector inletWS(splineWS.evaluate(c.y()), 0, 0);
		//Write BC
		inletWSBC[faceI] = inletWS;
	}

        Info<< "\nWriting R inlet boundary profile.\n" << endl;

        //Look over patch faces
        forAll(inletRBC, faceI)
        {
            //Get coordinates of face centers
            const vector& c = faceCenters[faceI];
            //Write values of converged R field to symm tensor inletR
            symmTensor inletR(splineRxx.evaluate(c.y()), splineRxy.evaluate(c.y()), 0, splineRyy.evaluate(c.y()), 0, splineRzz.evaluate(c.y()));
            //Write BC
            inletRBC[faceI] = inletR;
        }

	U.write();
	k.write();
	eps.write();
    	nut.write();
	WS.write();
    	//omega.write();
	R.write();

	Info<< "End\n" << endl;

	return 0;
}

// ************************************************************************* //
