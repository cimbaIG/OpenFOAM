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
    writeGradU

Description
    Utility for calculating and defining velocity gradient value at the patch
	top.

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

	Info<< "\nReading converged velocity data.\n" << endl;

	//Load converged U field
	CubicSpline splineU;
	std::ifstream fileU("convergedU.dat");
	std::string line;
	while (std::getline(fileU, line))
	{
	    std::stringstream iss(line);
	    double height, velocity;
	    iss >> height >> velocity;
	    splineU.addPoint(height, velocity);
	}
	fileU.close();

	Info<< "\nCalculating velocity gradient.\n" << endl;

	// Get index of top patch
	label topPatchID = mesh.boundaryMesh().findPatchID("top");
	
	//Get reference to boundary value
	const fvPatchVectorField& faceCenters = mesh.C().boundaryField()[topPatchID];
	fvPatchVectorField& topU = U.boundaryField()[topPatchID];
	vectorField& UgradTop = refCast<fixedGradientFvPatchVectorField>(topU).gradient();
	
	//Look over patch faces
	forAll(topU, faceI)
	{
		//Get coordinates of face centers
		const vector& c = faceCenters[faceI];

		tensor grad(tensor::zero);
        grad.xy() = splineU.evaluate(c.y(), 1);
		//Info << "grad.xy() " << grad.xy() << nl;

        const vector gradUTop = grad & vector(0, 1, 0);
        //Info << "gradUTop " << gradUTop << endl;
        //Info << "grad " << grad << endl;

        UgradTop = gradUTop;
	}
	U.write();

	Info<< "End\n" << endl;

	return 0;
}

// ************************************************************************* //
