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
    writeInletBCs

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

	Info<< "\nStarting applying initial conditions...\n" << endl;

	const volVectorField& centers = mesh.C();
	volScalarField& k = turbulence->k()();
	volScalarField& nut = turbulence->nut()();
	volScalarField& eps = turbulence->epsilon()();

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

	//- Write velocity U
	forAll(U, cellI)
	{
	    double height = centers[cellI].y();

	    U[cellI].x() = splineU.evaluate(height);
	    U[cellI].y() = 0.;
	    U[cellI].z() = 0.;
	}
	U.write();

	//- Write TKE k
	forAll(k, cellI)
	{
	    double height = centers[cellI].y();

	    k[cellI] = splineK.evaluate(height);
	}
	k.write();

	//- Write turbulent viscosity nut
	forAll(nut, cellI)
	{
	    double height = centers[cellI].y();

	    nut[cellI] = splineNut.evaluate(height);
	}
	nut.write();

	//- Write turbulent dissipation epsilon
	forAll(eps, cellI)
	{
	    double height = centers[cellI].y();

	    eps[cellI] = splineEps.evaluate(height);
	}
	eps.write();

	Info<< "End\n" << endl;

	return 0;
}

// ************************************************************************* //
