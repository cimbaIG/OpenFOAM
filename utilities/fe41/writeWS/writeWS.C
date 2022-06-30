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
#include "nearWallDist.H"

int main(int argc, char *argv[])
{
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
	#include "createFields.H"

	Info<< "\nStarting applying WS field...\n" << endl;

	const volVectorField& centers = mesh.C();
	
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

   Info<< "\nWriting wind source (WS) field.\n" << endl;

   volScalarField::GeometricBoundaryField d = nearWallDist(mesh);

   const fvPatchList& patches = mesh.boundary();

   forAll(patches, patchi)
   {
       const fvPatch& currPatch = patches[patchi];

       if (currPatch.isWall())
       {
           const scalarField& dist = d[patchi];

           const scalar minDist = min(dist);

           //- Wind source term field
           forAll(WS, cellI)
           {
               double zc = centers[cellI].y();

               if (zc <= 2*minDist)
               {
                   WS[cellI].x() = 0.;
                   WS[cellI].y() = 0.;
                   WS[cellI].z() = 0.;
               }
               else
               {
                   WS[cellI].x() = splineWS.evaluate(zc);
                   WS[cellI].y() = 0.;
                   WS[cellI].z() = 0.;
               }
               
               Info << "WS in cell " << cellI << " has value " << WS[cellI] << endl;
           }
           WS.write();
       }
   }

   /*forAll(WS, cellI)
   {
      double height = centers[cellI].y();

      WS[cellI].x() = splineWS.evaluate(height); 
      WS[cellI].y() = 0.;
      WS[cellI].z() = 0.;
   }
   WS.write();*/

	Info<< "End\n" << endl;

	return 0;
}

// ************************************************************************* //
