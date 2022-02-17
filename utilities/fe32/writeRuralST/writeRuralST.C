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
writeST

Description
Utility for setting WS field for rural terrain simulations.
The utility sets WS profile that is zero up to height h1,
and linear decreases from height h1 to the height of the domain.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "fixedGradientFvPatchFields.H"

#include "cubicSpline.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

   Info<< "\nStarting applying boundary conditions\n" << endl;

   // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

   const volVectorField& centers = mesh.C();
   
   //- Create velocity field from Power-Law velocity profile!
   forAll(U, cellI)
   {
      scalar height = centers[cellI].y();
      U[cellI].x() = Uref * Foam::pow(height/Zref,alpha);
      U[cellI].y() = 0.;
      U[cellI].z() = 0.;
   }

   // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

   //- Load Rxy measurements and create spline interpolation
   CubicSpline splineRxy;
   std::ifstream fileRxy("measurementsRxy.dat");
   std::string line;
   std::vector<scalar> height, Rxy, WsAvg;
   while (std::getline(fileRxy, line))
   {
      std::stringstream iss(line);
      double h, value;
      iss >> h >> value;
      splineRxy.addPoint(h, value); //ovo sam odkomentirao jer inače ne radi spline od podataka za Rxy
      height.push_back(h);
      Rxy.push_back(value);
   }
   fileRxy.close();

 //******************************************************************  
   std::ofstream outFileRxy("calculatedRxy.dat");
   for(label i=1;i<=1000;++i)
   {
       const scalar Rxy = splineRxy.evaluate(i * 0.001,0);
       outFileRxy << (i * 0.001) << " " << Rxy << nl;
   }
   outFileRxy.close();

   std::ofstream outFileRxyDerivative("calculatedRxyDerivative.dat");
   for(label i=1;i<=1000;++i)
   {
        const scalar Rxy = splineRxy.evaluate(i * 0.001, 1);
        outFileRxyDerivative << (i * 0.001) << " " << Rxy << nl;
   }
   outFileRxyDerivative.close();
//********************************************************************

//- Load velocity measurements and create spline interpolation
   CubicSpline splineVelocity;
   std::ifstream fileVelocity("measurementsU.dat");
   while (std::getline(fileVelocity, line))
   {
      std::stringstream iss(line);
      double height, velocity;
      iss >> height >> velocity;
      splineVelocity.addPoint(height, velocity);
   }
   fileVelocity.close();
  
   std::ofstream outFileVelocity("calculatedMeasurementsU.dat");
   for(label i=1;i<=1000;++i)
   {
        const scalar UX = splineVelocity.evaluate(i * 0.001);
        outFileVelocity << (i * 0.001) << " " << UX << nl;
   }
   outFileVelocity.close();
 
   //- Experimental Velocity field
   forAll(measurementsU, cellI)
   {
      double height = centers[cellI].y();
      measurementsU[cellI].x() = splineVelocity.evaluate(height);
      measurementsU[cellI].y() = 0.;
      measurementsU[cellI].z() = 0.;
   }
   measurementsU.write(); 

   //- Experimental R field
   forAll(measurementsR, cellI)
   {
      double height = centers[cellI].y();

      measurementsR[cellI].xy() = splineRxy.evaluate(height);
      measurementsR[cellI].xx() = 0.;
      measurementsR[cellI].yy() = 0.;
      measurementsR[cellI].zz() = 0.;
      measurementsR[cellI].yz() = 0.;
      measurementsR[cellI].xz() = 0.;
    }
   measurementsR.write();

   //- Wind source term field
   forAll(WS, cellI)
   {
      double height = centers[cellI].y();

      if (height <= 0.202) 
      {
          WS[cellI].x() = 0.; 
          WS[cellI].y() = 0.;
          WS[cellI].z() = 0.;
      }
      else
      {
          WS[cellI].x() = 0.9734;
          WS[cellI].y() = 0.;
          WS[cellI].z() = 0.;
      }
      Info << "WS in cell " << cellI << " has value " << WS[cellI] << endl;
   }
   WS.write();

   // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
   volScalarField& nut = turbulence->nut()();
   volScalarField& k = turbulence->k()();
   volScalarField& eps = turbulence->epsilon()();

   //- Load K measurements and create spline interpolation
   CubicSpline splineK;
   std::ifstream fileK("measurementsK.dat");
   while (std::getline(fileK, line))
   {
      std::stringstream iss(line);
      double height, value;
      iss >> height >> value;
      splineK.addPoint(height, value);
   }
   fileK.close();
   
   std::ofstream outFileK("calculatedK.dat");
   for(label i=1;i<=1000;++i)
   {
        const scalar kVal = splineK.evaluate(i * 0.001);
        outFileK << (i * 0.001) << " " << kVal << nl;
   }
   outFileK.close();

   forAll(U.boundaryField(), patchI)
   {

      if( mesh.boundary()[patchI].name() == "top" )
      {
         Info << "Patch name " << mesh.boundary()[patchI].name() << endl;
         Info << "Patch index " << patchI << endl;

         fvPatchVectorField& uTop = U.boundaryField()[patchI];
         vectorField& UgradTop = refCast<fixedGradientFvPatchVectorField>(uTop).gradient();
	 scalarField& kTop = k.boundaryField()[patchI];
         //scalarField& KgradTop = refCast<fixedGradientFvPatchScalarField>(k.boundaryField()[patchI]).gradient();
	 scalarField& nut_top = refCast<fixedValueFvPatchScalarField>(nut.boundaryField()[patchI]);
	 scalarField& eps_top = refCast<fixedValueFvPatchScalarField>(eps.boundaryField()[patchI]);
        
	 const vector& centre = mesh.C().boundaryField()[patchI][0];
        
         tensor grad(tensor::zero);
         //grad.xy() = ((Uref * alpha)/Foam::pow(Zref,alpha)) * Foam::pow(centre.y(),alpha-1.0);
         grad.xy() = splineVelocity.evaluate(centre.y(), 1);
         //Info << "UTop " << U[patchI].x() << nl;
	 Info << "grad.xy() " << grad.xy() << nl;
         const vector gradUTop = grad & vector(0, 1, 0);
         Info << "gradUTop " << gradUTop << endl;
         Info << "grad " << grad << endl;

         UgradTop = gradUTop;
         Info << "Uref " << Uref << endl;
         Info << "Zref " << Zref << endl;
         Info << "alpha " << alpha << endl;

         const scalar C_mu = 0.044;
         const scalar Nu = 1.511e-5;

         /*const scalar nu =
         readScalar
             (
                 transportProperties.lookup("nu")
             );*/

         const scalar nutTop = mag((splineRxy.evaluate(centre.y())/mag(gradUTop)) - Nu);
         nut_top = nutTop;
         Info << "RxyTop " << splineRxy.evaluate(centre.y()) << nl;
         Info << "nutTop " << nutTop << nl;

	 kTop = splineK.evaluate(centre.y());
         const scalar kGradTop = splineK.evaluate(centre.y(), 1);
         Info << "k " << splineK.evaluate(centre.y()) << nl;
         //KgradTop = kGradTop;
         Info << "Grad k " << kGradTop << endl;

         const scalar epsTop = mag(C_mu * Foam::pow(splineK.evaluate(centre.y(),0),2) / nutTop);
         eps_top = epsTop;
         Info << "epsTop " << epsTop << nl; 

// Check Rxy value on the top boundary patch     

         const scalar TauXY = -(nutTop) * mag(gradUTop);
         Info << "TauXY " << TauXY << nl;

      }
   }
   U.write();
   nut.write();
   k.write();
   eps.write();

   // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

   Info<< "End\n" << endl;
   return 0;
}


// ************************************************************************* //
