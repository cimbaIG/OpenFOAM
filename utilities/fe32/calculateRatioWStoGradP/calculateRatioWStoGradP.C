/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    calculateRatioWStoGradP

Description
    Utility for each written time step calulates the ratio between 
    wind source term WS ans pressure gradient gradP.

\*---------------------------------------------------------------------------*/

#   include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	#   include "setRootCase.H"
	#   include "createTime.H"
	#   include "createMesh.H"

	instantList timeDirs = timeSelector::select0(runTime, args);	

	forAll(timeDirs, timeI)
    	{
		runTime.setTime(timeDirs[timeI], timeI);

        	Info<< "Time = " << runTime.timeName() << "\n" << endl;

		//Reading field WS
		volVectorField WS
    		(
        		IOobject
        		(
            			"WS",
            			runTime.timeName(),
            			mesh,
            			IOobject::MUST_READ,
            			IOobject::AUTO_WRITE
        		),
        		mesh
    		);

		//Reading field p
		volScalarField p
    		(
        		IOobject
        		(
            			"p",
            			runTime.timeName(),
            			mesh,
            			IOobject::MUST_READ,
            			IOobject::AUTO_WRITE
        		),
        		mesh
    		);

                Info << "Calculating ratio between WSx and gradPx\n" << endl;

                volVectorField gradP = fvc::grad(p);

                Info << "maxWS = " << max(WS.component(vector::X)) << endl;
                Info << "minWS = " << min(WS.component(vector::X)) << endl;
              
                Info << "maxGradP_x = " << max(gradP.component(vector::X)) << endl;
                Info << "minGradP_x = " << min(gradP.component(vector::X)) << endl;

                //Creating and calculating field
                volScalarField ratioWSxtoGradPx
                (
                    IOobject
                    (
                        "ratioWSxtoGradPx",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    (WS.component(vector::X) / (gradP.component(vector::X) + stabilise(gradP.component(vector::X), 
                     dimensionedScalar("SMALL", dimensionSet(0,1,-2,0,0,0,0), SMALL))))
                );
              
                Info << "maxRatioWSxtoGradP_x = " << max(ratioWSxtoGradPx) << endl;
                Info << "minRatioWSxtoGradP_x = " << min(ratioWSxtoGradPx) << endl;
 
                ratioWSxtoGradPx.write();

                //Creating and calculating field
                volScalarField ratioMagWStoMagGradP
                (
                    IOobject
                    (
                        "ratioMagWStoMagGradP",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mag(WS) / ((mag(gradP) + stabilise(mag(gradP),
                    dimensionedScalar("SMALL", dimensionSet(0,1,-2,0,0,0,0), SMALL))))
                );

                Info << "maxRatiomagWStomagGradP = " << max(ratioMagWStoMagGradP) << endl;
                Info << "minRatiomagWStomagGradP = " << min(ratioMagWStoMagGradP) << endl;

                ratioMagWStoMagGradP.write();
	}

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
