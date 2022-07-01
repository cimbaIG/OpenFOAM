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
    frictionVelocity

Description
    Calculates and writes the friction velocity at the wall.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createMesh.H"

    IOdictionary ABLProperties
    (
        IOobject
        (
        "ABLProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        )
    );

    const scalar Cmu =
	readScalar(ABLProperties.lookup("Cmu"));

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;

        IOobject kheader
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check k exists
        if (kheader.headerOk())
        {
            mesh.readUpdate();

            Info<< "    Reading k" << endl;
            volScalarField k(kheader, mesh);

            Info<< "    Calculating friction velocity" << endl;

            volScalarField frictionVelocity
            (
                IOobject
                (
                    "frictionVelocity",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "frictionVelocity",
                    dimLength/dimTime,
                    0
                )
            );

            forAll(frictionVelocity.boundaryField(), patchi)
            {
                frictionVelocity.boundaryField()[patchi] =
                    Foam::pow(Cmu,0.25) * Foam::sqrt(k.boundaryField()[patchi]);
            } 

            frictionVelocity.write();
        }
        else
        {
            Info<< "    No k." << endl;
        }
    }

    Info<< "End" << endl;

    return 0;
}


// ************************************************************************* //
