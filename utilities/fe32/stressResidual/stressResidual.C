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
    stressResidual

Description
    Calculates and writes difference between measured and calculated stress.

    The -noWrite option just outputs the max/min values without writing
    the field.

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    IOobject Rxyheader
    (
        "Rxy",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject measurementsRheader
    (
        "measurementsR",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );


    if (Rxyheader.headerOk() && measurementsRheader.headerOk())
    {
        Info<< "    Reading Rxy" << endl;
        volTensorField Rxy(Rxyheader, mesh);

	Info<< "    Reading measurements R" << endl;
	volTensorField measurementsR(measurementsRheader, mesh);

        Info<< "    Calculating stress residual" << endl;
        volScalarField stressResidual
        (
            IOobject
            (
                "stressResidual",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            measurementsR.component(tensor::XY) - Rxy.component(tensor::YX)
        );

        Info<< "    stressResidual max/min : "
            << max(stressResidual).value() << " "
            << min(stressResidual).value() << endl;

        if (writeResults)
        {
            stressResidual.write();
        }
    }
    else
    {
        Info<< "    No Rxy or measurementsR" << endl;
    }

    Info<< "\nEnd\n" << endl;
}


// ************************************************************************* //
