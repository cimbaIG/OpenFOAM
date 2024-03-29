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

\*---------------------------------------------------------------------------*/

#include "fluidSolidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidSolidInterface> Foam::fluidSolidInterface::New
(
    dynamicFvMesh& fluidMesh,
    dynamicFvMesh& solidMesh
)
{
    word fsiTypeName;

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the fluid is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary fsiProperties
        (
            IOobject
            (
                "fsiProperties",
                fluidMesh.time().constant(),
                fluidMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        fsiProperties.lookup("fluidSolidInterface") >> fsiTypeName;
    }

    Info<< "Selecting fluidSolidInterface " << fsiTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(fsiTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "fluidSolidInterface::New(dynamicFvMesh&, dynamicFvMesh&)"
        )   << "Unknown fluidSolidInterface type " << fsiTypeName
            << endl << endl
            << "Valid fluidSolidInterface types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<fluidSolidInterface>(cstrIter()(fluidMesh, solidMesh));
}


// ************************************************************************* //
