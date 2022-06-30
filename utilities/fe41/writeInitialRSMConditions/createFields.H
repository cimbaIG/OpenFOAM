    Info << "Reading field p\n" << endl;
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

    Info << "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    /*Info << "Reading field omega\n" << endl;
    volScalarField omega
    (
        IOobject
        (
            "omega",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );*/

    Info << "Reading wind source term\n" << endl;
    volVectorField WS
    (
        IOobject
        (
            "WS",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("WS", dimLength / sqr(dimTime), vector::zero ),
        zeroGradientFvPatchVectorField::typeName
    );

    Info << "Reading field measurementsU\n" << endl;
    volVectorField measurementsU
    (
        IOobject
        (
            "measurementsU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("measurementsU", dimLength / dimTime, vector::zero ),
        zeroGradientFvPatchVectorField::typeName
    );

    Info << "Reading field measurementsR\n" << endl;
    volTensorField measurementsR
    (
        IOobject
        (
            "measurementsR",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("measurementsR", sqr(dimLength) / sqr(dimTime), tensor::zero ),
        zeroGradientFvPatchTensorField::typeName
    );

    IOdictionary transportProperties
    (
      IOobject
      (
      "transportProperties",
      runTime.constant(),
      mesh,
      IOobject::MUST_READ,
      IOobject::NO_WRITE
      )
    );

    dimensionedScalar nu 
    (
	transportProperties.lookup("nu")
    );

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

#   include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );
