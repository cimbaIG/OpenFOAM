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

    Info << "Reading field measurementsU\n" << endl;
    const volVectorField measurementsU
    (
        IOobject
        (
            "measurementsU",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info << "Reading field measurementsR\n" << endl;
    const volTensorField measurementsR
    (
        IOobject
        (
            "measurementsR",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info << "Reading wind source term" << endl;
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
        //dimensionedVector("WS", dimLength / sqr(dimTime), vector::zero ),
        //zeroGradientFvPatchVectorField::typeName
    );

        volTensorField Rxy
	(
            IOobject
            (
		"Rxy",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
            mesh
        );

        volScalarField stressResidual
	(
            IOobject
            (
		"stressResidual",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
            mesh
        );

        /*volVectorField WScorrection
	(
            IOobject
            (
		"WScorrection",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
            mesh,
            dimensionedVector("WScorrection", dimLength / sqr(dimTime), vector::zero ),
            zeroGradientFvPatchVectorField::typeName
        );*/

        /*volVectorField diffR
	(
            IOobject
            (
		"diffR",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
            mesh,
            dimensionedVector("diffR", sqr(dimLength) / sqr(dimTime), vector::zero ),
            zeroGradientFvPatchVectorField::typeName
        );*/

    IOdictionary ABLDict
    (
        IOobject
        (
        "ABLDict",
        runTime.system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
        )
    );

    label WStime = readLabel(ABLDict.subDict("correctWS").lookup("correctionTimeStep"));

    label nWScorrect = readLabel(ABLDict.subDict("correctWS").lookup("numberOfWScorrections"));

    label iterCounter = 1;

    label iterWScorrect = 1;

#   include "createPhi.H"

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, simple.dict(), pRefCell, pRefValue);
    mesh.schemesDict().setFluxRequired(p.name());

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );
