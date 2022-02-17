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

    scalar omega = readScalar(ABLDict.subDict("coriolisSetup").lookup("omega")); //omega = 0.086213
    scalar latitude = readScalar(ABLDict.subDict("coriolisSetup").lookup("latitude")); //latitude = 51
    scalar ug = readScalar(ABLDict.subDict("gradP").lookup("ug")); //ug = 15.5; gp = -2.07698
    scalar vg = readScalar(ABLDict.subDict("gradP").lookup("vg")); //ug = 0.;

    vector f_;
    f_.x() = 0.;
    f_.y() = 0.;
    f_.z() = 2.0*omega*Foam::sin(latitude*mathematicalConstant::pi/180.0);

    vector GradP_;
    GradP_.x() = f_.z()*vg;
    GradP_.y() = -f_.z()*ug;
    GradP_.z() = 0.;
    
    // Create Coriolis force vector
    Info << "Creating the Coriolis force vector, fCoriolis..." << endl;
    volVectorField fCoriolis
    (
        IOobject
        (
            "fCoriolis",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("fCoriolis",dimensionSet(0, 0, -1, 0, 0, 0, 0),f_)
    );  

    // Create pressure gradient force vector
    Info << "Creating the pressure gradient force vector, GradP..." << endl;
    volVectorField GradP
    (
        IOobject
        (
            "GradP",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("GradP",dimensionSet(0, 1, -2, 0, 0, 0, 0),GradP_)
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
