surfaceScalarField au
(
    IOobject
    (
        "phi",
        runTime().timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("1", dimless, 1)
);

dimensionedScalar alpha
(
    "alpha",
    dimless,
    1.0
);

{
    volVectorField HoA
    (
        IOobject
        (
            "HoA",
            runTime().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimVelocity, vector::zero),
        zeroGradientFvPatchVectorField::typeName
    );

    if
    (
        mesh.schemesDict().ddtScheme(U_.name())[0].wordToken()
     == fv::EulerDdtScheme<vector>::typeName
    )
    {
        au += fvc::interpolate(rAU)/runTime().deltaT();

        HoA = U_ + gradp_*rAU + (U_-U_.oldTime())*rAU/runTime().deltaT();
        if (1-alpha.value() > SMALL)
        {
            HoA -= (1-alpha)*(1+rAU/runTime().deltaT())
               *(U_.prevIter()-U_)/alpha;
        }
        HoA.correctBoundaryConditions();

        phi_ = (fvc::interpolate(HoA) & mesh.Sf());

        phi_ -= fvc::interpolate(rAU)*fvc::snGrad(p_)*mesh.magSf();
        phi_ += fvc::interpolate(rAU)*phi_.oldTime()/runTime().deltaT();
        if (1-alpha.value() > SMALL)
        {
            phi_ += (1-alpha)*au*phi_.prevIter()/alpha;
        }
        au = au + (1-alpha)*au/alpha;
        
        phi_ /= au;
    }
    else if
    (
        mesh.schemesDict().ddtScheme(U_.name())[0].wordToken()
     == fv::backwardDdtScheme<vector>::typeName
    )
    {
        scalar deltaT = runTime().deltaT().value();
        scalar deltaT0 = runTime().deltaT0().value();

        if (U_.oldTime().timeIndex() == U_.oldTime().oldTime().timeIndex())
        {
            deltaT0 = GREAT;
        }

        dimensionedScalar rDeltaT("rDeltaT", dimless/dimTime, 1.0/deltaT);
    
        scalar ct = 1 + deltaT/(deltaT + deltaT0); //+
        scalar c00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0)); //+
        scalar c0 = ct + c00; //-

        au = dimensionedScalar("1", dimless, 1)
          + fvc::interpolate(rAU)*ct*rDeltaT;

        HoA = U_ + gradp_*rAU
          + (ct*U_-c0*U_.oldTime()+c00*U_.oldTime().oldTime())*rAU*rDeltaT;
        if (1-alpha.value() > SMALL)
        {
            HoA -= (1-alpha)*(1+rAU*ct*rDeltaT)*(U_.prevIter()-U_)/alpha;
        }
        HoA.correctBoundaryConditions();

        phi_ = (fvc::interpolate(HoA) & mesh.Sf());

        phi_ -= fvc::interpolate(rAU)*fvc::snGrad(p_)*mesh.magSf();
        phi_ += fvc::interpolate(rAU)*c0*phi_.oldTime()*rDeltaT;
        phi_ -= fvc::interpolate(rAU)*c00*phi_.oldTime().oldTime()*rDeltaT;
        if (1-alpha.value() > SMALL)
        {
            phi_ += (1-alpha)*au*phi_.prevIter()/alpha;
        }

        au = au + (1-alpha)*au/alpha;
        
        phi_ /= au;
    }
    else if
    (
        mesh.schemesDict().ddtScheme(U_.name())[0].wordToken()
     == fv::steadyStateDdtScheme<vector>::typeName
    )
    {
        au = dimensionedScalar("1", dimless, 1);

        HoA = U_ + gradp_*rAU;
        if (1-alpha.value() > SMALL)
        {
            HoA -= (1-alpha)*(U_.prevIter()-U_)/alpha;
        }
        HoA.correctBoundaryConditions();

        phi_ = (fvc::interpolate(HoA) & mesh.Sf());
        phi_ -= fvc::interpolate(rAU)*fvc::snGrad(p_)*mesh.magSf();
        if (1-alpha.value() > SMALL)
        {
            phi_ += (1-alpha)*phi_.prevIter()/alpha;
        }
        au = 1 + (1-alpha)/alpha;
        
        phi_ /= au;
    }

    rAUf = fvc::interpolate(rAU)/au;
    
    forAll(phi_.boundaryField(), patchI)
    {
        if (U_.boundaryField()[patchI].fixesValue())
        {
            phi_.boundaryField()[patchI] =
            (
                U_.boundaryField()[patchI]
              & mesh.Sf().boundaryField()[patchI]
            );
        }
        else if (isA<slipFvPatchVectorField>(U_.boundaryField()[patchI]))
        {
            phi_.boundaryField()[patchI] = 0;
        }
        else if (isA<symmetryFvPatchVectorField>(U_.boundaryField()[patchI]))
        {
            phi_.boundaryField()[patchI] = 0;
        }
    }
}
