    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    UEqn().relax();

    eqnResidual = solve(UEqn() == (-fvc::grad(p) + WS)).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);
