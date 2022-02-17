    // Solve the momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
      + (fCoriolis ^ U)
      ==
      - GradP
    );

    // Momentum solution
    solve(relax(UEqn()) == -fvc::grad(p));
