    // Solve the density transport equation
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho) 
        + fvm::div(phi, rho)
        - fvm::laplacian(DK,rho)
     ==
        fvOptions(rho)
    );

    rhoEqn.relax();

    fvOptions.constrain(rhoEqn);

    rhoEqn.solve();

    fvOptions.correct(rho);

    deltaRho = rho/rhoRef*I;

