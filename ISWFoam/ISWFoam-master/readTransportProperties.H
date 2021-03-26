singlePhaseTransportModel laminarTransport(U, phi);
dimensionedScalar rhoRef("rhoRef", dimDensity, laminarTransport);
dimensionedVector CoriolisVector("CoriolisVector", laminarTransport);
