    IOdictionary momSourceDict
    (
	IOobject
	(
	    "momSourceDict",
	    runTime.constant(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::NO_WRITE,
	    false
	)
    );
	
    const scalar tankL(readScalar(momSourceDict.lookup("tankLength")));
    const scalar dampX(readScalar(momSourceDict.lookup("dampZoneOutletStartPoint")));
    const scalar dampX0(readScalar(momSourceDict.lookup("dampZoneInletEndPoint")));
    const scalar dampR(readScalar(momSourceDict.lookup("dampCoffe")));
  
    volScalarField thetaField
    (       
        IOobject
        (
            "thetaField",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("thetaField", dimensionSet(0,0,-1,0,0,0,0), 0.0)
    );

    scalar exp_1 = Foam::exp(scalar(1)) - 1.0;
    thetaField = dimensionedScalar("thetaField", dimensionSet(0,0,-1,0,0,0,0), 0.0);

    forAll(mesh.cells(),cellI)
    { 
	if(mesh.C()[cellI].component(0) >= dampX)
	{ 
		scalar d0 = mesh.C()[cellI].component(0) - dampX;
		thetaField[cellI] = (Foam::exp(Foam::magSqr(d0 / (tankL - dampX))) - 1.0) * dampR / exp_1; 
            //thetaField[cellI] = Foam::sqrt( (d0 / (tankL - dampX)) * (d0 / (tankL - dampX))) * dampR; 
	}
    }
    forAll(mesh.cells(),cellI)
    { 
	if(mesh.C()[cellI].component(0) <= dampX0)
	{ 
		scalar d0 = mesh.C()[cellI].component(0) - dampX0;
		thetaField[cellI] = (Foam::exp(Foam::magSqr(d0 / dampX0)) - 1.0) * dampR / exp_1; 
            //thetaField[cellI] = Foam::sqrt( (d0 / (tankL - dampX)) * (d0 / (tankL - dampX))) * dampR;
	}
    }

