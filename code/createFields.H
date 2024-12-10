Info<< "Reading pressure field\n" << endl;

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

Info<< "Reading displacement field\n" << endl;

volVectorField u
(
	IOobject
	(
		"u",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

Info<< "Reading porosity field\n" << endl;

volScalarField phi
(
	IOobject
	(
		"phi",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("phi", dimless, 1)
);

Info<< "Reading velocity field\n" << endl;

volVectorField v
(
	IOobject
	(
		"v",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

volTensorField sigma
(
	IOobject
	(
		"sigma",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedTensor("sigma", dimensionSet(1,-1,-2,0,0,0,0), tensor(0,0,0,0,0,0,0,0,0))
);

volTensorField Ksat
(
	IOobject
	(
		"Ksat",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

volTensorField Kh
(
	IOobject
	(
		"Kh",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	Ksat
);

volScalarField lambda
(
	IOobject
	(
		"lambda",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

volScalarField Gv
(
	IOobject
	(
		"Gv",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

volScalarField Mv_inv
(
	IOobject
	(
		"Mv_inv",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

volScalarField Sf
(
	IOobject
	(
		"Sf",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

volVectorField v_diff
(
	IOobject
	(
		"v_diff",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedVector("v_diff", dimensionSet(0,1,-1,0,0,0,0), vector(0,0,0))
);

volScalarField unit_p
(
	IOobject
	(
		"unit_p",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("unit_p", dimensionSet(1,-1,-2,0,0,0,0), 1.0)
);

volScalarField unit_umag
(
	IOobject
	(
		"unit_umag",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("unit_umag", dimensionSet(0,1,0,0,0,0,0), 1.0)
);
	
volScalarField divU
(
	IOobject
	(
		"divU",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("divU", dimensionSet(0,0,0,0,0,0,0), 0)
);

volTensorField gradU
(
	IOobject
	(
		"gradU",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedTensor("gradU", dimensionSet(0,0,0,0,0,0,0), tensor(0,0,0,0,0,0,0,0,0))
);

surfaceScalarField df
(
	IOobject
	(
		"df",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("df", dimensionSet(0,1,0,0,0,0,0), 0)
);

/*---------------------------------------------------------------------*/


