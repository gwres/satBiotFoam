/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    satBIOTFoam

Description
    Solves Consolidation and Pore Pressure Distribution

\*---------------------------------------------------------------------------*/
///////////////////////////////////////////////////////////////////////////////
IOdictionary domainProperties
(
	IOobject
	(
		"domainProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);
	
dimensionedScalar rho_f      ("rho_f",domainProperties);
dimensionedScalar ag		 ("ag",domainProperties);
word couplingParameter(domainProperties.lookup("couplingParameter"));
const scalar dims = domainProperties.lookupOrDefault<scalar>("dimensionality", 3.0);
const scalar phi_0 = domainProperties.lookupOrDefault<scalar>("phi_0", 0);
const scalar beta = domainProperties.lookupOrDefault<scalar>("beta", 4.4e-10);
const scalar beta1 = domainProperties.lookupOrDefault<scalar>("beta1", 1.0);
const scalar beta2 = domainProperties.lookupOrDefault<scalar>("beta2", 1.0);
dimensionedScalar alpha      ("alpha",domainProperties);
dimensionedScalar err_lim    ("err_lim",domainProperties);
///////////////////////////////////////////////////////////////////////////////
IOdictionary PicardProperties
(
	IOobject
	(
		"PicardProperties",
		runTime.constant(),
		mesh,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);
	
const label Npl = PicardProperties.lookupOrDefault<label>("Picard_lower_limit", 1);
const label Npu = PicardProperties.lookupOrDefault<label>("Picard_upper_limit", 1);
const label Npmax = PicardProperties.lookupOrDefault<label>("Max_Picard_iterations", 1);
const label Npstabmax = PicardProperties.lookupOrDefault<label>("Max_Stabilization_Counters", 1);
const scalar theta_f = PicardProperties.lookupOrDefault<scalar>("Time_factor", 1.0);
const scalar minDeltaT = PicardProperties.lookupOrDefault<scalar>("minDeltaT", 1.0);
///////////////////////////////////////////////////////////////////////////////
