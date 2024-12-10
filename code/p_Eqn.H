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
volScalarField C11 = (Mv_inv)+PI;
volScalarField C12 = PI;
dimensionedScalar C13 = alpha;
dimensionedScalar C14 = 1/(rho_f*ag);

volScalarField div_um = fvc::div(u_m);
p.oldTime()=p_n;

surfaceTensorField Kfc = fvc::interpolate(Kh,"interpolate_K");

while (piso.loop())
{
	volScalarField rAU(1.0/U_Eqn.A());
	HbyA = constrainHbyA(rAU*U_Eqn.H(), u, p);
	surfaceScalarField phiHbyA("phiHbyA",fvc::flux(HbyA));
	surfaceScalarField phiHbyA_n("phiHbyA",fvc::flux(HbyA_n));
	constrainPressure(p, u, phiHbyA, rAU);


	fvScalarMatrix p_Eqn
	(
		C11*fvm::ddt(p)
	  - C14*fvm::laplacian(Kh,p,"lapK")
	 
	 ==
		C12*(p_m - p_n)/(runTime.deltaT()) 
	  //~ - (C13/runTime.deltaT())*(fvc::div(u) - fvc::div(u_n))
	  - (C13/runTime.deltaT())*(fvc::div(phiHbyA) - fvc::div(phiHbyA_n))
	  - (C13/runTime.deltaT())*(fvm::laplacian(rAU,p) - fvc::laplacian(rAU,p_n))
	  + Sf
	);

	p_Eqn.relax();
	p_Eqn.solve();
	p.correctBoundaryConditions();

	u = HbyA + C13*rAU*fvc::grad(p);
	u.correctBoundaryConditions();
}

