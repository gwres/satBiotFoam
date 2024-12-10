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
    satBiotFoam

Description
    Solves Consolidation and Pore Pressure Distribution

\*---------------------------------------------------------------------------*/


/*-------------------------------Calculation of Boundary Flux--------------------------------*/
patchinFlux       = 0;
patchoutFlux      = 0;
forAll (mesh.boundaryMesh(), patchI)
{
	patchFlux = 0;
	forAll(v.boundaryField()[patchI], faceI)
	{
		vector Area_cell = mesh.Sf().boundaryField()[patchI][faceI];
		vector vel_cell = 0.5*(v_n.boundaryField()[patchI][faceI]+v.boundaryField()[patchI][faceI]);
		patchFlux += (vel_cell & Area_cell)*runTime.deltaTValue();
	}
	if(patchFlux<0)
	{
		patchinFlux += patchFlux;
	}
	else
	{
		patchoutFlux += patchFlux;
	}
}

inFlux += patchinFlux;
outFlux += patchoutFlux;
boundaryFlux = inFlux + outFlux;


/*-------------------------------Calculation of volume change--------------------------------*/
scalar V_n(0.0);
scalar V_np1(0.0);
scalar dV_p_np1(0.0);
volScalarField epsilon_n = fvc::div(u_n);
volScalarField epsilon_np1 = fvc::div(u);

forAll(mesh.cells(),cellI)
{
	V_n += (1+epsilon_n[cellI])*mesh.V()[cellI];
	V_np1 += (1+epsilon_np1[cellI])*mesh.V()[cellI];
	
	dV_p_np1 += beta*phi[cellI]*mesh.V()[cellI]*(p[cellI] - p_n[cellI]);
}

deltaVol_deform += V_np1-V_n;
deltaVol_press += dV_p_np1;
Info << "Influx: " << inFlux << ", Outflux: " << outFlux << ", Flux passing the boundaries: " << boundaryFlux << endl;
Info << "Initial Volume, Final Volume, Volume Change, New Volume Change: " << V_0 << "," << V_np1 << "," << V_np1-V_0 << ", " << deltaVol_deform << ", " <<endl;
Info << "Old Mass Balance Error: " << mag(1+(deltaVol_deform+deltaVol_press)/boundaryFlux)*100 << "%" << endl;
Info << "New Mass Balance Error: " << mag(1-(mag(deltaVol_deform)+mag(deltaVol_press))/(boundaryFlux))*100 << "%" << endl;
    
   
