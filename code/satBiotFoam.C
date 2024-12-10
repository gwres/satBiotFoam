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

#include "fvCFD.H"
#include "IFstream.H"
#include "OFstream.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{	
	std::clock_t t_start= std::clock(); 									//Start Time

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
	#include "createFields.H"
	#include "createTimeControls.H"
	#include "readTimeControls.H"
    pisoControl piso(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nsatBiotFoam: Solving Biot Model Equations for fully saturated regime\n" << endl;
    
    
    /*------------------------------------Calculating no of cells in domain------------------------------------*/
	int NC = mesh.nCells();												
	Info<< "\nNumber of Cells in mesh = \n" << NC << endl;
	
	const faceList& faces = mesh.faces();                   // Face to node
	const labelList& faceOwner = mesh.faceOwner();          // Face to owner cell
	const labelList& faceNeighbour = mesh.faceNeighbour();  // Face to neighbour cell
	const volVectorField& C = mesh.C();         // Cell center coordinates
	Info << "Face Info \n" << faces << endl;
	Info << "Face Owner \n" << faceOwner << endl;
	Info << "Face Neighbour \n" << faceNeighbour << endl;
	surfaceTensorField KD_f = 0*fvc::interpolate(Ksat);
	
	forAll(faces, face)
	{
	  if (mesh.isInternalFace(face)) // Internal face is found
	  {
		  Info << "=================================" << endl;
		  Info << "Face No: " << face << endl;
		  Info << "Face ID: " << faces[face] << endl;
		  Info << "Owners ID: " << faceOwner[face] << endl;
		  Info << "Neighbours ID: " << faceNeighbour[face] << endl;
		  Info << "Coordiantes of Owner: " << C[faceOwner[face]] << ", Coordinate of Neighbour: " << C[faceNeighbour[face]] << endl;
		  df[face] = mag(C[faceOwner[face]]-C[faceNeighbour[face]]);
		  Info << "Distance at Face: " << df[face] << endl;
		  KD_f[face] = (Ksat[faceOwner[face]] - Ksat[faceNeighbour[face]])/
		  //~ Foam::pow(0.5*((Ksat[faceOwner[face]]&&Ksat[faceOwner[face]]) + (Ksat[faceNeighbour[face]]&&Ksat[faceNeighbour[face]])),0.5);
		  Foam::pow((Ksat[faceOwner[face]]&&Ksat[faceNeighbour[face]]),0.5);
	  }
	}
	Info << "df:\n" << df << endl;
	Info << "KD_f:\n" << KD_f << endl;
		
    /*----------------------------------Calculating domain and initial volume---------------------------------*/
	scalar Vol_D(0.0);	
	scalar V_0(0.0);	
	volScalarField epsilon_0 = fvc::div(u);											
	forAll(mesh.cells(),cellI)
	{
		Vol_D += mesh.V()[cellI];
		V_0 += (1+epsilon_0[cellI])*mesh.V()[cellI];
	}
	Info << "Initial Volume: " << V_0 << endl;
	Info << "Domain Volume: " << Vol_D << endl;
    
    /*------------------------------------Reading user defined Dictionaries------------------------------------*/
    Info << "\nReading User Defined Dictionaries\n" << endl;
    #include "read_Dictionaries.H"
    Info << "Time Factor: " << theta_f << endl;
    Info << "Lower, Upper and Max Limit for Picard iterations are: " << Npl << Npu << Npmax << endl;
    
    /*------------------------------------Calculating Fixed-Stress Parameter-----------------------------------*/
    #include "calculateFSParameter.H"
    
	double rmse = 1.0;
	double eps_max = err_lim.value();
		
    /*-----------------------Opening files for writing Picard iterations and Mass Balance Error-------------------------*/
	fileName file1 = ("Iterations.csv");
    OFstream OS1(file1);
	fileName file2 = ("massBalance.csv");
    OFstream OS2(file2);
    OS1 << "Time" << "," << "Delta T" << "," << "Iterations" << "," << "Stabilizing Counter" << endl;	
	OS2 << "Time" << "," <<  "Mass Balance Error" << endl;
	
	/*----------------------------------Initialisation of the previous time-step fields----------------------------------*/
	volScalarField p_0 = p;	
	volScalarField p_n = p;	
	volScalarField p_nm1 = p;	
	volVectorField u_n = u;
	volVectorField u_nm1 = u;
	volVectorField v_n = v;
	volVectorField HbyA_n = u;
	volVectorField HbyA = u;
	
    /*----------------------------------Initialisation of the previous iteration fields----------------------------------*/
	volScalarField p_m = p;	
	volVectorField u_m = u;
	volVectorField u_mm1 = u;
	
	//Initializing hydraulic conductivity
	phi = phi_0;
	Info << "phi: \n" << phi << endl;
	
	scalar patchFlux = 0;
	scalar patchinFlux = 0;
	scalar patchoutFlux = 0;
	scalar inFlux = 0;
	scalar outFlux = 0;
	scalar boundaryFlux = 0;
	scalar deltaVol_deform = 0;
	scalar deltaVol_press = 0;

    scalar oldTime = 0.0;			//Initialization of previous time variable
	label oldTimeIndex = 0;	
	label sc = 0;
	label TotalIter = 0;
		
    while (runTime.loop())
    {
		#include "readTimeControls.H"
        Info << nl << "Time = " << runTime.timeName() << endl;
        
        label itr=0; 
        rmse=1.0;
        
		gradU = fvc::grad(u_n);
		volScalarField div_un = fvc::div(u_n);	
		volScalarField div_unm1 = fvc::div(u_nm1);	
        
        while(itr<Npmax)
        {	
			itr = itr +1;
						
			#include "U_Eqn.H"
			
			#include "p_Eqn.H"
						
			#include "ConvergenceCriteria.H"
			
			#include "porosityModel.H"
			
			v = -1/(rho_f*ag)*(Kh & fvc::grad(p));
			
			Info << nl << "Time: " << runTime.timeName() << ", Iteration = " << itr << ", RMSE: " << rmse << endl;
			
			p_m = p;
			u_m = u;
			
			if(rmse<eps_max && itr>2)
			{
				break;
			}
		}
		
		if(itr==Npmax)
		{
			Info << "Convergence failed. Reiterating previous time step." << endl;
			Info << "Current Time: " << runTime.timeName() << "Current Delta T: " << runTime.deltaT().value() << endl; 
			
			scalar dt_new = max((runTime.deltaT().value())/theta_f,	minDeltaT);
			
			runTime.setDeltaT(dt_new, false);
			
			runTime.setTime			//Set Runtime to previous time-level
			(
				oldTime, oldTimeIndex
			);
			
			Info << "New Time: " << runTime.timeName() << "New Delta T: " << runTime.deltaT().value() << endl;
			
			p_n = p;
			u_n = u;
			
			sc = 0;
		}
		else
		{
			/*----------------------------------Time-Stepping Alogrithm----------------------------------*/
			#include "setDeltaT.H"			
			
			/*----------------------------Calculation of Mass Balance Error------------------------------*/
			#include "calculateMBE.H"
			
			
			/*--------------------------------Defining previous time-step--------------------------------*/
			oldTime = runTime.value();
			oldTimeIndex = runTime.timeIndex();
			
			
			/*----------------Writing Picard iterations and Mass Balance Error to files------------------*/
			OS1 << runTime.timeName() << "," << runTime.deltaT().value() << "," << itr << "," << sc << endl;
			OS2 << runTime.timeName() << "," <<  mag(1-mag(deltaVol_deform)/boundaryFlux)*100 << "," << mag(1-(mag(deltaVol_deform)+mag(deltaVol_press))/(boundaryFlux))*100 << endl;	
			
			p_nm1 = p_n;
			p_n = p;
			u_nm1 = u_n;
			u_n = u;
			v_n = v;
			HbyA_n = HbyA;
			
			TotalIter = TotalIter + itr;
			
			
			#include "calculateStress.H"
						
			runTime.write();
			
		}
		
    }
    
    std::clock_t t_end = std::clock();	//End Time
	
	//Writing the Simulation Time
	scalar t_total;
	t_total = (t_end - t_start)/ (double) CLOCKS_PER_SEC;
	Info<< "\nTotal Simulation Time: " << t_total << nl << endl;
	Info<< "\nTotal Iteration Count: " << TotalIter << nl << endl;
	
    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
