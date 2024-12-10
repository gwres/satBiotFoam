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

volScalarField PI = dims*pow(alpha,2)/(dims*lambda+2*Gv);
volScalarField Kdr = lambda+2*Gv/dims;

Info<<"Calculating Couplin Parameter for Fixed-Stress Split in "<<dims<<" dimensions."<<endl;  
if (couplingParameter.match("classical"))
{
	Info<<"Using classical Fixed-Stress coupling parameter"<<endl;
	PI = pow(alpha,2)/Kdr;
}
else if (couplingParameter.match("lambda"))
{
	Info<<"Using lambda Fixed-Stress coupling parameter"<<endl;
	PI = pow(alpha,2)/Kdr;
}
else if (couplingParameter.match("optimal"))
{
	Info<<"Using optimal Fixed-Stress coupling parameter"<<endl;
	PI = pow(alpha,2)/Kdr;
}
    
   
