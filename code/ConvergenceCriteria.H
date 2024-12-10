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

double rmse_p = Foam::sqrt(sum(pow(((p - p_m)/(1e-6*unit_p+p_m)),2.0)).value()/NC);
Info << "RMSE_p: " << rmse_p << endl;
			
double rmse_u = Foam::sqrt(sum(pow(((mag(u) - mag(u_m))/(1e-6*unit_umag+mag(u_m))),2.0)).value()/NC);
Info << "RMSE_u: " << rmse_u << endl;

rmse = std::max({rmse_p,rmse_u});


