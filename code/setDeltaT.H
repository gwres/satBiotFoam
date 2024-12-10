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
if (adjustTimeStep)
{	
	if(itr<=Npl)
	{
		if(sc<Npstabmax)
		{
			sc += 1;
		}
		else
		{
			scalar dt_new = min(theta_f*runTime.deltaT().value(), maxDeltaT);
			runTime.setDeltaT(dt_new, false);
			sc = 0;
			Info << "New Delta T: " << runTime.deltaT().value() << endl;
		}
	}
	else if(itr>Npl && itr<=Npu)
	{
		sc=0;
	}
	else
	{
		scalar dt_new = max((1.0/theta_f)*runTime.deltaT().value(),minDeltaT);
		runTime.setDeltaT(dt_new, false);
		//~ if(itr==Npmax)
		//~ {
			//~ runTime.setTime			//Set Runtime to previous time-level
			//~ (
				//~ oldTime, oldTimeIndex
			//~ );
		//~ }
		sc = 0;
		Info << "New Delta T: " << runTime.deltaT().value() << endl;
	}
}
