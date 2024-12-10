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

volScalarField C21 = Gv;
volScalarField C22 = lambda+Gv;
dimensionedScalar C23 = alpha;

fvVectorMatrix U_Eqn
(
    fvm::laplacian(lambda+2*Gv,u)
  ==
	fvc::laplacian(lambda+Gv,u)
  - (lambda+Gv)*(fvc::grad(fvc::div(u)))
);
U_Eqn.relax();
solve(U_Eqn == C23*fvc::grad(p));
u.correctBoundaryConditions();
