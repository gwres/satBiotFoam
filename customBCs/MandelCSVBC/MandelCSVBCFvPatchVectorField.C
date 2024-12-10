/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
/*-----------------------------------------------------------------------*\
Class
    Foam::MandelCSVBCFvPatchVectorField

Group
    SatConsolFoam/customBCs

Description
    This boundary condition implements displacement 
    over right face for free flow zero traction condition.
\*-----------------------------------------------------------------------*/

#include "MandelCSVBCFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


namespace Foam
{
	//Constructors
	
	MandelCSVBCFvPatchVectorField::MandelCSVBCFvPatchVectorField
	(
		const fvPatch& p,
		const DimensionedField<vector, volMesh>& iF
	)
	:
		fixedGradientFvPatchVectorField(p, iF),
		dudn_()

		{}
	/*--------------------------------------------------------------------------*/
	MandelCSVBCFvPatchVectorField::MandelCSVBCFvPatchVectorField
	(
		const fvPatch& p,
		const DimensionedField<vector, volMesh>& iF,
		const dictionary& dict
	)
	:
		fixedGradientFvPatchVectorField(p, iF)
		
		{
			if (dict.found("dudn"))
			{
                dudn_ = Function1<scalar>::New("dudn", dict);
			}
			else
			{
				FatalIOErrorInFunction(dict)
				<< "Please supply dudn" << nl << exit(FatalIOError);
			}
			
			if (dict.found("gradient"))
			{
				gradient() = vectorField("gradient", dict, p.size());
				fixedGradientFvPatchVectorField::updateCoeffs();
				fixedGradientFvPatchVectorField::evaluate();
			}
			else
			{
				fvPatchField<vector>::operator=(patchInternalField());
				gradient() = vector(0,0,0);
			}
		}
	/*-------------------------------------------------------------------------------------*/	
	MandelCSVBCFvPatchVectorField::MandelCSVBCFvPatchVectorField
	(
		const MandelCSVBCFvPatchVectorField& ptf,
		const fvPatch& p,
		const DimensionedField<vector, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedGradientFvPatchVectorField(ptf, p, iF, mapper),
		//~ dudn_(ptf.dudn_.clone())
		dudn_()
		
		{}
	/*--------------------------------------------------------------------------*/
	MandelCSVBCFvPatchVectorField::MandelCSVBCFvPatchVectorField
	(
		const MandelCSVBCFvPatchVectorField& ptf
	)
	:
		fixedGradientFvPatchVectorField(ptf),
		//~ dudn_(ptf.dudn_.clone())
		dudn_()
		{}
	/*-------------------------------------------------------------------------------------*/
	MandelCSVBCFvPatchVectorField::MandelCSVBCFvPatchVectorField
	(
		const MandelCSVBCFvPatchVectorField& ptf,
		const DimensionedField<vector, volMesh>& iF
	)
	:
		fixedGradientFvPatchVectorField(ptf, iF),
		//~ dudn_(ptf.dudn_.clone())
		dudn_()
		{}

	//Member Functions

	void MandelCSVBCFvPatchVectorField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}
		
		const scalar t = db().time().timeOutputValue();
		
		Info << "Entering MandelCSVBC" << endl;
		Info << "dvdy entering top boundary: " << dudn_->value(t) << endl;
		//~ Info << "dvdy entering top boundary: " << dudn_ << endl;
		
		vectorField n = patch().nf();
        
		gradient() = (n*dudn_->value(t));
		//~ gradient() = (n);
		
		fixedGradientFvPatchVectorField::updateCoeffs();
	}

	void MandelCSVBCFvPatchVectorField::write(Ostream& os) const
	{
		fvPatchVectorField::write(os);
		if (dudn_.valid())
		{
			dudn_->writeData(os);
		}
		writeEntry("value", os);
	}

	makePatchTypeField(fvPatchVectorField, MandelCSVBCFvPatchVectorField);

} //End namespace Foam

