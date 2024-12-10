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
    Foam::pTractionBCFvPatchVectorField

Group
    satBiotFoam/customBCs

Description
    This boundary condition implements displacement 
    over right face for free flow zero Traction condition.
\*-----------------------------------------------------------------------*/

#include "pTractionBCFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


namespace Foam
{
	//Constructors
	
	pTractionBCFvPatchVectorField::pTractionBCFvPatchVectorField
	(
		const fvPatch& p,
		const DimensionedField<vector, volMesh>& iF
	)
	:
		fixedGradientFvPatchVectorField(p, iF),
		alpha_(p.size(), 0),
		traction_(p.size(), vector(0,0,0))

		{}
	/*--------------------------------------------------------------------------*/
	pTractionBCFvPatchVectorField::pTractionBCFvPatchVectorField
	(
		const fvPatch& p,
		const DimensionedField<vector, volMesh>& iF,
		const dictionary& dict
	)
	:
		fixedGradientFvPatchVectorField(p, iF),
		alpha_("alpha", dict, p.size()),
		traction_("traction", dict, p.size())
		
		{
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
	pTractionBCFvPatchVectorField::pTractionBCFvPatchVectorField
	(
		const pTractionBCFvPatchVectorField& ptf,
		const fvPatch& p,
		const DimensionedField<vector, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedGradientFvPatchVectorField(ptf, p, iF, mapper),
		alpha_(ptf.alpha_, mapper),
		traction_(ptf.traction_, mapper)
		
		{}
	/*--------------------------------------------------------------------------*/
	pTractionBCFvPatchVectorField::pTractionBCFvPatchVectorField
	(
		const pTractionBCFvPatchVectorField& ptf
	)
	:
		fixedGradientFvPatchVectorField(ptf),
		alpha_(ptf.alpha_),
		traction_(ptf.traction_)
		{}
	/*-------------------------------------------------------------------------------------*/
	pTractionBCFvPatchVectorField::pTractionBCFvPatchVectorField
	(
		const pTractionBCFvPatchVectorField& ptf,
		const DimensionedField<vector, volMesh>& iF
	)
	:
		fixedGradientFvPatchVectorField(ptf, iF),
		alpha_(ptf.alpha_),
		traction_(ptf.traction_)
		{}

	//Member Functions

	void pTractionBCFvPatchVectorField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}
		
		Info << "Entering TractionBC" << endl;
		
		const fvPatchField<scalar>& p =
        patch().lookupPatchField<volScalarField, scalar>("p");
                
		//~ const fvPatchField<scalar>& ph =
        //~ patch().lookupPatchField<volScalarField, scalar>("ph");
                
        const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("lambda");
        
		const fvPatchField<scalar>& Gv =
        patch().lookupPatchField<volScalarField, scalar>("Gv");
        
		const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>("gradU");
        
        vectorField n = patch().nf();

        vectorField gradient_XYZ = ((traction_ + alpha_*p*n) - (n & (Gv*gradU.T() - (Gv + lambda)*gradU)) 
						- (n*tr(gradU)*lambda))/(2.0*Gv + lambda);
        //~ vectorField gradient_XYZ = ((traction_ + alpha_*(p-ph)*n) - (n & (Gv*gradU.T() - (Gv + lambda)*gradU)) 
						//~ - (n*tr(gradU)*lambda))/(2.0*Gv + lambda);
        
        Info << "Traction applied on patch: " << traction_ << endl;
        
		gradient() = gradient_XYZ;
		
		//~ gradient() =
		//~ (
		//~ traction_ + (alpha_*p)*n
		//~ - (n & (Gv*gradU.T() - (Gv + lambda)*gradU))
		//~ - n*tr(gradU)*lambda
		//~ )/(2.0*Gv + lambda);
		
		fixedGradientFvPatchVectorField::updateCoeffs();
	}

	void pTractionBCFvPatchVectorField::write(Ostream& os) const
	{
		fvPatchVectorField::write(os);
		alpha_.writeEntry("alpha",os);	
		traction_.writeEntry("traction",os);
		//~ os.writeEntry("alpha", alpha_);	
		//~ os.writeEntry("Traction", Traction_);	
		writeEntry("value", os);
	}

	makePatchTypeField(fvPatchVectorField, pTractionBCFvPatchVectorField);

} //End namespace Foam

