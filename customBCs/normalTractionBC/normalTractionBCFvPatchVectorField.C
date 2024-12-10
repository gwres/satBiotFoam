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
    Foam::normalTractionBCFvPatchVectorField

Group
    SatConsolFoam/customBCs

Description
    This boundary condition implements displacement 
    over right face for free flow zero Traction condition.
\*-----------------------------------------------------------------------*/

#include "normalTractionBCFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


namespace Foam
{
	//Constructors
	
	normalTractionBCFvPatchVectorField::normalTractionBCFvPatchVectorField
	(
		const fvPatch& p,
		const DimensionedField<vector, volMesh>& iF
	)
	:
		fixedGradientFvPatchVectorField(p, iF),
		alpha_(p.size(), 0),
		traction_(p.size(), 0)

		{}
	/*--------------------------------------------------------------------------*/
	normalTractionBCFvPatchVectorField::normalTractionBCFvPatchVectorField
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
	normalTractionBCFvPatchVectorField::normalTractionBCFvPatchVectorField
	(
		const normalTractionBCFvPatchVectorField& ptf,
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
	normalTractionBCFvPatchVectorField::normalTractionBCFvPatchVectorField
	(
		const normalTractionBCFvPatchVectorField& ptf
	)
	:
		fixedGradientFvPatchVectorField(ptf),
		alpha_(ptf.alpha_),
		traction_(ptf.traction_)
		{}
	/*-------------------------------------------------------------------------------------*/
	normalTractionBCFvPatchVectorField::normalTractionBCFvPatchVectorField
	(
		const normalTractionBCFvPatchVectorField& ptf,
		const DimensionedField<vector, volMesh>& iF
	)
	:
		fixedGradientFvPatchVectorField(ptf, iF),
		alpha_(ptf.alpha_),
		traction_(ptf.traction_)
		{}

	//Member Functions

	void normalTractionBCFvPatchVectorField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}
		
		string patchName = patch().name();
		
		Info << "Entering normalTractionBC: " << patchName << endl;
		
		const fvMesh& mesh(patch().boundaryMesh().mesh());
		
		label patchID = mesh.boundaryMesh().findPatchID(patchName);

		const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
		
		scalar patchArea = 0.0;
		forAll(cPatch, faceI)
		{
			patchArea += mesh.magSf().boundaryField()[patchID][faceI];
			
		}
		
		Info << "Patch Area: " << patchArea << endl;
		
		const fvPatchField<vector>& u =
        patch().lookupPatchField<volVectorField, vector>("u");
        
		const fvPatchField<scalar>& p =
        patch().lookupPatchField<volScalarField, scalar>("p");
        
        const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("lambda");
        
		const fvPatchField<scalar>& Gv =
        patch().lookupPatchField<volScalarField, scalar>("Gv");
        
		const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>("gradU");
        
        vectorField n = patch().nf();
        //~ Info << (n & (n*traction_)) << endl;
        //~ Info << (n*p) << endl;
        
        //~ Info << "Normal Vector Field: " << n << endl;

        //~ vectorField gradient_XYZ = ((traction_ + alpha_*p)*n - (n & (Gv*gradU.T() - (Gv + lambda)*gradU)) 
						//~ - (n*divU*lambda))/(2.0*Gv + lambda);
        //~ vectorField gradient_XYZ = ((traction_ + alpha_*p)*n - (n & (Gv*gradU.T() - (Gv + lambda)*gradU)) 
						//~ - (n*tr(gradU)*lambda))/(2.0*Gv + lambda);
		//~ vectorField gradient_XYZ = (n*(traction_+alpha_*p-lambda*tr(gradU)) - (n&(Gv*gradU.T()-(lambda+Gv)*gradU)))/(lambda+2*Gv);
		//~ vectorField gradient_XYZ = (n*(TT+alpha_*p-lambda*tr(gradU)) - (n&(Gv*gradU.T()-(lambda+Gv)*gradU)))/(lambda+2*Gv);
		//~ vectorField gradient_XYZ = (n*(traction_+alpha_*p) - (n&(Gv*gradU.T())) - n*(lambda*divU))/(Gv);
		//~ vectorField gradient_XYZ = (n*(traction_+alpha_*p) - (n&(Gv*gradU.T())))/(Gv);
		//~ vectorField gradient_XYZ = (n*(traction_+alpha_*p) + (n&((lambda+Gv)*gradU)) - (n&(Gv*gradU.T())))/(lambda+2*Gv);
        
        //~ Info << gradient_XYZ << endl;
        //~ Info << "Traction applied on patch: " << traction_ << endl;
        
        gradient() =
		(
		(traction_ + alpha_*p)*n
		- (n & (Gv*gradU.T() - (Gv + lambda)*gradU))
		- n*tr(gradU)*lambda
		)/(2.0*Gv + lambda);
        
		//~ gradient() = gradient_XYZ;
		
		fixedGradientFvPatchVectorField::updateCoeffs();
	}

	void normalTractionBCFvPatchVectorField::write(Ostream& os) const
	{
		fvPatchVectorField::write(os);
		alpha_.writeEntry("alpha",os);	
		traction_.writeEntry("traction",os);
		//~ os.writeEntry("alpha", alpha_);	
		//~ os.writeEntry("Traction", Traction_);	
		writeEntry("value", os);
	}

	makePatchTypeField(fvPatchVectorField, normalTractionBCFvPatchVectorField);

} //End namespace Foam

