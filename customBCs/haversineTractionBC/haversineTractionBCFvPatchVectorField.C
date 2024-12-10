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
    Foam::haversineTractionBCFvPatchVectorField

Group
    SatConsolFoam/customBCs

Description
    This boundary condition implements displacement 
    over right face for free flow zero traction condition.
\*-----------------------------------------------------------------------*/

#include "haversineTractionBCFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


namespace Foam
{
	//Constructors
	
	haversineTractionBCFvPatchVectorField::haversineTractionBCFvPatchVectorField
	(
		const fvPatch& p,
		const DimensionedField<vector, volMesh>& iF
	)
	:
		fixedGradientFvPatchVectorField(p, iF),
		alpha_(p.size(), 0),
		traction_(p.size(), vector(0,0,0)),
		//~ alpha_(0),
		//~ traction_(vector(0,0,0)),
		timePeriod_(0.0)

		{}
	/*--------------------------------------------------------------------------*/
	haversineTractionBCFvPatchVectorField::haversineTractionBCFvPatchVectorField
	(
		const fvPatch& p,
		const DimensionedField<vector, volMesh>& iF,
		const dictionary& dict
	)
	:
		fixedGradientFvPatchVectorField(p, iF),
		alpha_("alpha", dict, p.size()),
		traction_("traction", dict, p.size()),
		//~ alpha_(scalarField(dict.lookup("alpha"))),
		//~ traction_(vectorField(dict.lookup("traction"))),
		timePeriod_(readScalar(dict.lookup("timePeriod")))
		
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
	haversineTractionBCFvPatchVectorField::haversineTractionBCFvPatchVectorField
	(
		const haversineTractionBCFvPatchVectorField& ptf,
		const fvPatch& p,
		const DimensionedField<vector, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedGradientFvPatchVectorField(ptf, p, iF, mapper),
		alpha_(ptf.alpha_, mapper),
		traction_(ptf.traction_, mapper),
		timePeriod_(ptf.timePeriod_)
		
		{}
	/*--------------------------------------------------------------------------*/
	haversineTractionBCFvPatchVectorField::haversineTractionBCFvPatchVectorField
	(
		const haversineTractionBCFvPatchVectorField& ptf
	)
	:
		fixedGradientFvPatchVectorField(ptf),
		alpha_(ptf.alpha_),
		traction_(ptf.traction_),
		timePeriod_(ptf.timePeriod_)
		{}
	/*-------------------------------------------------------------------------------------*/
	haversineTractionBCFvPatchVectorField::haversineTractionBCFvPatchVectorField
	(
		const haversineTractionBCFvPatchVectorField& ptf,
		const DimensionedField<vector, volMesh>& iF
	)
	:
		fixedGradientFvPatchVectorField(ptf, iF),
		alpha_(ptf.alpha_),
		traction_(ptf.traction_),
		timePeriod_(ptf.timePeriod_)
		{}

	//Member Functions

	void haversineTractionBCFvPatchVectorField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}
		
		Info << "Entering haversineTractionBC" << endl;
		
		const fvPatchField<scalar>& p =
        patch().lookupPatchField<volScalarField, scalar>("p");
        
		const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("lambda");
        
		const fvPatchField<scalar>& Gv =
        patch().lookupPatchField<volScalarField, scalar>("Gv");
        
		scalar t = this->db().time().value();
        
		const fvPatchField<tensor>& gradU =
        patch().lookupPatchField<volTensorField, tensor>("gradU");
        
        vectorField n = patch().nf();
        
        vectorField updated_traction = traction_;
        
        updated_traction = traction_*Foam::pow(Foam::sin(constant::mathematical::pi*t/timePeriod_),2);
        
        vectorField gradient_XYZ = ((updated_traction + alpha_*p*n) - (n & (Gv*gradU.T() - (Gv + lambda)*gradU))
						- (n*tr(gradU)*lambda))/(2.0*Gv + lambda);
        
        Info << "Traction applied on patch: " << updated_traction << endl;
        
		gradient() = gradient_XYZ;
		
		fixedGradientFvPatchVectorField::updateCoeffs();
	}

	void haversineTractionBCFvPatchVectorField::write(Ostream& os) const
	{
		fvPatchVectorField::write(os);
		alpha_.writeEntry("alpha",os);	
		traction_.writeEntry("traction",os);	
		//~ os.writeEntry("alpha", alpha_);	
		//~ os.writeEntry("traction", traction_);
		os.writeEntry("timePeriod", timePeriod_);	
		writeEntry("value", os);
	}

	makePatchTypeField(fvPatchVectorField, haversineTractionBCFvPatchVectorField);

} //End namespace Foam

