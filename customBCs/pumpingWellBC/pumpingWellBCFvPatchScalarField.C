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
    Foam::pumpingWellBCFvPatchScalarField

Group
    SatConsolFoam/customBCs

Description
    This boundary condition implements displacement 
    over right face for free flow zero traction condition.
\*-----------------------------------------------------------------------*/

#include "pumpingWellBCFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"


namespace Foam
{
	//Constructors
	
	pumpingWellBCFvPatchScalarField::pumpingWellBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(p, iF),
		Qw_(p.size(), 0)

		{}
	/*--------------------------------------------------------------------------*/
	pumpingWellBCFvPatchScalarField::pumpingWellBCFvPatchScalarField
	(
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const dictionary& dict
	)
	:
		fixedGradientFvPatchScalarField(p, iF),
		Qw_("Qw", dict, p.size())
		
		{
			if (dict.found("gradient"))
			{
				gradient() = scalarField("gradient", dict, p.size());
				fixedGradientFvPatchScalarField::updateCoeffs();
				fixedGradientFvPatchScalarField::evaluate();
			}
			else
			{
				fvPatchField<scalar>::operator=(patchInternalField());
				gradient() = 0;
			}
		}
	/*-------------------------------------------------------------------------------------*/	
	pumpingWellBCFvPatchScalarField::pumpingWellBCFvPatchScalarField
	(
		const pumpingWellBCFvPatchScalarField& ptf,
		const fvPatch& p,
		const DimensionedField<scalar, volMesh>& iF,
		const fvPatchFieldMapper& mapper
	)
	:
		fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
		Qw_(ptf.Qw_, mapper)
		
		{}
	/*--------------------------------------------------------------------------*/
	pumpingWellBCFvPatchScalarField::pumpingWellBCFvPatchScalarField
	(
		const pumpingWellBCFvPatchScalarField& ptf
	)
	:
		fixedGradientFvPatchScalarField(ptf),
		Qw_(ptf.Qw_)
		{}
	/*-------------------------------------------------------------------------------------*/
	pumpingWellBCFvPatchScalarField::pumpingWellBCFvPatchScalarField
	(
		const pumpingWellBCFvPatchScalarField& ptf,
		const DimensionedField<scalar, volMesh>& iF
	)
	:
		fixedGradientFvPatchScalarField(ptf, iF),
		Qw_(ptf.Qw_)
		{}

	//Member Functions

	void pumpingWellBCFvPatchScalarField::updateCoeffs()
	{
		if (updated())
		{
			return;
		}
		
		string patchName = patch().name();
		
		Info << "Entering pumpingWellBC - " << patchName << endl;
		
		const fvMesh& mesh(patch().boundaryMesh().mesh());
		
		label patchID = mesh.boundaryMesh().findPatchID(patchName);

		const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
		
		scalar patchArea = 0.0;
		forAll(cPatch, faceI)
		{
			patchArea += mesh.magSf().boundaryField()[patchID][faceI];
		}
        
        const fvPatchField<tensor>& Kh =
        patch().lookupPatchField<volTensorField, tensor>("Kh");
        
        const vectorField v_B = patch().nf()*Qw_/patchArea;
        
        gradient() = -(997*9.81)*patch().nf()&(v_B & inv(Kh));
		
		fixedGradientFvPatchScalarField::updateCoeffs();
	}

	void pumpingWellBCFvPatchScalarField::write(Ostream& os) const
	{
		fvPatchScalarField::write(os);
		Qw_.writeEntry("Qw",os);
		writeEntry("value", os);
	}

	makePatchTypeField(fvPatchScalarField, pumpingWellBCFvPatchScalarField);

} //End namespace Foam

