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

#ifndef pTractionBCFvPatchVectorField_H
#define pTractionBCFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

namespace Foam
{
class pTractionBCFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
	{	
		// Private data
		scalarField alpha_;
		vectorField traction_;
		    	
    	public:

			//Runtime type information
			TypeName("pTractionBC");

			//Constructors

				//Construct from patch and internal field
				pTractionBCFvPatchVectorField
				(
					const fvPatch&,
					const DimensionedField<vector, volMesh>&
				);

				//Construct from patch, internal field and dictionary
				pTractionBCFvPatchVectorField
				(
					const fvPatch&,
					const DimensionedField<vector, volMesh>&,
					const dictionary&
				);

				//Construct by mapping given TractionBCFvPatchVectorField onto a new patch
				pTractionBCFvPatchVectorField
				(
					const pTractionBCFvPatchVectorField&,
					const fvPatch&,
					const DimensionedField<vector, volMesh>&,
					const fvPatchFieldMapper&
				);

				//Construct as copy
				pTractionBCFvPatchVectorField
				(
					const pTractionBCFvPatchVectorField&
				);

				//Construct and return a clone
				virtual tmp<fvPatchVectorField> clone() const
				{
					return tmp<fvPatchVectorField>
					(
						new pTractionBCFvPatchVectorField(*this)
					);
				}

				//Construct as copy setting internal field reference
				pTractionBCFvPatchVectorField
				(
					const pTractionBCFvPatchVectorField&,
					const DimensionedField<vector, volMesh>&
				);

				//Construct and return a clone setting internal field reference
				virtual tmp<fvPatchVectorField> clone
				(
					const DimensionedField<vector, volMesh>& iF
				) const
				{
					return tmp<fvPatchVectorField>
					(
						new pTractionBCFvPatchVectorField(*this, iF)
					);
				}

			//Member functions

				//Update Coefficients
				virtual void updateCoeffs();

				//Write
				virtual void write(Ostream&) const;			
				
			//- access			
			virtual const scalarField& alpha() const
			{
					return alpha_;
			}
			
			virtual const vectorField& traction() const
			{
					return traction_;
			}
	};

} // End namespace Foam

#endif
