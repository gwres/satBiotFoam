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

#ifndef pumpingWellBCFvPatchScalarField_H
#define pumpingWellBCFvPatchScalarField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

namespace Foam
{
class pumpingWellBCFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
	{	
		// Private data
		scalarField Qw_;
		    	
    	public:

			//Runtime type information
			TypeName("pumpingWellBC");

			//Constructors

				//Construct from patch and internal field
				pumpingWellBCFvPatchScalarField
				(
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&
				);

				//Construct from patch, internal field and dictionary
				pumpingWellBCFvPatchScalarField
				(
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&,
					const dictionary&
				);

				//Construct by mapping given pumpingWellBCFvPatchScalarField onto a new patch
				pumpingWellBCFvPatchScalarField
				(
					const pumpingWellBCFvPatchScalarField&,
					const fvPatch&,
					const DimensionedField<scalar, volMesh>&,
					const fvPatchFieldMapper&
				);

				//Construct as copy
				pumpingWellBCFvPatchScalarField
				(
					const pumpingWellBCFvPatchScalarField&
				);

				//Construct and return a clone
				virtual tmp<fvPatchScalarField> clone() const
				{
					return tmp<fvPatchScalarField>
					(
						new pumpingWellBCFvPatchScalarField(*this)
					);
				}

				//Construct as copy setting internal field reference
				pumpingWellBCFvPatchScalarField
				(
					const pumpingWellBCFvPatchScalarField&,
					const DimensionedField<scalar, volMesh>&
				);

				//Construct and return a clone setting internal field reference
				virtual tmp<fvPatchScalarField> clone
				(
					const DimensionedField<scalar, volMesh>& iF
				) const
				{
					return tmp<fvPatchScalarField>
					(
						new pumpingWellBCFvPatchScalarField(*this, iF)
					);
				}

			//Member functions

				//Update Coefficients
				virtual void updateCoeffs();

				//Write
				virtual void write(Ostream&) const;			
				
			//- access			
			virtual const scalarField& Qw() const
			{
					return Qw_;
			}
	};

} // End namespace Foam

#endif
