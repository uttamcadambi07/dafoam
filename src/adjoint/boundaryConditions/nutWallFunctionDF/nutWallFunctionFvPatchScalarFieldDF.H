
/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

    This file is a modified combination of the OpenFOAM's source codes
    src/TurbulenceModels/turbulenceModels/derivedFvPatchFields/wallFunctions/
    nutWallFunctions/nutUSpaldingWallFunction/nutUSpaldingWallFunctionFvPatchScalarField.H,
    src/TurbulenceModels/turbulenceModels/derivedFvPatchFields/wallFunctions/
    nutWallFunctions/nutUWallFunction/nutUWallFunctionFvPatchScalarField.H, and
    src/TurbulenceModels/turbulenceModels/derivedFvPatchFields/wallFunctions/
    nutWallFunctions/nutkWallFunction/nutkWallFunctionFvPatchScalarField.H.

    OpenFOAM: The Open Source CFD Toolbox

    Copyright (C): 2011-2016 OpenFOAM Foundation

    OpenFOAM License:

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

    Description:

        NOTE: It's a simplified nut wall function usable with any turbulence model
	(only available for flat plates for now).

        This boundary condition provides a wall constraint on the turbulent viscosity,
	i.e. \c nut, based on velocity, i.e. \c U. Using the differential form of the
	log law to obtain utau (like in the nutUWallFunction), and a similar relation
	for nutw as in nutkWallFunction

    Usage
        Example of the boundary condition specification:
        \verbatim
        <patchName>
        {
            // Mandatory entries (unmodifiable)
            type		nutWallFunctionDF;

            // Optional entries (unmodifiable)
            kappa		0.4;
            E		        5;

            // Optional (inherited) entries
            ...
        }
        \endverbatim

        where the entries mean:
        \table
		Property	| Description			| Type   | Req'd  | Dflt
		type		| Type name: nutWallFunctionDF	| word   | yes    | -
		E		| Log Law Constant =exp(kappa*B)| scalar | no     | 5
           	kappa		| von Karman Constant		| scalar | no     | 0.4
        \endtable

        The inherited entries are elaborated in:
          - \link nutWallFunctionFvPatchScalarField.H \endlink

    Note
        - Suffers with highly curved surfaces or cases with strong PG since
        the implementation for buffer and subviscous layer has not been included
\*---------------------------------------------------------------------------*/

#ifndef nutWallFunctionFvPatchScalarFieldDF_H
#define nutWallFunctionFvPatchScalarFieldDF_H

#include "nutWallFunctionFvPatchScalarField.H"
#include "DAModel.H"
#include "DAOption.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nutWallFunctionFvPatchScalarFieldDF Declaration
\*---------------------------------------------------------------------------*/

class nutWallFunctionFvPatchScalarFieldDF
    : public nutWallFunctionFvPatchScalarField
{
protected:
// Protected Member Functions

    //- Calculate and return the yPlus values
    virtual tmp<scalarField> CalcYPlus(const scalarField& magGradU) const;

    //- Calculate the turbulence viscosity
    virtual tmp<scalarField> calcNut() const;

//    /// Foam::fvMesh object
//    const fvMesh& mesh_;
//
//    /// DAOptoin object
//    const DAOption& daOption_;
volScalarField etaWallDV_;
public:
    //- Runtime type information
    TypeName("nutWallFunctionDF");

// Constructors

        //- Construct from patch and internal field
        nutWallFunctionFvPatchScalarFieldDF
        (
            const fvPatch&,
	    const fvMesh&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutWallFunctionFvPatchScalarFieldDF
        (
            const fvPatch&,
            const fvMesh&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nutWallFunctionFvPatchScalarFieldDF
        //  onto a new patch
        nutWallFunctionFvPatchScalarFieldDF
        (
            const nutWallFunctionFvPatchScalarFieldDF&,
            const fvPatch&,
            const fvMesh&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutWallFunctionFvPatchScalarFieldDF
        (
            const nutWallFunctionFvPatchScalarFieldDF&,
            const fvMesh&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutWallFunctionFvPatchScalarFieldDF(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutWallFunctionFvPatchScalarFieldDF
        (
            const nutWallFunctionFvPatchScalarFieldDF&,
            const fvMesh&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
	    const fvMesh& mesh,
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new nutWallFunctionFvPatchScalarFieldDF(*this, mesh, iF)
            );
        }

// Member functions

    // Evaluation functions

	//- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
