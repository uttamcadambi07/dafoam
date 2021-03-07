/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

    Description:
        Fixed traction boundary condition for the standard linear elastic, 
        fixed coefficient displacement equation.

\*---------------------------------------------------------------------------*/

#ifndef tractionDisplacementFvPatchVectorField_H
#define tractionDisplacementFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class tractionDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class tractionDisplacementFvPatchVectorField
    : public fixedGradientFvPatchVectorField
{

    // Private Data

    vectorField traction_;
    scalarField pressure_;

public:
    //- Runtime type information
    TypeName("tractionDisplacement");

    // Constructors

    //- Construct from patch and internal field
    tractionDisplacementFvPatchVectorField(
        const fvPatch&,
        const DimensionedField<vector, volMesh>&);

    //- Construct from patch, internal field and dictionary
    tractionDisplacementFvPatchVectorField(
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const dictionary&);

    //- Construct by mapping given
    //  tractionDisplacementFvPatchVectorField onto a new patch
    tractionDisplacementFvPatchVectorField(
        const tractionDisplacementFvPatchVectorField&,
        const fvPatch&,
        const DimensionedField<vector, volMesh>&,
        const fvPatchFieldMapper&);

    //- Construct as copy
    tractionDisplacementFvPatchVectorField(
        const tractionDisplacementFvPatchVectorField&);

    //- Construct and return a clone
    virtual tmp<fvPatchVectorField> clone() const
    {
        return tmp<fvPatchVectorField>(
            new tractionDisplacementFvPatchVectorField(*this));
    }

    //- Construct as copy setting internal field reference
    tractionDisplacementFvPatchVectorField(
        const tractionDisplacementFvPatchVectorField&,
        const DimensionedField<vector, volMesh>&);

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchVectorField> clone(
        const DimensionedField<vector, volMesh>& iF) const
    {
        return tmp<fvPatchVectorField>(
            new tractionDisplacementFvPatchVectorField(*this, iF));
    }

    // Member functions

    // Access

    virtual const vectorField& traction() const
    {
        return traction_;
    }

    virtual vectorField& traction()
    {
        return traction_;
    }

    virtual const scalarField& pressure() const
    {
        return pressure_;
    }

    virtual scalarField& pressure()
    {
        return pressure_;
    }

    // Mapping functions

    //- Map (and resize as needed) from self given a mapping object
    virtual void autoMap(
        const fvPatchFieldMapper&);

    //- Reverse map the given fvPatchField onto this fvPatchField
    virtual void rmap(
        const fvPatchVectorField&,
        const labelList&);

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Write
    virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //