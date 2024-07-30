/*---------------------------------------------------------------------------*\
    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

    This file is a modified combination of the OpenFOAM's source codes
    src/TurbulenceModels/turbulenceModels/derivedFvPatchFields/wallFunctions/nutWallFunctions
    /nutUSpaldingWallFunction/nutUSpaldingWallFunctionFvPatchScalarField.C,
    src/TurbulenceModels/turbulenceModels/derivedFvPatchFields/wallFunctions/nutWallFunctions
    /nutUWallFunction/nutUWallFunctionFvPatchScalarField.C, and
    src/TurbulenceModels/turbulenceModels/derivedFvPatchFields/wallFunctions/nutWallFunctions
    /nutkWallFunction/nutWallFunctionFvPatchScalarFieldDF.C.

i.e. a simplified nut wall function usable with any turbulence model (only available for flat plates for now).

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

\*---------------------------------------------------------------------------*/
#include "nutWallFunctionFvPatchScalarFieldDF.H"
#include "DATurbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
addToRunTimeSelectionTable(nutWallFunctionFvPatchScalarField, nutWallFunctionFvPatchScalarFieldDF, dictionary);
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



tmp<scalarField> nutWallFunctionFvPatchScalarFieldDF::calcNut() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = fvPatchField<scalar>::db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
	(
            turbulenceModel::propertiesName,
            internalField().group()
	)
    );

    // Compute the gradient of the velocity field
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));
    const scalarField& y = turbModel.y()[patchi];

    fvPatchScalarField etaPatch_ = etaWallDV_.boundaryField()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    tmp<scalarField> tyPlus = CalcYPlus(magGradU);
    const scalarField& yPlus = tyPlus();

    tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
    scalarField& nutw = tnutw.ref();

    forAll(yPlus, facei)
    {
        // Set nutw based on yPlus threshold
        if (yPlus[facei] < 35)
        {
            nutw[facei] = 0;
        }
        else
        {
            nutw[facei] = nuw[facei] * ((yPlus[facei] * kappa_ / max(log(E_ * yPlus[facei]/etaPatch_[facei]),1)) - 1.0);
        }
    }
    return tnutw;
}

tmp<scalarField> nutWallFunctionFvPatchScalarFieldDF::CalcYPlus
(
    const scalarField& magGradU
) const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = fvPatchField<scalar>::db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
	(
            turbulenceModel::propertiesName,
            internalField().group()
	)
    );

    const scalarField& y = turbModel.y()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();
    scalar kappa=0.4;

    const scalarField uTau = kappa * y * magGradU;

    tmp<scalarField> tyPlus(new scalarField(patch().size()));
    scalarField& yPlus = tyPlus.ref();

    forAll(yPlus, facei)
    {
        yPlus[facei] = y[facei] * uTau[facei] / nuw[facei];
    }

    return tyPlus;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
nutWallFunctionFvPatchScalarFieldDF::nutWallFunctionFvPatchScalarFieldDF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName("turbulenceProperties", iF.group()),
            iF.time().constant(),
            iF.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    runTime_(iF.time()),
    mesh_(iF.mesh()),
    etaWallDV_(
        IOobject(
            "etaWallDV",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE),
        mesh_,
        dimensionedScalar("etaWallDV", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0),
        "zeroGradient")
{}

nutWallFunctionFvPatchScalarFieldDF::nutWallFunctionFvPatchScalarFieldDF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName("turbulenceProperties", iF.group()),
            iF.time().constant(),
            iF.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    runTime_(iF.time()),
    mesh_(iF.mesh()),
    etaWallDV_
    (
        IOobject
        (
            "etaWallDV",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("etaWallDV", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0),
        "zeroGradient"
    )
{}

nutWallFunctionFvPatchScalarFieldDF::nutWallFunctionFvPatchScalarFieldDF
(
    const nutWallFunctionFvPatchScalarFieldDF& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    IOdictionary(ptf),
    runTime_(ptf.runTime_),
    mesh_(ptf.mesh_),
    etaWallDV_(ptf.etaWallDV_)
{}

nutWallFunctionFvPatchScalarFieldDF::nutWallFunctionFvPatchScalarFieldDF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const word& propertiesName
)
:
    nutWallFunctionFvPatchScalarField(p,iF),
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName(propertiesName, alphaRhoPhi.group()),
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    runTime_(U.time()),
    mesh_(U.mesh()),
    etaWallDV_(
        IOobject(
            "etaWallDV",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE),
        mesh_,
        dimensionedScalar("etaWallDV", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0),
        "zeroGradient")
{}

nutWallFunctionFvPatchScalarFieldDF::nutWallFunctionFvPatchScalarFieldDF
(
    const nutWallFunctionFvPatchScalarFieldDF& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf),
    IOdictionary(wfpsf),
    runTime_(wfpsf.runTime_),
    mesh_(wfpsf.mesh_),
    etaWallDV_(wfpsf.etaWallDV_)
{}

nutWallFunctionFvPatchScalarFieldDF::nutWallFunctionFvPatchScalarFieldDF
(
    const nutWallFunctionFvPatchScalarFieldDF& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF),
    IOdictionary(wfpsf),
    runTime_(wfpsf.runTime_),
    mesh_(wfpsf.mesh_),
    etaWallDV_(wfpsf.etaWallDV_)
{}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<scalarField> nutWallFunctionFvPatchScalarFieldDF::yPlus() const
{
    const label patchi = patch().index();
    const turbulenceModel& turbModel = fvPatchField<scalar>::db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));
//    const scalarField& y = turbModel.y()[patchi];


    return CalcYPlus(magGradU);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(
    fvPatchScalarField,
    nutWallFunctionFvPatchScalarFieldDF);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
