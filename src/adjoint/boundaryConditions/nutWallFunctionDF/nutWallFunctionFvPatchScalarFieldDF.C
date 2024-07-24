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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



tmp<scalarField> nutWallFunctionFvPatchScalarFieldDF::calcNut() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
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
    const fvMesh& mesh = turbModel.U().mesh();

    //const DATurbulenceModel& constDaTurb = mesh.thisDb().lookupObject<DATurbulenceModel>("DATurbulenceModel");
    //DATurbulenceModel& daTurb = const_cast<DATurbulenceModel&>(constDaTurb);

//    DATurbulenceModel& daTurb = const_cast<DATurbulenceModel&>(daModelPtr_->getDATurbulenceModel());

//    Info<< "Reading field eta\n" << endl; volScalarField eta ( IOobject ( "eta", mesh.time().timeName(), mesh,
//    IOobject::MUST_READ, IOobject::AUTO_WRITE ), mesh );

//    volScalarField teta= db().lookupObject<volScalarField>("eta");
//    volScalarField teta = daTurb.eta();
//    fvPatchScalarField eta_ = teta.boundaryField()[patchi];
    //largeClass& returnMe = tReturnMe();
    tmp<volScalarField> tetaWallDV_ = initializeEtaDV(mesh);
    const volScalarField& etaDV_ = tetaWallDV_();
    fvPatchScalarField etaPatch_ = etaDV_.boundaryField()[patchi];
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

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
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

    const scalarField uTau = sqrt(nuw * magGradU);

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
    nutWallFunctionFvPatchScalarField(p, iF)
{}

nutWallFunctionFvPatchScalarFieldDF::nutWallFunctionFvPatchScalarFieldDF
(
    const nutWallFunctionFvPatchScalarFieldDF& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)

{}


nutWallFunctionFvPatchScalarFieldDF::nutWallFunctionFvPatchScalarFieldDF
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict)

{ 
    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
	(
            turbulenceModel::propertiesName,
            internalField().group()
	)
    );
    const fvMesh& mesh = turbModel.U().mesh();
/*
    this->initializeEtaDV(mesh);//volScalarField& etaDV = etaDV_();
    */
}



nutWallFunctionFvPatchScalarFieldDF::nutWallFunctionFvPatchScalarFieldDF
(
    const nutWallFunctionFvPatchScalarFieldDF& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf)
{}


nutWallFunctionFvPatchScalarFieldDF::nutWallFunctionFvPatchScalarFieldDF
(
    const nutWallFunctionFvPatchScalarFieldDF& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF)

{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<scalarField> nutWallFunctionFvPatchScalarFieldDF::yPlus() const
{
    const label patchi = patch().index();
    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));
    const scalarField& y = turbModel.y()[patchi];


    return CalcYPlus(magGradU);
}
/*
tmp<fvMesh> nutWallFunctionFvPatchScalarFieldDF::getFvMesh() const
{
    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
	(
            turbulenceModel::propertiesName,
            internalField().group()
	)
    );
    const fvMesh& mesh = turbModel.U().mesh();
    return mesh;
    
}
*/

tmp<volScalarField> nutWallFunctionFvPatchScalarFieldDF::initializeEtaDV(const fvMesh& mesh) const
{
    // use .ref() when returning with the new keyword
    // add code to get the mesh. "New"
    tmp<volScalarField> tetaWallDV_(new volScalarField(
          IOobject(
              "etaWallDV",
              mesh.time().timeName(),
              mesh,
              IOobject::MUST_READ_IF_MODIFIED,
              IOobject::AUTO_WRITE),
          mesh,
          dimensionedScalar("etaWallDV", dimensionSet(0, 0, 0, 0, 0, 0, 0), 1.0)));
    //volScalarField& etaWallDV_ = tetaWallDV_();
    return tetaWallDV_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(
    fvPatchScalarField,
    nutWallFunctionFvPatchScalarFieldDF);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
