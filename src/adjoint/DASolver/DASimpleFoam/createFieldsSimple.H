Info << "Reading field p\n"
     << endl;
pPtr_.reset(
    new volScalarField(
        IOobject(
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE),
        mesh));
volScalarField& p = pPtr_();

Info << "Reading field U\n"
     << endl;
UPtr_.reset(
    new volVectorField(
        IOobject(
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE),
        mesh));
volVectorField& U = UPtr_();

/// user defined
Info << "Reading field uData\n"
     << endl;
uDataPtr_.reset(
    new volVectorField(
        IOobject(
            "uData",
            runTime.timeName(),
            target,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE),
        target));
volVectorField& uData = uDataPtr_();


#include "createPhiPython.H"

// create alpha porosity term
alphaPorosityPtr_.reset(
    new volScalarField(
        IOobject(
            "alphaPorosity",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("alphaPorosity", dimensionSet(0, 0, -1, 0, 0, 0, 0), pTraits<scalar>::zero),
        zeroGradientFvPatchField<scalar>::typeName));

// actuator related stuff
Info << "Creating source term. " << endl;
fvSourcePtr_.reset(
    new volVectorField(
        IOobject(
            "fvSource",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE),
        mesh,
        dimensionedVector("fvSource", dimensionSet(0, 1, -2, 0, 0, 0, 0), vector::zero),
        zeroGradientFvPatchField<vector>::typeName));

setRefCell(p, simple.dict(), pRefCell_, pRefValue_);
mesh.setFluxRequired(p.name());

laminarTransportPtr_.reset(
    new singlePhaseTransportModel(U, phi));
singlePhaseTransportModel& laminarTransport = laminarTransportPtr_();

turbulencePtr_.reset(
    incompressible::turbulenceModel::New(U, phi, laminarTransport));

// create MRF
MRFPtr_.reset(new IOMRFZoneListDF(mesh));

