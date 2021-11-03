/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v2

\*---------------------------------------------------------------------------*/

#include "DAMotionTranslationCoupled.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAMotionTranslationCoupled, 0);
addToRunTimeSelectionTable(DAMotion, DAMotionTranslationCoupled, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAMotionTranslationCoupled::DAMotionTranslationCoupled(
    const dynamicFvMesh& mesh,
    const DAOption& daOption)
    : DAMotion(
        mesh,
        daOption)
{
    M_ = daOption_.getAllOptions().subDict("rigidBodyMotion").getScalar("mass");
    C_ = daOption_.getAllOptions().subDict("rigidBodyMotion").getScalar("damping");
    K_ = daOption_.getAllOptions().subDict("rigidBodyMotion").getScalar("stiffness");
    y0_ = daOption_.getAllOptions().subDict("rigidBodyMotion").lookupOrDefault<scalar>("y0", 0.0);
    V0_ = daOption_.getAllOptions().subDict("rigidBodyMotion").lookupOrDefault<scalar>("V0", 0.0);
    scalarList dirList;
    daOption_.getAllOptions().subDict("rigidBodyMotion").readEntry<scalarList>("direction", dirList);
    direction_[0] = dirList[0];
    direction_[1] = dirList[1];
    direction_[2] = dirList[2];
    daOption_.getAllOptions().subDict("rigidBodyMotion").readEntry<wordList>("patchNames", patchNames_);
}

void DAMotionTranslationCoupled::correct()
{
    volVectorField& cellDisp =
        const_cast<volVectorField&>(mesh_.thisDb().lookupObject<volVectorField>("cellDisplacement"));

    scalar dT = mesh_.time().deltaT().value();

    vector force = this->getForce(mesh_);
    scalar yForce = force & direction_;

    // Euler method to solve the mass-spring-damper model
    y_ = y0_ + dT * V0_;
    V_ = V0_ + dT * (yForce - C_ * V0_ - K_ * y0_) / M_;

    y0_ = y_;
    V0_ = V_;

    forAll(patchNames_, idxI)
    {
        const word& patchName = patchNames_[idxI];
        label patchI = mesh_.boundaryMesh().findPatchID(patchName);

        forAll(cellDisp.boundaryField()[patchI], faceI)
        {
            cellDisp.boundaryFieldRef()[patchI][faceI] = y_ * direction_;
        }
    }

    // print information
    Info << "yForce: " << yForce << "  y: " << y_ << "  V: " << V_ << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //