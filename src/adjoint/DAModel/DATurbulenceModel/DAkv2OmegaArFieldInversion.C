/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

    This file is modified from OpenFOAM's source code
    src/TurbulenceModels/turbulenceModels/RAS/kOmega/kOmega.C

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

#include "DAkv2OmegaArFieldInversion.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(DAkv2OmegaArFieldInversion, 0);
addToRunTimeSelectionTable(DATurbulenceModel, DAkv2OmegaArFieldInversion, dictionary);
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DAkv2OmegaArFieldInversion::DAkv2OmegaArFieldInversion
(
    const word& modelType,
    const fvMesh& mesh,
    const DAOption& daOption)
    : DATurbulenceModel(modelType, mesh, daOption),
    // kv2OmegaAr parameters
    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            this->coeffDict_,
            4.04
        )
    ),
    AS_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "AS",
            this->coeffDict_,
            2.12
        )
    ),
    Anu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Anu",
            this->coeffDict_,
            3.8
        )
    ),
    ABP_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ABP",
            this->coeffDict_,
            0.2
        )
    ),
    ANAT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ANAT",
            this->coeffDict_,
            200
        )
    ),
    ATS_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "ATS",
            this->coeffDict_,
            200
        )
    ),
    CBPcrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CBPcrit",
            this->coeffDict_,
            1.5
        )
    ),
    CNC_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CNC",
            this->coeffDict_,
            0.1
        )
    ),
    CNATcrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CNATcrit",
            this->coeffDict_,
            1450
        )
    ),
    CINT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CINT",
            this->coeffDict_,
            0.95
        )
    ),
    CTScrit_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CTScrit",
            this->coeffDict_,
            1000
        )
    ),
    CRNAT_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CRNAT",
            this->coeffDict_,
            0.02
        )
    ),
    C11_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C11",
            this->coeffDict_,
            3.4e-6
        )
    ),
    C12_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C12",
            this->coeffDict_,
            1.0e-10
        )
    ),
    CR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CR",
            this->coeffDict_,
            0.32
        )
    ),
    CalphaTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CalphaTheta",
            this->coeffDict_,
            0.035
        )
    ),
    CSS_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CSS",
            this->coeffDict_,
            3.0
        )
    ),
    Ctau1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Ctau1",
            this->coeffDict_,
            4360
        )
    ),
    Cw1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw1",
            this->coeffDict_,
            0.44
        )
    ),
    Cw2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            this->coeffDict_,
            0.92
        )
    ),
    CwR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CwR",
            this->coeffDict_,
            1.15
        )
    ),
    Clambda_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Clambda",
            this->coeffDict_,
            2.495
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    PrTheta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "PrTheta",
            this->coeffDict_,
            0.85
        )
    ),
    sigmaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaK",
            this->coeffDict_,
            1
        )
    ),
    sigmaW_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaW",
            this->coeffDict_,
            1.17
        )
    ),
    sigmaW2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaW2",
            this->coeffDict_,
            1.856
        )
    ),
    sigmaAr_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaAr",
            this->coeffDict_,
            10
        )
    ),
    CA1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CA1",
            this->coeffDict_,
            1.0
        )
    ),
    CA2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CA2",
            this->coeffDict_,
            1.0
        )
    ),
    // Augmented variables for adjoint residuals 
    omega_(const_cast<volScalarField&>(
        mesh_.thisDb().lookupObject<volScalarField>("omega"))),
    omegaRes_(
        IOobject
        (
            "omegaRes",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh,    
#ifdef CompressibleFlow
          dimensionedScalar("omegaRes", dimensionSet(1, -3, -2, 0, 0, 0, 0), 0.0),
#endif
#ifdef IncompressibleFlow
          dimensionedScalar("omegaRes", dimensionSet(0, 0, -2, 0, 0, 0, 0), 0.0),
#endif
            zeroGradientFvPatchField<scalar>::typeName),
    omegaResRef_(
          IOobject(
              "omegaResRef",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          omegaRes_),
      omegaResPartDeriv_(
          IOobject(
              "omegaResPartDeriv",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          omegaRes_),
      omegaRef_(
          IOobject(
              "omegaRef",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          omega_),
    Ar_(const_cast<volScalarField&>(
        mesh_.thisDb().lookupObject<volScalarField>("Ar"))),
    ArRes_(
        IOobject
        (
            "ArRes",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
#ifdef CompressibleFlow
          dimensionedScalar("ArRes", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0),
#endif
#ifdef IncompressibleFlow
          dimensionedScalar("ArRes", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0),
#endif
            zeroGradientFvPatchField<scalar>::typeName),
    ArResRef_(
          IOobject(
              "ArResRef",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          ArRes_),
      ArResPartDeriv_(
          IOobject(
              "ArResPartDeriv",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          ArRes_),
      ArRef_(
          IOobject(
              "ArRef",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          Ar_),
    v2_(const_cast<volScalarField&>(
        mesh_.thisDb().lookupObject<volScalarField>("v2"))),
    v2Res_(
        IOobject
        (
            "v2Res",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
#ifdef CompressibleFlow
          dimensionedScalar("v2Res", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0),
#endif
#ifdef IncompressibleFlow
          dimensionedScalar("v2Res", dimensionSet(0, 2, -3, 0, 0, 0, 0), 0.0),
#endif
            zeroGradientFvPatchField<scalar>::typeName),
    v2ResRef_(
          IOobject(
              "v2ResRef",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          v2Res_),
      v2ResPartDeriv_(
          IOobject(
              "v2ResPartDeriv",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          v2Res_),
      v2Ref_(
          IOobject(
              "v2Ref",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          v2_),
    k_(const_cast<volScalarField&>(
        mesh_.thisDb().lookupObject<volScalarField>("k"))),
    kRes_(
        IOobject
        (
            "kRes",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
#ifdef CompressibleFlow
          dimensionedScalar("kRes", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0),
#endif
#ifdef IncompressibleFlow
          dimensionedScalar("kRes", dimensionSet(0, 2, -3, 0, 0, 0, 0), 0.0),
#endif
            zeroGradientFvPatchField<scalar>::typeName),
    kResRef_(
          IOobject(
              "kResRef",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          kRes_),
      kResPartDeriv_(
          IOobject(
              "kResPartDeriv",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          kRes_),
      kRef_(
          IOobject(
              "kRef",
              mesh.time().timeName(),
              mesh,
              IOobject::NO_READ,
              IOobject::NO_WRITE),
          k_),
    
    // field inversion parameters
   betaFieldInversion_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("betaFieldInversion"))),
      surfaceFriction_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("surfaceFriction"))),
      surfaceFrictionData_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("surfaceFrictionData"))),
      pData_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("pData"))),
      UData_(const_cast<volVectorField&>(
          mesh.thisDb().lookupObject<volVectorField>("UData"))),
      USingleComponentData_(const_cast<volScalarField&>(
          mesh.thisDb().lookupObject<volScalarField>("USingleComponentData"))),     
     y_(mesh_.thisDb().lookupObject<volScalarField>("yWall")) 
{
    // initialize printInterval_ we need to check whether it is a steady state
    // or unsteady primal solver
    IOdictionary fvSchemes(
        IOobject(
            "fvSchemes",
            mesh.time().system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false));
    word ddtScheme = word(fvSchemes.subDict("ddtSchemes").lookup("default"));
    if (ddtScheme == "steadyState")
    {
        printInterval_ =
            daOption.getAllOptions().lookupOrDefault<label>("printInterval", 100);
    }
    else
    {
        printInterval_ =
            daOption.getAllOptions().lookupOrDefault<label>("printIntervalUnsteady", 500);
    }

    // calculate the size of omegaWallFunction faces
    label nWallFaces = 0;
    forAll(omega_.boundaryField(), patchI)
    {
        if (omega_.boundaryField()[patchI].type() == "omegaWallFunction"
            && omega_.boundaryField()[patchI].size() > 0)
        {
            forAll(omega_.boundaryField()[patchI], faceI)
            {
                nWallFaces++;
            }
        }
    }

    // initialize omegaNearWall
    omegaNearWall_.setSize(nWallFaces);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// kv2OmegaAr member functions
tmp<volScalarField> DAkv2OmegaArFieldInversion::fv(const volScalarField& Ret) const
{
    return tmp<volScalarField>(new volScalarField(
        "fv",
        1.0 - exp(-sqrt(Ret)/Anu_)
    ));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::fINT() const
{
    return tmp<volScalarField>(new volScalarField(
        "fINT",
        (
            min
            (
                v2_ / (CINT_* ( k_ + this->kMin_)),
                dimensionedScalar("1.0", dimless, 1.0)
            )
        )
    ));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::fSS(const volScalarField& W) const
{
    return tmp<volScalarField>(new volScalarField(
        "fSS",
        exp( -sqr(CSS_ * this->nu() * W / (v2_ + this->kMin_) ) )
    ));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::Cmu(const volScalarField& S) const
{
    return tmp<volScalarField>(new volScalarField(
        "Cmu",
        1.0/(A0_ + AS_*(S/max(omega_,this->omegaMin_)))
    ));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::betaTS(const volScalarField& ReW) const
{
    return tmp<volScalarField>(new volScalarField(
            "betaTS",
            scalar(1) - exp( -sqr( max(ReW - CTScrit_, scalar(0)) ) /ATS_ )
        ));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::lambdaT() const
{
  
    return tmp<volScalarField>(new volScalarField(
            "lambdaT",
            sqrt(this->v2_) / max(this->omega_, this->omegaMin_)
    ));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::lambdaEff(const volScalarField& lambdaT) const
{
  
    return tmp<volScalarField>(new volScalarField(
            "lambdaEff",
            // min( this->Clambda_ * y_, lambdaT)
            min( this->Clambda_ * y_, lambdaT)*(1.0+this->CA1_*pow(this->Ar_,this->CA2_))
            // min( this->Clambda_ * y_ *1.2, lambdaT) // pass
        ));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::fTaul
(
    const volScalarField& lambdaEff,
    const volScalarField& v2l,
    const volScalarField& W
) const
{
    const dimensionedScalar vMin("ROOTVSMALL",  dimLength*inv(dimTime), ROOTVSMALL);

    return tmp<volScalarField>(new volScalarField(
        "fTaul",
        scalar(1)
        - exp
        (
            -Ctau1_ * v2l
            /
            sqr( max( lambdaEff * W, vMin ) )
        )
    ));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::alphaT
(
    const volScalarField& lambdaEff,
    const volScalarField& fv,
    const volScalarField& v2s
) const
{
    return tmp<volScalarField>(new volScalarField(
        "alphaT",
        fv * betaStar_ * sqrt(v2s) * lambdaEff
    ));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::phiBP(const volScalarField& W) const
{
    const dimensionedScalar wMin("ROOTVSMALL", inv(dimTime), ROOTVSMALL);

    return tmp<volScalarField>(new volScalarField(
        "phiBP",
        min
        (
            max
            (
                v2_/ ( this->nu() * max(W, wMin) ) - CBPcrit_, 
                scalar(0)
            ),
            scalar(50.0)
        )
    ));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::phiNAT
(
    const volScalarField& ReW,
    const volScalarField& fNATcrit
) const
{
    const dimensionedScalar small("ROTVSMALL", dimless, ROOTVSMALL);

    return tmp<volScalarField>(new volScalarField(
        "phiNAT",
        max
        (
            ReW - CNATcrit_ / max(fNATcrit, small),
            scalar(0)
        )
    ));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::D(const volScalarField& k) const
{
    return 2.0*this->nu()*magSqr(fvc::grad(sqrt(k)));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::F1() const
{
    const volScalarField CDkOmega(
        "CDkOmega",
        max(
            (scalar(2) * this->sigmaW2_)*
            (fvc::grad(this->k_) & fvc::grad(this->omega_))/this->omega_,
            dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
        )
    );

    const volScalarField arg1(
        "arg1",
        min(
            max(
                sqrt(this->v2_)/(this->omega_*this->y_),
                500.0 * this->nu() * this->betaStar_ / (sqr(this->y_) * this->omega_)
            ),
            (scalar(4.0) * this->sigmaW2_) * this->k_ / (CDkOmega * sqr(this->y_))
            ));

    return tmp<volScalarField>( new volScalarField(
        "F1",
        tanh(pow4(arg1))
    ));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::W(const volTensorField& gradU) const
{

    return 
        ::sqrt(2.0)*mag(skew(gradU));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::S2(const volTensorField& gradU) const
{

    return 
        2.0*magSqr(symm(gradU));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::fW(const volScalarField& lambdaT_,const volScalarField& lambdaEff_ ) const
{
    return 
        min(
            pow(
                lambdaEff_
                /max(lambdaT_,dimensionedScalar("SMALL", dimLength, ROOTVSMALL)),
                2.0/3.0
            ),
            dimensionedScalar("1.0", dimless, 1.0)
        );
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::v2s(const volScalarField& W ,const volScalarField& fW) const
{
    volScalarField r(fSS(W) * fW * v2_);

    return r;
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::ReT(const volScalarField& fW) const
{
    return 
        sqr(fW) * v2_ / (this->nu() * max(omega_, omegaMin_));
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::nuTs(const volScalarField& fv_,const volScalarField& S2,const volScalarField& v2s,const volScalarField& lambdaEff_) const
{
    return
        fv_ * fINT() * Cmu(sqrt(S2)) * sqrt(v2s) * lambdaEff_;
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::v2l(const volScalarField& v2s) const
{
    return 
        v2_ - v2s;
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::ReW(const volScalarField& W) const
{
    return
        sqr(y_) * W / this->nu();
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::nuTl(
    const volScalarField& lambdaEff_,
    const volScalarField& v2l,
    const volScalarField& W,
    const volScalarField& ReW,
    const volScalarField& v2s,
    const volScalarField& S2) const
{
    return 
        min
      (
          C11_* fTaul(lambdaEff_, v2l, W) * W * sqr(lambdaEff_) 
          * sqrt(v2l) * lambdaEff_ /this->nu()
          + C12_ * betaTS(ReW) * pow4(lambdaEff_/Clambda_) * sqr(W) / this->nu()
          ,
          0.5*(k_ - v2s)/max(sqrt(S2), omegaMin_)
      );
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::F1star(const volScalarField& W) const
{
    return 
        scalar(1) - (scalar(1)-F1())*fSS(W);
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::RBP(const volScalarField& W,const volScalarField& fW,const dimensionedScalar& small) const
{
    return 
        CR_ * (1.0 - exp(-phiBP(W)/ABP_))*omega_
      /max(fW, small);
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::fNATcrit() const
{
    return 
        1.0 - exp(-CNC_*sqrt(k_)*y_/this->nu());
}

tmp<volScalarField> DAkv2OmegaArFieldInversion::RNAT(const volScalarField& ReW,const volScalarField& fNATcrit,const volScalarField& W) const
{
    return 
        CRNAT_ * (1.0 - exp(-phiNAT(ReW, fNATcrit)/ANAT_))*W;
}

// Augmented functions
void DAkv2OmegaArFieldInversion::correctModelStates(wordList& modelStates) const
{
    /*
    Description
        Update the name in modelStates based on the selected physical model at runTimle

    Example:
        for kv2OmegaAr, modelStates = {"k", "omega", "Ar", "v2"}    
    
    */

    // For kv2OmegaAr, we need to replace nut with k, omega, Ar, v2

    forAll(modelStates,idxI)
    {
        word stateName = modelStates[idxI];
        if (stateName == "nut")
        {
            modelStates[idxI] = "omega";
            modelStates.append("k");
            modelStates.append("Ar");
            modelStates.append("v2");
        }
    }
}

void DAkv2OmegaArFieldInversion::correctNut()
{
    const dimensionedScalar small("ROTVSMALL", dimless, ROOTVSMALL);

    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU = tgradU();
    const volScalarField W(this->W(gradU));
    const volScalarField S2(this->S2(gradU));

    const volScalarField lambdaT_ = this->lambdaT();

    const volScalarField lambdaEff_ = this -> lambdaEff(lambdaT_);

    const volScalarField fW(this->fW(lambdaT_, lambdaEff_));

    const volScalarField v2s(this -> v2s(W, fW));

    tmp<volScalarField> ReT(this -> ReT(fW));
    
    const volScalarField fv_(this -> fv(ReT));

    const volScalarField nuTs(this -> nuTs(fv_,S2, v2s, lambdaEff_));
    

    const volScalarField v2l("v2l", v2_ - v2s);
    const volScalarField ReW("ReW", sqr(y_) * W / this->nu() );

    const volScalarField nuTl(this -> nuTl(lambdaEff_, v2l, W, ReW, v2s, S2));
    // Re-calculate turbulent viscosity
    nut_ = nuTs + nuTl;

    dimensionedScalar small2 = dimensionedScalar("small2", (dimLength*dimLength/dimTime), 1e-30);
    bound(nut_,small2);

    nut_.correctBoundaryConditions();

    this->correctAlphat();

    return;
}

void DAkv2OmegaArFieldInversion::correctBoundaryConditions()
{
    /*
    Description:
        Update turbulence variable boundary values
    */

    // correct the BCs for the perturbed fields
    // kqWallFunction is a zero-gradient BC
    k_.correctBoundaryConditions();
}

void DAkv2OmegaArFieldInversion::correctOmegaBoundaryConditions()
{
    /*
    Description:
        this is a special treatment for omega BC because we cant directly call omega.
        correctBoundaryConditions() because it will modify the internal omega and G that 
        are right next to walls. This will mess up adjoint Jacobians
        To solve this issue,
        1. we store the near wall omega before calling omega.correctBoundaryConditions()
        2. call omega.correctBoundaryConditions()
        3. Assign the stored near wall omega back
        4. Apply a zeroGradient BC for omega at the wall patches
        *********** NOTE *************
        this treatment will obviously downgrade the accuracy of adjoint derivative since it is 
        not 100% consistent with what is used for the flow solver; however, our observation 
        shows that the impact is not very large.
    */

    // save the perturbed omega at the wall
    this->saveOmegaNearWall();
    // correct omega boundary conditions, this includes updating wall face and near wall omega values,
    // updating the inter-proc BCs
    omega_.correctBoundaryConditions();
    // reset the corrected omega near wall cell to its perturbed value
    this->setOmegaNearWall();
}

void DAkv2OmegaArFieldInversion::saveOmegaNearWall()
{
    /*
    Description:
        Save the current near wall omega values to omegaNearWall_
    */
    label counterI = 0;
    forAll(omega_.boundaryField(), patchI)
    {
        if (omega_.boundaryField()[patchI].type() == "omegaWallFunction"
            and omega_.boundaryField()[patchI].size() > 0)
        {
            const UList<label>& faceCells = mesh_.boundaryMesh()[patchI].faceCells();
            forAll(faceCells, faceI)
            {
                //Info<<"Near Wall cellI: "<<faceCells[faceI]<<endl;
                omegaNearWall_[counterI] = omega_[faceCells[faceI]];
                counterI++;
            }
        }
    }
    return;
}

void DAkv2OmegaArFieldInversion::setOmegaNearWall()
{
    /*
    Description:
        Set the current near wall omega values from omegaNearWall_
        Here we also apply a zeroGradient BC to the wall faces
    */
    label counterI = 0;
    forAll(omega_.boundaryField(), patchI)
    {
        if (omega_.boundaryField()[patchI].type() == "omegaWallFunction"
            && omega_.boundaryField()[patchI].size() > 0)
        {
            const UList<label>& faceCells = mesh_.boundaryMesh()[patchI].faceCells();
            forAll(faceCells, faceI)
            {
                omega_[faceCells[faceI]] = omegaNearWall_[counterI];
                // zeroGradient BC
                omega_.boundaryFieldRef()[patchI][faceI] = omega_[faceCells[faceI]];
                counterI++;
            }
        }
    }
    return;
}

void DAkv2OmegaArFieldInversion::updateIntermediateVariables()
{
    /*
    Description:
        Update nut based on nuTilda. Note: we need to update nut and its BC since we 
        may have perturbed other turbulence vars that affect the nut values
    */

    this->correctNut();
}

void DAkv2OmegaArFieldInversion::correctStateResidualModelCon(List<List<word>>& stateCon) const
{
    /*
    Description:
        Update the original variable connectivity for the adjoint state 
        residuals in stateCon. Basically, we modify/add state variables based on the
        original model variables defined in stateCon.
    
    stateResCon: the connectivity levels for a state residual, defined in Foam::DAJacCon

    Example:
        for kv2OmegaAr, stateCon = {{"U","p”，"k", "omega", "Ar", "v2"}, {"U", "p"}}

    */
    
    label stateConSize = stateCon.size();
    forAll(stateCon, idxI)
    {
        label addUCon = 0;
        forAll(stateCon[idxI], idxJ)
        {
            word conStateName = stateCon[idxI][idxJ];
            if (conStateName == "nut")
            {
                stateCon[idxI][idxJ] = "omega";
                stateCon[idxI].append("k");
                stateCon[idxI].append("Ar");
                stateCon[idxI].append("v2");
                addUCon = 1;
            }
        }
        // add U for the current level and level+1 if it is not there yet
        label isU;
        if (addUCon == 1)
        {
            // first add U for the current level
            isU = 0;
            forAll(stateCon[idxI], idxJ)
            {
                word conStateName = stateCon[idxI][idxJ];
                if (conStateName == "U")
                {
                    isU = 1;
                }
            }
            if (!isU)
            {
                stateCon[idxI].append("U");
            }

            // now add U for level+1 if idxI is not the largest level
            // if idxI is already the largest level, we have a problem
            if (idxI != stateConSize - 1)
            {
                isU = 0;
                forAll(stateCon[idxI + 1], idxJ)
                {
                    word conStateName = stateCon[idxI + 1][idxJ];
                    if (conStateName == "U")
                    {
                        isU = 1;
                    }
                }
                if (!isU)
                {
                    stateCon[idxI + 1].append("U");
                }
            }
            else
            {
                FatalErrorIn(
                    "In DAStateInfo, nut shows in the largest connectivity level! "
                    "This is not supported!")
                    << exit(FatalError);
            }
        }
    }
}

void DAkv2OmegaArFieldInversion::addModelResidualCon(HashTable<List<List<word>>>& allCon) const
{
    /*
    Description:
        Add the model residual connectivity to stateCon. This is used to 
        update the adjoint residual connectivity in DAStateInfo::regStates_.
        Basically, we add the model variables to the stateCon if they are not there yet.
    
    allCon: the connectivity levels for all residuals, defined in Foam::DAJacCon

    Example:
        for kv2OmegaAr, allCon = {{"U","p”，"k", "omega", "Ar", "v2"}, {"U", "p"}}

    */

   word pName;

   if (mesh_.thisDb().foundObject<volScalarField>("p"))
   {
         pName = "p";
    }
    else if (mesh_.thisDb().foundObject<volScalarField>("p_rgh"))
    {
         pName = "p_rgh";
    }
    else
    {
        FatalErrorIn(
            "Neither p nor p_rgh was found in mesh.thisDb()!"
            "addModelResidualCon failed to setup turbulence residuals!")
            << exit(FatalError);
    }

    // NOTE: for compressible flow, it depends on rho so we need to add T and p
#ifdef IncompressibleFlow
    allCon.set(
        "omegaRes",
        {
            {"U", "omega", "k", "v2" ,"Ar","phi"}, // lv0
            {"U", "omega", "k", "v2" ,"Ar"}, // lv1
            {"U", "omega", "k", "v2" ,"Ar"} // lv2
        });
    allCon.set(
        "kRes",
        {
            {"U", "omega", "k", "v2" ,"Ar","phi"}, // lv0
            {"U", "omega", "k", "v2" ,"Ar"}, // lv1
            {"U", "omega", "k", "v2" ,"Ar"} // lv2
        });
    allCon.set(
        "ArRes",
        {
            {"U", "omega", "k", "v2" ,"Ar","phi"}, // lv0
            {"U", "omega", "k", "v2" ,"Ar"}, // lv1
            {"U", "omega", "k", "v2" ,"Ar"} // lv2
        });
    allCon.set(
        "v2Res",
        {
            {"U", "omega", "k", "v2" ,"Ar","phi"}, // lv0
            {"U", "omega", "k", "v2" ,"Ar"}, // lv1
            {"U", "omega", "k", "v2" ,"Ar"} // lv2
        });
#endif

#ifdef CompressibleFlow
    allCon.set(
        "omegaRes",
        {
            {"U", "omega", "k", "v2" ,"Ar","phi","T",pName}, // lv0
            {"U", "omega", "k", "v2" ,"Ar","T",pName}, // lv1
            {"U", "omega", "k", "v2" ,"Ar","T",pName} // lv2
        });
    allCon.set(
        "kRes",
        {
            {"U", "omega", "k", "v2" ,"Ar","phi","T",pName}, // lv0
            {"U", "omega", "k", "v2" ,"Ar","T",pName}, // lv1
            {"U", "omega", "k", "v2" ,"Ar","T",pName} // lv2
        });
    allCon.set(
        "ArRes",
        {
            {"U", "omega", "k", "v2" ,"Ar","phi","T",pName}, // lv0
            {"U", "omega", "k", "v2" ,"Ar","T",pName}, // lv1
            {"U", "omega", "k", "v2" ,"Ar","T",pName} // lv2
        });
    allCon.set(
        "v2Res",
        {
            {"U", "omega", "k", "v2" ,"Ar","phi","T",pName}, // lv0
            {"U", "omega", "k", "v2" ,"Ar","T",pName}, // lv1
            {"U", "omega", "k", "v2" ,"Ar","T",pName} // lv2
        });
#endif
}

void DAkv2OmegaArFieldInversion::correct()
{
    /*
    Descroption:
        Solve the residual equations and update the state. This function will be called 
        by the DASolver. It is needed because we want to control the output frequency
        of the residual convergence every 100 steps. If using the correct from turbulence
        it will output residual every step which will be too much of information.
    */
        // We set the flag solveTurbState_ to 1 such that in the calcResiduals function
    // we will solve and update nuTilda
    solveTurbState_ = 1;
    dictionary dummyOptions;
    this->calcResiduals(dummyOptions);
    // after it, we reset solveTurbState_ = 0 such that calcResiduals will not
    // update nuTilda when calling from the adjoint class, i.e., solveAdjoint from DASolver.
    solveTurbState_ = 0;
}

void DAkv2OmegaArFieldInversion::calcResiduals(const dictionary& options)
{
    /*
    Descroption:
        If solveTurbState_ == 1, this function solve and update k and omega, and 
        is the same as calling turbulence.correct(). If solveTurbState_ == 0,
        this function compute residuals for turbulence variables, e.g., nuTildaRes_

    Input:
        options.isPC: 1 means computing residuals for preconditioner matrix.
        This essentially use the first order scheme for div(phi,nuTilda)

        p_, U_, phi_, etc: State variables in OpenFOAM
    
    Output:
        kRes_/omegaRes_: If solveTurbState_ == 0, update the residual field variable

        k_/omega_: If solveTurbState_ == 1, update them
    */

    // Copy and modify based on the "correct" function

    label printToScreen = this->isPrintTime(mesh_.time(), printInterval_);

    word divKScheme = "div(phi,k)";
    word divOmegaScheme = "div(phi,omega)";
    word divArScheme = "div(phi,Ar)";
    word divv2Scheme = "div(phi,v2)";

    label isPC = 0;

    if (!solveTurbState_)
    {
        isPC = options.getLabel("isPC");

        if (isPC)
        {
            divKScheme = "div(pc)";
            divOmegaScheme = "div(pc)";
            divArScheme = "div(pc)";
            divv2Scheme = "div(pc)";
        }
    }

    const dimensionedScalar small("ROTVSMALL", dimless, ROOTVSMALL);

    tmp<volTensorField> tgradU(fvc::grad(U_));
    const volTensorField& gradU = tgradU();
    const volScalarField W(this->W(gradU));
    const volScalarField S2(this->S2(gradU));


    const volScalarField lambdaT_ = this->lambdaT();

    // ******************************************** //

    tmp <fvScalarMatrix> ArEqn
    (
        fvm::ddt(phase_, rho_, Ar_)
        + fvm::div(phaseRhoPhi_, Ar_, divArScheme)
        ==
        fvm::laplacian(phase_*rho_*DArEff(), Ar_)     
    );

    ArEqn.ref().relax();

    if (solveTurbState_)
    {
            //  get the solver performance info such as initial
            //  and final residuals
            SolverPerformance<scalar> solverAr = solve(ArEqn);
            if(printToScreen)
            {
                Info<< "Ar Initial residual:" << solverAr.initialResidual() << endl
                    << "     Final residual:" << solverAr.finalResidual() << endl;
            }

            DAUtility::boundVar(allOptions_, Ar_, printToScreen);
    }
    else
    {
        // calculate residuals
        ArRes_ = ArEqn() & Ar_;
        // need to normalize residuals
        normalizeResiduals(ArRes);
    }

    // length scale need to be corrected by Ar
    const volScalarField lambdaEff_ = this -> lambdaEff(lambdaT_);

    // ******************************************** //

    const volScalarField fW(this -> fW(lambdaT_, lambdaEff_));

    const volScalarField v2s(this -> v2s(W, fW));

    tmp<volScalarField> ReT(this -> ReT(fW));
    
    const volScalarField fv_(this -> fv(ReT));

    const volScalarField nuTs(this -> nuTs(fv_,S2, v2s, lambdaEff_));
    

    const volScalarField v2l("v2l", v2_ - v2s);
    const volScalarField ReW("ReW", sqr(y_) * W / this->nu() );

    const volScalarField nuTl(this -> nuTl(lambdaEff_, v2l, W, ReW, v2s, S2));

    const volScalarField F1star(this -> F1star(W));

    const volScalarField alphaT_ = alphaT(lambdaEff_, fv_, v2s);

    const volScalarField RBP(this -> RBP(W, fW, small));

    const volScalarField fNATcrit(this -> fNATcrit());

    const volScalarField RNAT(this -> RNAT(ReW, fNATcrit, W));

    if (solveTurbState_)
    {
        // Update omega and G at the wall
        omega_.boundaryFieldRef().updateCoeffs();
    }
    else
    {
        // NOTE instead of calling omega_.boundaryFieldRef().updateCoeffs();
        // here we call our self-defined boundary conditions
        this->correctOmegaBoundaryConditions();
    }

    tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(phase_, rho_, omega_)
            + fvm::div(phaseRhoPhi_, omega_, divOmegaScheme)
            - fvm::laplacian(phase_*rho_*DomegaEff(alphaT_), omega_)
            ==
            phase_ * rho_ * (
                Cw1_ * omega_ / max(v2_, kMin_) * nuTs * S2 * betaFieldInversion_
                - fvm::SuSp(
                    (1.0 - CwR_/max(fW,small)) * (k_ - v2_) * (RBP + RNAT)
                    /max(v2_, kMin_)
                    , omega_)
                - fvm::Sp(Cw2_*sqr(fW)*omega_, omega_)
                + betaStar_ * 2 * (1.0 - F1star) * sigmaW2_ / max(omega_, omegaMin_) *
                (fvc::grad(k_) & fvc::grad(omega_) )
            )
        );
    
    omegaEqn.ref().relax();
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());

    if (solveTurbState_)
    {

        // get the solver performance info such as initial
        // and final residuals
        SolverPerformance<scalar> solverOmega = solve(omegaEqn);
        if (printToScreen)
        {
            Info << "omega Initial residual: " << solverOmega.initialResidual() << endl
                 << "        Final residual: " << solverOmega.finalResidual() << endl;
        }

        DAUtility::boundVar(allOptions_, omega_, printToScreen);
    }
    else
    {
        // reset the corrected omega near wall cell to its perturbed value
        this->setOmegaNearWall();

        // calculate residuals
        omegaRes_ = omegaEqn() & omega_;
        // need to normalize residuals
        normalizeResiduals(omegaRes);
    }

    // ******************************************** //
    tmp<fvScalarMatrix> v2Eqn
        (
            fvm::ddt(phase_, rho_, v2_)
            + fvm::div(phaseRhoPhi_, v2_)
            - fvm::laplacian(phase_*rho_*DkEff(alphaT_), v2_)
            ==
            phase_ * rho_ * (
                nuTs * S2 
                + (RBP + RNAT) * k_ 
                - fvm::Sp(RBP + RNAT, v2_)
                - fvm::Sp(omega_ + D(v2_)/max(v2_,kMin_), v2_)
            )
    );

    v2Eqn.ref().relax();
    v2Eqn.ref().boundaryManipulate(v2_.boundaryFieldRef());

    if(solveTurbState_)
    {
        // get the solver performance info such as initial
        // and final residuals
        SolverPerformance<scalar> solverv2 = solve(v2Eqn);
        if(printToScreen)
        {
            Info<< "v2 Initial residual:" << solverv2.initialResidual() << endl
                << "     Final residual:" << solverv2.finalResidual() << endl;
        }

        DAUtility::boundVar(allOptions_, v2_, printToScreen);
    }
    else
    {
        // calculate residuals
        v2Res_ = v2Eqn() & v2_;
        // need to normalize residuals
        normalizeResiduals(v2Res);
    }

    // ******************************************** //
    tmp<fvScalarMatrix> kEqn
        (
            fvm::ddt(phase_, rho_, k_)
            + fvm::div(phaseRhoPhi_, k_)
            - fvm::laplacian(phase_*rho_*DkEff(alphaT_), k_)
            ==
            phase_*rho_*(
                nut_ * S2
                - fvm::Sp((omega_*min(k_,v2_) + D(k_))/max(k_,kMin_), k_)
            )
        );
    kEqn.ref().relax();
    kEqn.ref().boundaryManipulate(k_.boundaryFieldRef());

    if(solveTurbState_)
    {
        // get the solver performance info such as initial
        // and final residuals
        SolverPerformance<scalar> solverk = solve(kEqn);
        if(printToScreen)
        {
            Info<< "k Initial residual:" << solverk.initialResidual() << endl
                << "     Final residual:" << solverk.finalResidual() << endl;
        }

        DAUtility::boundVar(allOptions_, k_, printToScreen);
        
        this->correctNut();
    }

    else
    {
        // calculate residuals
        kRes_ = kEqn() & k_;
        // need to normalize residuals
        normalizeResiduals(kRes);
    }

    return;
}

// ******************************************** //

} // namespace Foam

// ******************************************** //

