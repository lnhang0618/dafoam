/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

    Description:
    Child class for the kv2OmegaAr model

    This file is modified from OpenFOAM's source code
    src/TurbulenceModels/turbulenceModels/RAS/kOmega/kOmega.H

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

#ifndef DAkv2OmegaArFieldInversion_H
#define DAkv2OmegaArFieldInversion_H

#include "DATurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class DAkv2OmegaArFieldInversion Declaration
\*---------------------------------------------------------------------------*/

class DAkv2OmegaArFieldInversion
    : public DATurbulenceModel
{

protected:

    // Protected data
    
    /// \name kv2OmegaAr parameters
    //@{
    dimensionedScalar A0_;
    dimensionedScalar AS_;
    dimensionedScalar Anu_;
    dimensionedScalar ABP_;
    dimensionedScalar ANAT_;
    dimensionedScalar ATS_;
    dimensionedScalar CBPcrit_;
    dimensionedScalar CNC_;
    dimensionedScalar CNATcrit_;
    dimensionedScalar CINT_;
    dimensionedScalar CTScrit_;
    dimensionedScalar CRNAT_;
    dimensionedScalar C11_;
    dimensionedScalar C12_;
    dimensionedScalar CR_;
    dimensionedScalar CalphaTheta_;
    dimensionedScalar CSS_;
    dimensionedScalar Ctau1_;
    dimensionedScalar Cw1_;
    dimensionedScalar Cw2_;
    dimensionedScalar CwR_;
    dimensionedScalar Clambda_;
    dimensionedScalar betaStar_;
    dimensionedScalar PrTheta_;
    dimensionedScalar sigmaK_;
    dimensionedScalar sigmaW_;
    dimensionedScalar sigmaW2_;
    dimensionedScalar sigmaAr_;
    dimensionedScalar CA1_;
    dimensionedScalar CA2_;
    //@}

    /// \name Augmented variables for adjoint residuals
    ///@{
    volScalarField& omega_;
    volScalarField omegaRes_;
    volScalarField omegaResRef_;
    volScalarField omegaResPartDeriv_;
    volScalarField omegaRef_;
    volScalarField& Ar_;
    volScalarField ArRes_;
    volScalarField ArResRef_;
    volScalarField ArResPartDeriv_;
    volScalarField ArRef_;
    volScalarField& v2_;
    volScalarField v2Res_;
    volScalarField v2ResRef_;
    volScalarField v2ResPartDeriv_;
    volScalarField v2Ref_;
    volScalarField& k_;
    volScalarField kRes_;
    volScalarField kResRef_;
    volScalarField kResPartDeriv_;
    volScalarField kRef_;
    //@}

        /// cell-center omega values near the wall, this is to fix the issue that the
    /// omegaWallFunction will try to modify omega values for the cells near walls
    /// this will cause issue for FD-based partial derivatives, so here we basically
    /// implement a zeroGradient BC for near wall omega.
    scalarList omegaNearWall_;

    /// whether to solve for turb states
    label solveTurbState_ = 0;

    /// time step interval to print residual
    label printInterval_;
        
    /// fields
    volScalarField& betaFieldInversion_;

    /// a surface friction 'field' when using skin friction data for field inversion
    volScalarField surfaceFriction_;

    /// the reference field for surfaceFriction
    volScalarField surfaceFrictionData_; 

    /// the reference field for surfacePressure
    volScalarField pData_; 
    
    /// reference field (e.g. velocity for DNS)
    volVectorField UData_; 

    /// the reference pressure field data
    volScalarField USingleComponentData_;

    //- Wall distance
    const volScalarField& y_;

    /// \name kv2OmegaAr functions
    //@{
    tmp<volScalarField> fv(const volScalarField& Ret) const;

    tmp<volScalarField> fINT() const;

    tmp<volScalarField> fSS(const volScalarField& W) const;

    tmp<volScalarField> Cmu(const volScalarField& S) const;

    tmp<volScalarField> betaTS(const volScalarField& ReW) const;

    tmp<volScalarField> lambdaT() const;

    tmp<volScalarField> lambdaEff(const volScalarField& lambdaT) const;

    tmp<volScalarField> fTaul
    (
        const volScalarField& lambdaEff,
        const volScalarField& v2l,
        const volScalarField& W
    ) const;

    tmp<volScalarField> alphaT
    (
        const volScalarField& lambdaEff,
        const volScalarField& fv,
        const volScalarField& v2s
    ) const;

    tmp<volScalarField> phiBP(const volScalarField& W) const;

    tmp<volScalarField> phiNAT
    (
        const volScalarField& ReOmega,
        const volScalarField& fNatCrit
    ) const;


    tmp<volScalarField> D(const volScalarField& k) const;

    tmp<volScalarField> F1() const;

    tmp<volScalarField> W(const volTensorField& gradU) const;

    tmp<volScalarField> S2(const volTensorField& gradU) const;

    tmp<volScalarField> fW(const volScalarField& lambdaT_,const volScalarField& lambdaEff_) const;

    tmp<volScalarField> v2s(const volScalarField& W ,const volScalarField& fW) const;

    tmp<volScalarField> ReT(const volScalarField& fW) const;

    tmp<volScalarField> nuTs(const volScalarField& fv_,const volScalarField& S2,const volScalarField& v2s ,const volScalarField& lambdaEff_) const;

    tmp<volScalarField> v2l(const volScalarField& v2s) const;

    tmp<volScalarField> ReW(const volScalarField& W) const;

    tmp<volScalarField> nuTl(
        const volScalarField& lambdaEff_,
        const volScalarField& v2l,
        const volScalarField& W,
        const volScalarField& ReW,
        const volScalarField& v2s,
        const volScalarField& S2) const;

    tmp<volScalarField> F1star(const volScalarField& W) const;

    tmp<volScalarField> RBP(const volScalarField& W, const volScalarField& fW,const dimensionedScalar& small) const;

    tmp<volScalarField> fNATcrit() const;

    tmp<volScalarField> RNAT(const volScalarField& ReW, const volScalarField& fNATcrit,const volScalarField& W) const;


    //- Return the effective diffusivity for k
    tmp<volScalarField> DkEff(const volScalarField& alphaT) const
    {
        return tmp<volScalarField>
        (
            new volScalarField("DkEff", alphaT/sigmaK_ + this->nu())
        );
    }

    //- Return the effective diffusivity for omega
    tmp<volScalarField> DomegaEff(const volScalarField& alphaT) const
    {
        return tmp<volScalarField>
        (
            new volScalarField("DomegaEff", alphaT/sigmaW_ + this->nu())
        );
    }

    //- Return the effective diffusivity for roughness amplification
    tmp<volScalarField> DArEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DArEff",
                    (this->nut_ + this->nu()) * betaFieldInversion_
                )
            );
        }
    tmp<volScalarField> epsilon() const
        {
            return tmp<volScalarField>
                (
                    new volScalarField(
                        "epsilon", 
                        omega_ * min(k_, v2_) + D(k_)
                    ));
        }
    

public:
    TypeName("kv2OmegaArFieldInversion");
    // Constructor

    //- Construct from components
    DAkv2OmegaArFieldInversion(
        const word &modelType,
        const fvMesh& mesh,
        const DAOption& daOption);

    //- Destructor
    virtual ~DAkv2OmegaArFieldInversion() 
    {
    }

    // Member functions

    /// update the turbulence state for DAStateInfo::regStates_
    virtual void correctModelStates(wordList& modelStates) const;

    /// update nut based on other turbulence variables and update the BCs
    virtual void correctNut();

    /// update turbulence variable boundary values
    virtual void correctBoundaryConditions();

    /// update any intermediate variables that are dependent on state variables and are used in calcResiduals
    virtual void updateIntermediateVariables();

    /// update the original variable connectivity for the adjoint state residuals in stateCon
    virtual void correctStateResidualModelCon(List<List<word>>& stateCon) const;

    /// add the model residual connectivity to stateCon
    virtual void addModelResidualCon(HashTable<List<List<word>>>& allCon) const;

    /// compute the turbulence residuals
    virtual void calcResiduals(const dictionary& options);

    /// solve the residual equations and update the state
    virtual void correct();

    /// save near wall omega values to omegaNearWall_
    void saveOmegaNearWall();

    /// set omegaNearWall_ to near wall omega values
    void setOmegaNearWall();

    /// specially treatment to correct epsilon BC
    void correctOmegaBoundaryConditions();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
