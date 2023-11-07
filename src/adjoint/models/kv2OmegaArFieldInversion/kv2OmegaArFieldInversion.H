/*---------------------------------------------------------------------------*\

    DAFoam  : Discrete Adjoint with OpenFOAM
    Version : v3

    Description:
        kv2Omega_Ar turbulence model with FieldInversion fields.
    
    Note: this is just a dummy turbulence model implementation, just to
    create the variable k, v2, Ar, omega. We shouldn't call any function in this class
    Instead, we should call the correct function in DAModel/DATurbulenceModel

\*---------------------------------------------------------------------------*/
#ifndef kv2OmegaArFieldInversion_H
#define kv2OmegaArFieldInversion_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class kv2OmegaArFieldInversion Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kv2OmegaArFieldInversion
    : public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

    // Disallow default bitwise copy construct and assignment
    kv2OmegaArFieldInversion(const kv2OmegaArFieldInversion&);
    void operator=(const kv2OmegaArFieldInversion&);

protected:

    // Protected data

    // Fields

    // turbulence variables for kv2Omega_Ar
    volScalarField omega_;
    volScalarField k_;
    volScalarField v2_;
    volScalarField Ar_;

    // A beta field multiplying to the production term
    volScalarField betaFieldInversion_;

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
    //  Note: different to wall distance in parent RASModel
    //  which is for near-wall cells only
    const volScalarField& y_;

    virtual void correctNut();

public:
    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    //- Runtime type information
    TypeName("kv2OmegaArFieldInversion");

    // Constructors

    //- Construct from components
    kv2OmegaArFieldInversion(
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const transportModel& transport,
        const word& propertiesName = turbulenceModel::propertiesName,
        const word& type = typeName);

    //- Destructor
    virtual ~kv2OmegaArFieldInversion() {}

    // Member functions

    //- Re-read model coefficients if they have changed
    virtual bool read();

    // -Return the turbulence model Fields
    virtual tmp<volScalarField> k() const;
    virtual tmp<volScalarField> epsilon() const;

    //- Solver the turbulence equations and correct the turbulence viscosity
    virtual void correct();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // namespace RASModels
} // namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
# include "kv2OmegaArFieldInversion.C"
#endif 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif