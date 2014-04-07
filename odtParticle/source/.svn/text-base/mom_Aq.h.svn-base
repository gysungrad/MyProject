/**
 * @file mom_Aq.h
 * Header file for class mom 
 */
#include "diffuser.h"
#include "ETA_Aq.h"
#include "odtline.h"
#include "pdgen.h"
#ifndef MOM_Aq_H    
#define MOM_Aq_H    

///////////////////////////////////////////////////////////////////////////////

/** Class for implementing particle moment equations
 *  
 *  @author David O. Lignell
 */

class mom_Aq : public MOM { 

  ////////////////////// DATA MEMBERS /////////////////////

    public:
    std::vector<double>          Source;        ///< source term with # of moments
  //  int                      *np;        ///< Number of gridpoints
    int                       nmom;         ///< number of moments
    vector<double>           gamma;        ///< Surface energy (J/cm^2)
    vector<double>         gamma_o;        ///< Surface energy (J/cm^2) at 25 C
    vector<vector<double> >gammaAC;        ///< Surface energy between polymorphs (J/cm^2)
    vector<vector<double> >     mp;        ///< m for heterogenous nucleation
    vector<double>            mu_0;        ///< initial value of the moments
    bool                      Lhet;        ///< Flag for heterogenous nucleation
    bool                   Lvanoss;        ///< Flag for temperature dependant interfacial energies
    double                     M3i;        ///< proportional to the inital vol of particles
    double                   r1   ;        ///< Nuclation size (cm)
    vector<double>           ro   ;        ///< Nuclation parameter (cm)
    double                   kb   ;        ///< Boltzman constant (J/K*(# of molecules))
    vector<double>           rhos ;        ///< Molar density of the solid (gmol/cm^3)
    double                   Diff ;        ///< Diffusivity of Liquid species to be precipiated
    double                   Diff0;        ///< Diffusion of CaCO3 at 25 C
    vector<double>           A    ;        ///< Constant for nucleation (unitless) (5-15-2012 not currently used)
    double                   B    ;        ///< Nucleation birth rate units of # of Nucleates/cm^3/s (sean has discussed the discrepency)  
    vector<double>           Ceq  ;        ///< Equilibrium Concentration of solid phase  
    double                   K    ;        ///< Toy Problem Rate constant (5-15-2012 not currently used)
    double                   SS   ;        ///< Super Satueration ratio
    double                   Jm   ;        ///< Nucleation Constant 1/cm^2*s (seans discrepency) (5-15-2012 not currently used)
    double                    Temp;        ///< Temperature
    double                Aqchem_S;        ///< Source term for the Aqueous species
    double                       R;        ///< Gas Constant
    double                     MWCC;       ///< Molecular Weight of Calcium Carbonate
    int                      npoly;        ///< Number of Polymorphs
    pda                          x;        ///< PD Algorithm Object 
    int                        np2;        ///< Number of weights or Abscissa 
    bool                      Ldqmom ;     ///< DQMOM on/off switch
    bool                      Lostwald;    ///< Ostwald ripening switch
    double                    R_gas;       ///< gas Constant
    double              CC_activity;       ///< Activity of CaCO3_o

    vector<double>              Ksp;        ///< Solubility Constant, K_sp  for polymorphs
    vector<double>              ap;        ///< slope of surface energy components for van oss 
    vector<double>              bp;        ///< intercept of surface energy components for van oss 

    vector<double>          gp_lw;        ///< Non-polar component Van Oss Equation 
    vector<double>           gp_p;        ///< Positive component in Van Oss Equation
    vector<double>           gp_m;        ///< Negative component in Van Oss Equation

    double                    unit;
    double                    unit2;
    double                    unit3;
    ETA_Aq                    *ETApx;
    odtline                      *odtlp;   ///<
    double  get_moment_source(vector<double>   mu, double Caq, double H, double CC_activity);        ///<source terms
    bool  check(vector<double>   mu);        ///<source terms
    void  EvalPropAtTemp(double H);
    void  AdaptiveODE(double dtStep, int &RD,std::vector<double> &dxML, std::vector<std::vector<double> > &FlxEta);
    void  get_interfacial_energy(double Caq);
    void  transport_moments_and_eta( double dt ,std::vector<double> dXML) ;


         void computeFluxes(std::vector<double> &dd) ;
         void computeSourceTerms() ; 
         void setVolumeFraction(std::vector<double> &volFrac) ;

    private:

        ////////////////////// MEMBER FUNCTIONS  /////////////////////
double Sb( double radius, int k);

        ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:
    mom_Aq(ETA_Aq *ETAp,  odtline *odtlp);                ///< Constructor Function
};

#endif

