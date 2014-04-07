/**
 * @file ETA_Aq.h
 * Header file for class ETA 
 */

#ifndef ETA_Aq_H    
#define ETA_Aq_H    
#include "odtline.h"
#include <vector>
#include <string>

///////////////////////////////////////////////////////////////////////////////

/** Class for implementing additional transported variables (e.g. mixture fraction)
 *  
 *  @author David O. Lignell
 */

class ETA_Aq : public  ETA  {

  ////////////////////// DATA MEMBERS /////////////////////

    public:

    int                               nElem;          ///< # of elements
    int                               nEqSp;          ///< # of equilibrium species
    int                               nAqSp;          ///< # of aqueous species
    int                               nEqSolnVars;    ///< # of variables in the nonlinear solve
    int                               nEqRxns;        ///< # of equilibrium reactions
    std::vector<std::string>          aqSpNames;      ///< species names
    std::vector<int>                  z_charges;      ///< ionic charges on all species
    std::vector<double>               gamma;          ///< molar activity coefficients 
    std::vector<std::vector<double> > A_rxnCoef;      ///< Coefficient matrix for eq rxns
    std::vector<std::vector<double> > Q_rxnCoef;      ///< Q of QR decomp of A
    std::vector<std::vector<double> > R_rxnCoef;      ///< R of QR decomp of A
    std::vector<double>               Klog10;         ///< log10(Keq)
    std::vector<double>               Conc;           ///< concetrations (moles/L) of all species
    std::vector<double>               mwel;           ///< elemental molecular weights: C, H, O, Ca, Na, Cl
    std::vector<double>               mwsp;           ///< elemental molecular weights: C, H, O, Ca, Na, Cl
    std::vector<std::vector<double> > Yspel;          ///< mass frac of element in species [sp][el]
    std::vector<std::vector<double> > allConc;        ///< concentrations (moles/L) everywhere [sp][loc]
    std::vector<std::vector<double> > J;              ///< Analytical Jacobian by Derek
    std::vector<std::vector<int> >    moleElSp;       ///< moles of element in species [sp][el]
    double                         Total;             ///< Total Density, used in the sawada case using seans density
    double                            Cc;             ///< concentration of carbon (moles/L)
    double                            Ch;             ///< concentration of hydrogen (moles/L)
    double                            Co;             ///< concentration of oxygen (moles/L)
    double                            Cca;            ///< concentration of calcium (moles/L)
    double                            Cna;            ///< concentration of sodium (moles/L)
    double                            Ccl;            ///< concentration of chlorine (moles/L)
    bool                            Llewis;           ///< True for preset Lewis number (defined in ETA)
    double                              TA;           ///< Temperature of current stream (placeholder)
    double                            T0;             ///< temperature of stream 0 (K)
    double                            h0;             ///< enthalpy of stream 0 (J/kg)
    double                            cp;             ///< water heat capacity (J/kg*K)
    std::vector<double>              xOld;            ///< equilibrium solution vector
    std::vector<double>              xOldi;           ///< old solution to previous time steps in the cell with index zero
    std::vector<double>               ISc;            ///< Coefficients for Ionci stength approximation 
    std::vector<double>               CaCO3;          ///< CaCO3 concentration (mol/L) 
    std::vector<double>      CaCO3_gamma;           ///< CaCO3 concentration (mol/L) 
    std::vector<double>                b;             ///< RHS function for vector oringinally used in linear solve which was in the newton solve
    double                             a;             ///< RHS function dummy variable, constant within the newton solve
    double                             B;             ///< RHS function dummy variable, constant within the newton solve
    double                             c;             ////< RHS function dummy variable, constant within the newton solve
    double                             d;             ///< RHS function dummy variable, constant within the newton solve
    double                             e;             ///< RHS function dummy variable, constant within the newton solve
    double                             f;             ///< RHS function dummy variable, constant within the newton solve


   

    ////////////////////// MEMBER FUNCTIONS  /////////////////////

    public:
    void Maxwell_Stefan_Diffusion(std::vector<double> &dx2, std::vector<std::vector<double> > &allJi, std::vector<double> &Temp); 
    void setEta(std::vector<double> &mixf, std::vector<std::vector<double> > &eta);
    void getTempVecFromH(std::vector<std::vector<double> > &eta, std::vector<double> &T);
    void getEtaDiffusiveFluxes(std::vector<std::vector<double> > &eta, double rho,
                               std::vector<double> &dx, 
                               std::vector<std::vector<double> > &flxEta,
                               bool Lmom);
    void   equilibrateAll(std::vector<std::vector<double> > &eta, double rho);
    bool   equilibrate(std::vector<double> &Yelem, double rho, double T, bool Lfirst);
    void   Analytical_Jac(std::vector<double> x, double T) ;
    double getTgivenH(double h);
    void QR(std::vector<std::vector<double> > &a,std::vector<double> &b,std::vector<double> &x);
    void Correct_rho(double &rho, double Yna, double Ycl);



       void computeFluxes(vector<double> &dd) ;
       void computeSourceTerms() ; 
       void updateOdtLineVecs(bool updateRho = false) ;
       void updateOdtLineVecs(int & i, bool updateRho = false);

    private:
    void   Compute_RHS_Constants();
    void   elemMassFrac_to_elemConc(std::vector<double> &Yelem, double rho, double T );
    void   setActivityCoefficients(double T, double rho);
    void   setActivityCoefficients2(double T,double x);
    void   equilFrhs(std::vector<double> &x, std::vector<double> &y);
    void   equilFrhsw(std::vector<double> &x, std::vector<double> &y);
    void   equilJac(std::vector<double> x1, std::vector<std::vector<double> > &Jac);
    void   setKlog10(double T);
    //double getTgivenH(double h);
    void   getDiffusivities(std::vector<double> &T, 
                            std::vector<std::vector<double> > &Di);
    double getViscosity(double T);
    double getThermalConductivity(double T);

    ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:
        
    ETA_Aq(odtline *x);
    ETA_Aq();
           
};

#endif

