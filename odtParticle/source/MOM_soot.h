/**
 * @file MOM_soot.h
 * Header file for class MOM_soot 
 */

#ifndef MOM_SOOT_H    
#define MOM_SOOT_H    

#include <vector>
#include <string>
#include "odtline.h"
#include "odtParam.h"
#include "MOM.h"

#ifdef DOCANTERA
#ifndef CANTERA21
#include "Cantera.h"
#endif
#include "IdealGasMix.h"
#endif

#include "cantera_shell_functions.h"

class odtline;
class diffuser;

///////////////////////////////////////////////////////////////////////////////

/** Class for implementing particle moment equations
 *  
 *  @author David O. Lignell
 */

class MOM_soot : public MOM {

  ////////////////////// DATA MEMBERS /////////////////////
  
    private:
        
        double                            rho_s;        ///< soot density
        double                            kb;           ///< Boltzman constant
        double                            Na;           ///< Avagadro's number
        double                            Ca;           ///< agglomeration rate constant
        double                            Cmin;         ///< minimum carbons in a soot particle
        double                            MWc;          ///< molecular weight of carbon
        int                               iO2;          ///< these are all references to species

        int                               nmixf;        ///< number of mixf values in look-up table
        int                               nhl;          ///< number of heat loss values in look-up table
        std::vector<double>               hl;           ///< vector of heat loss values found in look up table
        std::vector<double>               mixF;         ///< vector of mixture fraction values found in look up table
        std::vector<std::vector<double> > Y_C2H2_Table; ///< Table of Y_C2H2 at given mixf and hl
                                                        ///  Also includes sensible heat for given mixf

        std::vector<std::vector<double> >      mom_source;   ///<source term for soot  
        //std::vector<std::vector<double> >      gasSp_source; ///<source term for gas species  
  
        

 ////////////////////// MEMBER FUNCTIONS  /////////////////////

    private:

        void resizeVar();
        void getSootSourceTermsPt1mom(double &Ys_source, const int ipos);
        void getSootSourceTermsPt2mom(double &mhat0_source,
                                      double &mhat1_source,
                                      //std::vector<double> &gSp_source,
                                      const int ipos);
        void getSootSourceTermsPt3mom(double &mhat0_source,
                                      double &mhat1_source,
                                      double &mhat2_source,
                                      const int ipos);
        void get_soot_rr(std::vector<double> &soot_rr, double temp, int ipos);
        void getC2H2(const double mf, const double h, double &YC2H2, int posi);
        void readYc2h2Table(std::string fname);
        void interpOnePt(double x1, double y1, double x2, double y2, 
                         double actualx, double &ansy);
        void getIndices(double value, std::vector<double> vec, int &up, int &dn);

    public:
        
        void computeFluxes(std::vector<double> &dd);
        void computeSourceTerms();

        void setVolumeFraction(std::vector<double> &volFrac);
          


 ////////////////// CONSTRUCTOR FUNCTIONS /////////////////////

    public:

        MOM_soot(odtline *odtlp, std::string ftable); 

};

#endif
