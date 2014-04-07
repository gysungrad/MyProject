/**
 * @file MOM_default.h
 * Header file for class MOM_default 
 */

#ifndef MOM_DEFAULT_H    
#define MOM_DEFAULT_H    

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

class MOM_default : public MOM {

  ////////////////////// DATA MEMBERS /////////////////////
  
    public:
        
 ////////////////////// MEMBER FUNCTIONS  /////////////////////

        void computeFluxes(std::vector<double> &dd);
        void computeSourceTerms();
        void setVolumeFraction(std::vector<double> &volFrac){};


 ////////////////// CONSTRUCTOR FUNCTIONS /////////////////////

    public:

        MOM_default(odtline *odtlp); 

};

#endif
