/**
 * @file ETA_default.h
 * Header file for class ETA_default 
 */

#ifndef ETA_DEFAULT_H    
#define ETA_DEFAULT_H    

#include <vector>
#include <string>
#include "odtline.h"
#include "odtParam.h"
#include "ETA.h"

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

/** Class for implementing generic eta equations
 *  
 *  @author David O. Lignell
 */

class ETA_default : public ETA {

  ////////////////////// DATA MEMBERS /////////////////////
  
    public:
        
 ////////////////////// MEMBER FUNCTIONS  /////////////////////

        void computeFluxes(std::vector<double> &dd);
        void computeSourceTerms();
        void updateOdtLineVecs(bool updateRho = false) { }
        void updateOdtLineVecs(int & i, bool updateRho = false) { }


 ////////////////// CONSTRUCTOR FUNCTIONS /////////////////////

    public:

        ETA_default(odtline *odtlp); 

};

#endif
