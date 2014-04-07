/*
 * @file MOM.h
 * Header file for class eta
 */

#ifndef MOM_H
#define MOM_H

#include <vector>
#include <iostream>
#include "odtParam.h"
#include "processor.h"
#include "radiation.h"
#include "table.h"

///////////////////////////////////////////////////////////////////////////////

/*
 * * Class for implementing additional transported variables (e.g. particle moments)
 *
 *  @author David O. Lignell
 */


class odtline;
class diffuser;
class MOM {
    
    ////////////////////// DATA MEMBERS /////////////////////

    public:

        odtline                          *odtl;         ///< A backpointer to the odtline
        diffuser                         *diff;         ///< pointer to the diffuser object for convenience

        ////////////////////// MEMBER FUNCTIONS  /////////////////////

        virtual void computeFluxes(std::vector<double> &dd) = 0;
        virtual void computeSourceTerms() = 0; 
        virtual void setVolumeFraction(std::vector<double> &volFrac) = 0;

};
#endif
