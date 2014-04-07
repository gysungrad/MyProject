/*
 * @file ETA.h
 * Header file for class eta
 */

#ifndef ETA_H
#define ETA_H

///////////////////////////////////////////////////////////////////////////////

/*
 * * Class for implementing additional transported variables (e.g. mixture fraction)
 *
 *  @author David O. Lignell
 */

#include <vector>
#include <iostream>
#include "odtParam.h"
#include "processor.h"
#include "radiation.h"
#include "table.h"

using namespace std;

class odtline;
class diffuser;

class ETA {
    
    ////////////////////// DATA MEMBERS /////////////////////

    public:

        table                            *luTable;      ///< lookup table
        odtline                          *odtl;         ///< A backpointer to the odtline
        diffuser                         *diff;         ///< pointer to the diffuser object for convenience

        ////////////////////// MEMBER FUNCTIONS  /////////////////////

	virtual ~ETA() {}
        virtual void computeFluxes(vector<double> &dd) = 0;
        virtual void computeSourceTerms() = 0; 
        virtual void updateOdtLineVecs(bool updateRho = false) = 0;
        virtual void updateOdtLineVecs(int & i, bool updateRho = false)=0;

};
#endif
