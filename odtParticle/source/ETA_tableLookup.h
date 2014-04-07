/*
 * @file ETA_tableLookup.h
 * Header file for class ETA_tableLooup
 */

#ifndef ETA_TABLELOOKUP_H
#define ETA_TableLookup_H

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
class ETA_tableLookup : public ETA{
    
    ////////////////////// DATA MEMBERS /////////////////////

    public:



        ////////////////////// MEMBER FUNCTIONS  /////////////////////

        void computeFluxes(std::vector<double> &dd);
        void computeSourceTerms(); 

        void updateOdtLineVecs(bool updateRho = false);
        void updateOdtLineVecs(int & i, bool updateRho = false);

        ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

        ETA_tableLookup(odtline *newOdtl, table *luTable_);
};
#endif
