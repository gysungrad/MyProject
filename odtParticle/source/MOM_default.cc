/**
 * @file MOM_default.cc
 * Source file for class MOM_default
 */

#include "MOM_default.h"
#include "processor.h"
#include "diffuser.h"
#include <cmath>
#include <cstdlib>

using namespace std;

extern processor proc;

///////////////////////////////////////////////////////////////////////////////
/** Constructor function.
 *
*/

MOM_default::MOM_default(odtline *odtlp) {
    
    odtl         = odtlp;
}

///////////////////////////////////////////////////////////////////////////////


/** Update fluxes at cell faces
 *
 *  @param dd    - input  - the grad part of grad phi (multiplier
 */

void MOM_default::computeFluxes(std::vector<double> &dd){

    int i, im, k;

    std::vector<double> visc_f(odtl->ngrd);
    for(i=0; i<odtl->ngrd; i++){
        visc_f[i] = diff->linearInterpToFace(i, *odtl, odtl->molec);
    }
    
    //////////////////// interior faces ////////////////////

    for(i=1, im=0; i<odtl->ngrd; i++, im++) {
        for(k=0; k<odtl->odtP->Imom; k++) {

            diff->flxProp[diff->iptMom+k][i] = -dd[i]*visc_f[i]*(odtl->mom[k][i] - odtl->mom[k][im]);

        }

    } 

    //////////////////// boundary faces ////////////////////

    //-------------- periodic 

    if (diff->odtP->bcType == odtParam::BCTYPE_PERIODIC) { 

        i = 0;
        im = odtl->ngrd - 1;

        for(k=0; k<odtl->odtP->Imom; k++) {

            diff->flxProp[diff->iptMom+k][i]          = -dd[0]*visc_f[0]*( odtl->mom[k][i]-odtl->mom[k][im] );
            diff->flxProp[diff->iptMom+k][odtl->ngrd] = diff->flxProp[diff->iptMom+k][i];

        }
    }

    //----------------- wall on left outflow on right or outflow, etc.

    else if (diff->odtP->bcType == odtParam::BCTYPE_WALL    ||
             diff->odtP->bcType == odtParam::BCTYPE_OUTFLOW ||
             diff->odtP->bcType == odtParam::BCTYPE_INLET   ||
             diff->odtP->bcType == odtParam::BCTYPE_WALLOUTFLOW) { 

        for (int k = 0; k < odtl->odtP->Imom; k++) {
            diff->flxProp[k][0] = 0.0;
            diff->flxProp[k][odtl->ngrd] = 0.0;
        }
    }

    //------------------ unknown bc

    else { 
        *proc.ostrm << endl << "****** ERROR IN DIFFUSER, BCTYPE = " << diff->odtP->bcType << " unknown" << endl;
        exit(0);
    }


} // end function

///////////////////////////////////////////////////////////////////////////////

/** Compute source terms for default moments (currently zero)
 */

void MOM_default::computeSourceTerms() {
    
    return;
}


