/**
 * @file diffuser.cc
 * Source file for class diffuser
 */

#include "diffuser.h"
#include "processor.h"
#include "assert.h"
#include <string>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#include <vector>

using namespace std;

extern processor proc;

///////////////////////////////////////////////////////////////////////////////
// Declare the global function prototype so it can be used in this source file

void getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rr);

///////////////////////////////////////////////////////////////////////////////

/** Constructor function to use.
 * 
 * Class operates on lines performing diffusion step of ODT
 *
 * @param odtlp           \input odtline object, set pointer with
 * @param odtpp           \input parameters object, set pointer with
 * @param partp           \input particles object, set pointer with
 */

diffuser::diffuser(odtline *odtlp, odtParam &odtpp,particles* partp) :

//---------- Constructor initialization list

brxr(odtlp, &odtpp),
dmpr(odtpp.tEnd),
rad()
 {
    //------------------------------------------------------------------------

    odtP = &odtpp;
    odtl = odtlp;
    part = partp;

    neqns = odtl->nprops;
    k1rhs = vector<vector<double> >(neqns, vector<double>(odtl->ngrd, 0.0));
    rhsTrn = vector<vector<double> >(neqns, vector<double>(odtl->ngrd, 0.0));
    rhsPartSrc = vector<vector<double> >(neqns, vector<double>(odtl->ngrd, 0.0));
    rhsSrc = vector<vector<double> >(neqns, vector<double>(odtl->ngrd, 0.0));
    flxProp = vector<vector<double> >(neqns, vector<double>(odtl->ngrd, 0.0));


    iptUvel = odtl->LhasVel ? 0 : -1;
    iptVvel = odtl->LhasVel ? 1 : -1;
    iptWvel = odtl->LhasVel ? 2 : -1;

    iptEnth = -1;
    iptYspc = -1;
    iptMom = -1;
    iptEta = -1;

    dxML = vector<double> (odtl->ngrd);
    
    timeFS = 0.0;
    
    int iptDumb = 3; // dummy for setting ipt*
    
    if (odtP->LhasTemp) {
        iptTemp = 3;
        iptDumb++;
    }
    
    if (odtP->Lrxn) {
        iptEnth = odtl->LhasVel ? 3 : 0;
        iptYspc = odtl->LhasVel ? 4 : 1;
        iptDumb = iptYspc + odtl->nspc;
        visc_f = vector<double>(odtl->ngrdf, odtP->visc_0);
        lambda_f = vector<double>(odtl->ngrdf, 0.02);
        DmixYs_f = vector<vector<double> >(odtl->nspc, vector<double>(odtl->ngrdf, odtP->visc_0 / odtP->rho_0));
        hsp = vector<vector<double> >(odtl->nspc, vector<double>(odtl->ngrd, 0.0));
        rrSpc = vector<vector<double> >(odtl->nspc, vector<double>(odtl->ngrd, 0.0));

        iN2 = odtl->gas->speciesIndex("N2");
        iN2 = (iN2 > 0) ? iN2 : odtl->gas->speciesIndex("n2");
        assert(iN2 > 0); // error if N2 not in system dolcheck: generalize this
    } else {
        visc_f = vector<double>(odtl->ngrdf, odtP->visc_0);
    }

    if (odtP->Iradiation > 0) {
        rad = radiation(odtpp, odtl->gas);
        radSource_G = vector<double>(odtl->ngrd, 0.0);
        if (part)
            radSource_P = vector<double>(part->nPart, 0.0);
    }

    if (odtP->Ieta) {
        iptEta = iptDumb;
        iptDumb = iptEta + odtl->neta;
    }
    if (odtP->Imom) {
        iptMom = iptDumb;
        iptDumb = iptMom + odtl->nmom;
    }
    Gvel.resize(odtl->ngrdf);

    if (odtP->Iparticles) {
        int isize = (odtP->Lprxn) ? 6 : 4;
        P1rhs = vector<vector<double> >(isize, vector<double>(part->nPart, 0.0));
    }

    // gas Source------------------------
    if (odtP->Iparticles) {
        gSource = new vector<vector<double> >(part->iPtMass + 1);
        for (int kk = 0; kk < part->iPtMass + 1; kk++) { //iPtMass+1 for size because mass is the last vector
            (*gSource)[kk] = vector<double>(odtl->ngrd, 0);
        }
        //----------- This resize and for loop is similar to when odtP->Lprxn

        lambda_f.resize(odtl->ngrdf);

        for (int ii = 0; ii < odtl->ngrdf; ii++) {
            if (odtP->ItableLookup)
                lambda_f[ii] = odtlp->etaTools->luTable->getValAtGridPoint(ii, "LAMBDA_INDEX");
            else
                lambda_f[ii] = odtl->tran->thermalConductivity();
        }

    } else {
        gSource = 0;
    }

#ifdef IMPLICIT
    A      = vector<vector<double> >(neqns, vector<double>(odtl->ngrd));
    B      = vector<vector<double> >(neqns, vector<double>(odtl->ngrd));
    C      = vector<vector<double> >(neqns, vector<double>(odtl->ngrd));
    rhsImp = vector<vector<double> >(neqns, vector<double>(odtl->ngrd));
#endif
}


//for subdomain decomposition
void diffuser::set_diffuser(odtline *odtlp)
{
    odtl = odtlp;
//    brxr(odtlp, &odtpp);
//    rad();


    resetVarSizes();
    

}


///////////////////////////////////////////////////////////////////////////////

/** User interface and driver.
 * Integrates (diffuses) the solution from time to time+dtDiffuse
 * Simple explicit Euler integration, but the class should be easily
 * expandable for more complex solvers.
 *
 * @param dtDiffuse \input time duration to diffuse the line.
 * @param time      \input current run time in simulation realization.
 * @param odtStats  \inout stats object to update.
 * @param iStat     \input number of iStat period. Needed for stats-object update
 */

void diffuser::diffuseOdtLine(double dtDiffuse, double time, stats &odtStats, const int &iStat) {

#ifdef NEWSTATS
    uBulk = odtStats.getUBulk(0, iStat);
#else
    uBulk = 1;
#endif
    
    resetVarSizes();
    setGridSize();
    setDiffusivities_HSP_RR();
    setTimeStep();
#ifdef NEWSTATS
    odtStats.initBStats(*odtl); // check for new phases
#endif
#ifdef IMPLICIT
    double dtStep = dtDiffuse;
#else
    double dtStep = dtStepCFL;
#endif

    *proc.ostrm << "\n# DIFFUSING: dtDiffuse, dtStep, M, N = " << dtDiffuse << "  " << dtStep << " " << dtDiffuse/dtStep << " " << odtl->ngrd;
    (*proc.ostrm).flush();

    ///////////////////
    /*
     * Approximate dtAdvecMax: max time step allowed in opposedJets and adjust dtStepCFL if needed.
     * Only an approximation because it is done with velocity before updating to new time step.
     */
    if (odtl->probType == 4) {

        // set variables
        int ngrd = odtl->ngrd;
        double dtAdvec;
        double dtAdvecMax = 100.0;
        vector<double> vf(ngrd + 1, 0.0);

        // calculate velocity at faces
        vf[0] = odtP->uBClo;
        for (int i = 1; i < ngrd; i++)
            vf[i] = odtl->uvel[i - 1] + (odtl->posf[i] - odtl->pos[i - 1]) / (odtl->pos[i] - odtl->pos[i - 1]) * (odtl->uvel[i] - odtl->uvel[i - 1]);
        vf[ngrd] = odtP->uBChi;

        for (int i = 0; i < ngrd; i++) {
            dtAdvec = (odtl->posf[i + 1] - odtl->posf[i]) / (vf[i] - vf[i + 1]);
            if (dtAdvec < dtAdvecMax && dtAdvec > 0.0)
                dtAdvecMax = dtAdvec;
        }

        if (dtAdvecMax < dtStepCFL) {
            dtStepCFL = dtAdvecMax / 4.0; //factor 2 is only for safety
            dtStep = dtStepCFL;
        }

        *proc.ostrm << endl << "dtAdvecMax, dtStepCFL, dtStep :" << dtAdvecMax << "   " << dtStepCFL << "   " << dtStep << endl;
    }
    
#ifdef IMPLICIT // low accuracy, first part of 2nd order integration
  #ifdef NEWSTATS
    if (!odtP->Llem) odtStats.BStats(*odtl, dtDiffuse/2.0, iStat);
  #endif
#else
  #ifdef NEWSTATS
    #ifndef NEWSTATSH
    if (!odtP->Llem) odtStats.BStats(*odtl, dtDiffuse/2.0, iStat);
    #endif
  #endif
#endif
    
    for (double timeD = 0.0; timeD < dtDiffuse; timeD += dtStep) {
        adaptGridsIfNeeded();

        if (odtP->Lspatial) {
            setTimeStep();
            dtStep = dtStepCFL;
        }

        if (odtP->Lprxn)
            part->getRelativeVel(time + timeD);
        
        if(odtP->Iparticles) {
            if(part->PeddyType == 2 || part->PeddyType == 3) {
                checkActiveEddyeffect(time+timeD);
                //             getEddyUWvel(time+timeD);
                updatePartIxn(time+timeD);
            }
        }

#ifndef IMPLICIT
        if (timeD + dtStepCFL > dtDiffuse) dtStep = dtDiffuse - timeD;
#endif

	    diffuseSingleStep(dtStep, time+timeD);
        // diffuseSingleStep(dtStep);

        if (odtP->Iparticles && part->Lhistories) {
            if (odtP->Lrxn) {
                odtl->setMixfVec();
                odtl->setTempVec();
                odtl->setChiVec();
                odtl->setDiffZVec();
            }
            part->storeHistories(time + timeD);
        }

        if(!odtP->Lsubdomain) outputAtSetTimes(time, timeD);

#ifndef IMPLICIT
  #ifdef NEWSTATSH // high accuracy, first order integration
        if (!odtP->Llem) odtStats.BStats(*odtl, dtStep, iStat);
  #endif
  #ifndef NEWSTATS // high accuracy
        //if (!odtP->Llem) odtStats.updateMeans(*odtl, time+timeD, dtStep);
  #endif
#endif
    }
    
#ifdef IMPLICIT // lower accuracy, second part of 2nd order integration
  #ifdef NEWSTATS
    if (!odtP->Llem) odtStats.BStats(*odtl, dtDiffuse/2.0, iStat);
  #else
    if (!odtP->Llem) odtStats.updateMeans(*odtl, time + dtDiffuse, dtDiffuse);
  #endif
#else
  #ifdef NEWSTATS
    #ifndef NEWSTATSH
    if (!odtP->Llem) odtStats.BStats(*odtl, dtDiffuse/2.0, iStat);
    #endif
  #endif
  #ifndef NEWSTATS // lower accuracy
    if (!odtP->Llem) odtStats.updateMeans(*odtl, time + dtDiffuse, dtDiffuse);
  #endif
#endif
    timeFS += dtDiffuse;
    
    if(odtP->Lsubdomain)
      return;

    ////////////////////// chop outflow grid

    if (odtP->Lrxn || odtP->IetaType == odtP->IETATYPE_TABLELOOKUP || odtP->IetaType == odtP->IETATYPE_DEFAULT)
        if (odtP->bcType == 3 || odtP->bcType == 4 || odtP->bcType == 5)
            chopOutflowGrid();

    ////////////////////// update domainLength if needed

    if (odtl->Ldomain != odtP->domainLength) // todo: get rid of this
        odtP->resetDomainLength(odtl->Ldomain);
}

///////////////////////////////////////////////////////////////////////////////

void diffuser::getMomTotal(odtline &line) {

    line.momGas   = 0.;
    line.momPart  = 0.;
    line.momTotal = 0.;
    line.ergGas   = 0.;
    line.ergPart  = 0.;
    line.ergTotal = 0.;

    for(int i=0; i<line.ngrd; i++) {
        double dxML = line.posf[i + 1] - line.posf[i];
        line.momGas = line.momGas + line.rho[i] * dxML * (line.uvel[i]+line.vvel[i]+line.wvel[i]); 
        line.ergGas = line.ergGas + 0.5 * line.rho[i] * dxML * (line.uvel[i] * line.uvel[i] + line.vvel[i] * line.vvel[i] + line.wvel[i] * line.wvel[i]);
    }

    if(odtP->Iparticles) {
        for (int i = 0; i < part->nPart; i++) {
            double pMass = part->nInPseudoPart[i]*part->pDens0[i]*4./3.*3.14159*part->pRadi[i]*part->pRadi[i]*part->pRadi[i];
            line.momPart = line.momPart + pMass * (part->uvel[i] + part->vvel[i] + part->wvel[i]);
            line.ergPart = line.ergPart + 0.5 * pMass * (part->uvel[i] * part->uvel[i] + part->vvel[i] * part->vvel[i] + part->wvel[i] * part->wvel[i]);
        }

        line.momTotal = line.momGas + line.momPart;
        line.ergTotal = line.ergGas + line.ergPart;
    }
//     cout << endl << "~~~~~~~~~~~" << endl;
//     cout << endl << "line grid " << line.ngrd << endl;
     cout << endl << fixed << setprecision(12) << "mom gas " << line.momGas << endl;
     cout << endl << fixed << setprecision(12) << "mom Part " << line.momPart << endl;
     cout << endl << fixed << setprecision(12) << "mom total " << line.momTotal << endl;
//     cout << endl << fixed << setprecision(12) << "erg gas " << line.ergGas << endl;
//     cout << endl << fixed << setprecision(12) << "erg part " << line.ergPart << endl;
     cout << endl << fixed << setprecision(12) << "erg total " << line.ergTotal << endl;
}

///////////////////////////////////////////////////////////////////////////////

///subdomain decomposition

void diffuser::diffuseOdtLine(double dtDiffuse, double time) {


    resetVarSizes();
    setGridSize();
    setDiffusivities_HSP_RR();
    setTimeStep();

    double dtStep = dtStepCFL;


    *proc.ostrm << "\n# DIFFUSING: dtDiffuse, dtStep, M = " << dtDiffuse << "  " << dtStep << " " << dtDiffuse/dtStep;
    (*proc.ostrm).flush();

    ///////////////////
    /*
     * Approximate dtAdvecMax: max time step allowed in opposedJets and adjust dtStepCFL if needed.
     * Only an approximation because it is done with velocity before updating to new time step.
     */
    if (odtl->probType == 4) {

        // set variables
        int ngrd = odtl->ngrd;
        double dtAdvec;
        double dtAdvecMax = 100.0;
        vector<double> vf(ngrd + 1, 0.0);

        // calculate velocity at faces
        vf[0] = odtP->uBClo;
        for (int i = 1; i < ngrd; i++)
            vf[i] = odtl->uvel[i - 1] + (odtl->posf[i] - odtl->pos[i - 1]) / (odtl->pos[i] - odtl->pos[i - 1]) * (odtl->uvel[i] - odtl->uvel[i - 1]);
        vf[ngrd] = odtP->uBChi;

        for (int i = 0; i < ngrd; i++) {
            dtAdvec = (odtl->posf[i + 1] - odtl->posf[i]) / (vf[i] - vf[i + 1]);
            if (dtAdvec < dtAdvecMax && dtAdvec > 0.0)
                dtAdvecMax = dtAdvec;
        }

        if (dtAdvecMax < dtStepCFL) {
            dtStepCFL = dtAdvecMax / 4.0; //factor 2 is only for safety
            dtStep = dtStepCFL;
        }

        cout << endl << "dtAdvecMax, dtStepCFL, dtStep :" << dtAdvecMax << "   " << dtStepCFL << "   " << dtStep << endl;
    }
    for (double timeD = 0.0; timeD <= dtDiffuse; timeD += dtStepCFL) {

        adaptGridsIfNeeded();

        if (odtP->Lspatial) {
            setTimeStep();
            dtStep = dtStepCFL;
        }

        if (odtP->Lprxn)
            part->getRelativeVel(time + timeD);
	
	
        if (timeD + dtStepCFL > dtDiffuse) dtStep = dtDiffuse - timeD;
	
        if ( dtStep > 0 )  //will get zero on very first requested dump time
        // diffuseSingleStep(dtStep);
	        diffuseSingleStep(dtStep, time+timeD);
	
	
        if (odtP->Iparticles && part->Lhistories) {
            part->storeHistories(time + timeD);
        }

        if(!odtP->Lsubdomain) outputAtSetTimes(time, timeD);
	
	setGridSize();
	
	
    }



    if(odtP->Lsubdomain)
      return;

    ////////////////////// chop outflow grid
    
    
    if (odtP->Lrxn || odtP->IetaType == odtP->IETATYPE_TABLELOOKUP || odtP->IetaType == odtP->IETATYPE_DEFAULT)
        if (odtP->bcType == 3 || odtP->bcType == 4 || odtP->bcType == 5)
            chopOutflowGrid();

    ////////////////////// update domainLength if needed

    if (odtl->Ldomain != odtP->domainLength) // todo: get rid of this
        odtP->resetDomainLength(odtl->Ldomain);
}

///////////////////////////////////////////////////////////////////////////////

/** Output results at times specified in the dumpTimes.inp file
 * @param time    \input current run time in simulation realization.
 * @param timeD   \input current dtCFL step
 */

void diffuser::outputAtSetTimes(double time, double timeD) {

    if (dmpr.n_time > 0 && dmpr.checkDumpTime(time + timeD)) {
        *proc.ostrm << "\n# WRITING DMP FILE " << proc.dataDir + "dmp*" + dmpr.currentTimeFile << endl;

        odtline bkp_odtl = *odtl;
        particles *bkp_part;
        if (odtP->Iparticles) {
            bkp_part = new particles(*part);
        }
        
        if (odtP->Lrxn || odtP->IetaType == odtP->IETATYPE_TABLELOOKUP || odtP->IetaType == odtP->IETATYPE_DEFAULT)
            if (odtP->bcType == 3 || odtP->bcType == 4 || odtP->bcType == 5)
                chopOutflowGrid();
        
        odtl->outputProperties(proc.dataDir + "dmp_odtl_" + dmpr.currentTimeFile);
        if (odtP->Iparticles) {
            part->outputProperties(proc.dataDir + "dmp_part_" + dmpr.currentTimeFile);
            *part = *bkp_part;
            delete bkp_part;
        }

        *odtl = bkp_odtl;
        
        resetVarSizes();
        setGridSize();
    }
}

///////////////////////////////////////////////////////////////////////////////

/**Function resets the sizes of arrays, which change due to mesh adaption.
 * Does not re-initialize these arrays though.
 */

void diffuser::resetVarSizes() {

    k1rhs.resize(neqns); // redundant
    rhsTrn.resize(neqns);
    rhsSrc.resize(neqns);
    rhsPartSrc.resize(neqns);
    flxProp.resize(neqns);
    for (int i = 0; i < neqns; i++) {
        k1rhs[i].resize(odtl->ngrd);
        rhsTrn[i].resize(odtl->ngrd);
        rhsSrc[i].resize(odtl->ngrd);
        rhsPartSrc[i].resize(odtl->ngrd);
        flxProp[i].resize(odtl->ngrdf);
    }

    dxML.resize(odtl->ngrd);

    if (odtP->Iradiation > 0) {
        radSource_G = vector<double>(odtl->ngrd, 0.0);
        if (part)
            radSource_P = vector<double>(part->nPart, 0.0);
    }

    if (odtP->Lrxn) {

        odtl->temp.resize(odtl->ngrd); // auxiliary variable, set size here

        visc_f.resize(odtl->ngrdf);
        lambda_f.resize(odtl->ngrdf);
        for (int k = 0; k < odtl->nspc; k++) {
            DmixYs_f[k].resize(odtl->ngrdf);
            hsp[k].resize(odtl->ngrdf);
            rrSpc[k].resize(odtl->ngrdf);
        }

        invGamma.resize(odtl->ngrdf);
        mmw.resize(odtl->ngrdf); // not needed
    } else {
        visc_f.resize(odtl->ngrdf, odtP->visc_0);
        if (odtP->LhasTemp){
            lambda_f.resize(odtl->ngrdf);
        }
    }
    
    Gvel.resize(odtl->ngrdf);
    
    // todo: resize P1rhs if necessary, but currently the number of particles does not change

}
///////////////////////////////////////////////////////////////////////////////

/** Set the timestep based on the cfl criteria (smallest cell).
 * Make sure dx is up to date before calling this function
 * Step is based on a Fourier limit on the momentum field:
 *                         \cond
 * nu*dt/dx/dx <= 1/2.     \endcond                                     \fun{\nu*dt/dx/dx \le 1/2}.
 *
 * Here, \fun{\nu} is taken as \fun{\nu_{max}} from the odtline
 * and should be (for these purposes) the same as if use the comb line.
 */

void diffuser::setTimeStep() {

    if (odtP->Lspatial) {
        double velMin = *min_element(odtl->uvel.begin(), odtl->uvel.end());
        if (velMin <= 0.0) {
            cout << "\n**** Error setting velMin for Lspatial, negative or zero value" << endl;
            cout << "\nvelMin = " << velMin << endl;
            exit(0);
        }
    }

    double coef = 0.0;
    double dmb;
    for (int i = 0; i < odtl->ngrd; i++) {
        dmb = odtl->molec[i] / odtl->rho[i] / dxML[i] / dxML[i];
        if (odtP->Lspatial)
            dmb /= odtl->uvel[i];
        if (dmb > coef)
            coef = dmb;
    }

    dtStepCFL = 1.0 / coef;

    if (odtP->Lrxn) {
        for (int k = 0; k < odtl->nspc; k++) {
            coef = 0.0;
            for (int i = 0; i < odtl->ngrd; i++) {
                dmb = DmixYs_f[k][i] / dxML[i] / dxML[i];
                if (odtP->Lspatial)
                    dmb /= odtl->uvel[i];
                if (dmb > coef)
                    coef = dmb;
            }
            dtStepCFL = min(dtStepCFL, 1.0 / coef);
        }
    }

    if (odtP->LhasTemp) {
        coef = 0.0;
        for(int i=0; i<odtl->ngrd; i++) {
            dmb = lambda_f[i] / dxML[i] / dxML[i];
            if(dmb > coef)
                coef = dmb;
        }
        dtStepCFL = min(dtStepCFL, 1.0 /coef);
    }

    dtStepCFL *= 0.5; // diffusive step limit

    //--------- particle step limit

    if (odtP->Iparticles) {
        for (int i = 0; i < part->nPart; i++) {
            if (!part->pActive[i]) continue;
            part->set_TauP(i);
            // dmb = part->TauP[i] / part->f[i];
            dmb = part->TauP[i];
            dtStepCFL = min(dtStepCFL, dmb);
        }
    }

    //--------- Apply cfl factor

    dtStepCFL *= odtP->diffCFL;

}

///////////////////////////////////////////////////////////////////////////////

/**Set the cell sizes vector: dx */

void diffuser::setGridSize() {

    for (int i = 0; i < odtl->ngrd; i++) {
        dxML[i] = odtl->posf[i + 1] - odtl->posf[i];
    }

}

///////////////////////////////////////////////////////////////////////////////

/**Compute face values of viscosity, thermal conductivity, species diffusivities.
 * Also compute cell centered values of species enthalpies.
 * Also compute cell centered reaction rates.
 * Also compute cell centered gamma.
 * todo: these yi, Di, hh, rr can be consolodated into one dummy
 * todo: create an odtline function that resets the gas object state
 */

void diffuser::setDiffusivities_HSP_RR() {


    if (!odtP->Lrxn) {
        if (odtP->LconstProp) {
            visc_f = vector<double>(odtl->ngrdf, odtP->visc_0);
            if (odtP->LhasTemp){
                lambda_f = vector<double>(odtl->ngrdf, odtP->lambda_0);
            }
        } else {
            if (odtP->ItableLookup){
                for (int i = 0; i < odtl->ngrd; i++) {
                    lambda_f[i] = odtl->lambda[i];
                }
                interpDiffusionCoeffToFaces(*odtl, lambda_f);
                odtl->etaTools->updateOdtLineVecs();
            }
            for (int i = 0; i < odtl->ngrd; i++)
                visc_f[i] = odtl->molec[i];
            if (odtP->LhasTemp)
                for (int i = 0; i < odtl->ngrd; i++)
                    lambda_f[i] = odtl->lambda[i];
            
            interpDiffusionCoeffToFaces(*odtl, visc_f);
            if (odtP->LhasTemp) {
                interpDiffusionCoeffToFaces(*odtl, lambda_f);
            }
        }
        if (!odtP->Lfalko)
            setRadHeatSources();
        return;
    }
    
    if (odtl->probType == 4) { //for opposedJets, calculate also Le_s, Pr, Sc_s

        //resize vectors
        odtl->Pr.resize(odtl->ngrdf, 0.0);
        odtl->Le_s.resize(odtl->nspc);
        odtl->Sc_s.resize(odtl->nspc);

        for (int k = 0; k < odtl->nspc; k++) {
            odtl->Le_s[k].resize(odtl->ngrdf, 0.0);
            odtl->Sc_s[k].resize(odtl->ngrdf, 0.0);
        }
    }

    //--------------------------------------------------------------------------------

    double GasConstant = 8314.47215; // J/kmol*K

    //---------- get diffusion coefficients at cell centers

    vector<double> yi(odtl->nspc); // working arrays
    vector<double> Di(odtl->nspc);
    vector<double> hh(odtl->nspc);
    vector<double> rr(odtl->nspc);

    // dolcheck this looks redundant, delete it
    visc_f.resize(odtl->ngrdf); // visc_f is on odtl below

    double cp; // for unity Le

    for (int i = 0; i < odtl->ngrd; i++) {
        odtl->getYspVecAtPt(i, yi);
        odtl->gas->setState_PY(odtl->pres, &yi[0]);
        odtl->gas->setState_HP(odtl->enth[i], odtl->pres, 1.E-10);
        //odtl->gas->setState_HP(odtl->enth[i], odtl->pres);
        odtl->temp[i] = odtl->gas->temperature();
        invGamma[i] = 1.0 - GasConstant / odtl->gas->cp_mole();
        odtl->gas->getEnthalpy_RT(&hh[0]);
#ifdef PROBLEMSPECIFICRR
        getProblemSpecificRR(odtl->gas->density(), odtl->temp[i], odtl->pres, &yi[0], &rr[0]);
#else
        odtl->gas->getNetProductionRates(&rr[0]);
#endif
        cp = odtl->gas->cp_mass(); // for unity Le
        visc_f[i] = odtl->tran->viscosity();
        //odtl->molec[i] = visc_f[i];
        lambda_f[i] = odtl->tran->thermalConductivity();
        odtl->tran->getMixDiffCoeffs(&Di[0]);
        for (int k = 0; k < odtl->nspc; k++) {
            //DmixYs_f[k][i] = lambda_f[i]/odtl->gas->density()/cp;   // for unity Le doldb
            DmixYs_f[k][i] = Di[k];
            hsp[k][i] = hh[k] * odtl->temp[i] * GasConstant /
                    odtl->gas->molecularWeight(k);
            rrSpc[k][i] = rr[k] * odtl->gas->molecularWeight(k);
        }

        if (odtl->probType == 4) { //calculate also Pr, Le_s, Sc_s
            // Pr
            odtl->Pr[i] = odtl->gas->density() * visc_f[i] * cp / lambda_f[i];

            // Le_s, Sc_s
            for (int k = 0; k < odtl->nspc; k++) {
                odtl->Le_s[k][i] = lambda_f[i] / cp / DmixYs_f[k][i] / odtl->gas->density();
                odtl->Sc_s[k][i] = odtl->Le_s[k][i] * odtl->Pr[i];
            }
        }
    }

    setRadHeatSources();

    //---------- interpolate to faces harmonically

    interpDiffusionCoeffToFaces(*odtl, visc_f);
    interpDiffusionCoeffToFaces(*odtl, lambda_f);
    for (int k = 0; k < odtl->nspc; k++)
        interpDiffusionCoeffToFaces(*odtl, DmixYs_f[k]);
    
}


///////////////////////////////////////////////////////////////////////////////

/** Forward Euler integration. 
 *
 *  @param dtStep \input single explicit step size as part of integration.
 */

void diffuser::diffuseSingleStep(double &dtStep, double time) {

    ////////////////////// save the original lines
    Class_dtStep = dtStep; //DH

    odtline *odtl0 = 0;
    particles *part0 = 0;


    if (odtP->Lstrang) { // strang splitting of diffusion and chemistry in species equation with implicit diffusion and chemistry (should be used for premixed to run stable)

        strangSplitting(dtStep, time);

        if (odtl0) delete odtl0;
        if (part0) delete part0;
        return;
    }

    if (odtP->LsecondOrder) {
        odtl0 = new odtline(*odtl);
        if (odtP->Iparticles)
            part0 = new particles(*part);
    }

    ////////////////////////////////////////////////////////////////////////////////////////// 
    //////////////////////////////////////////////////////////////////////////////////////////

    rhsf(k1rhs, P1rhs, dtStep, time);

    //-----------save uOld for spatial advancement

    vector<double> uOld;
    if (odtP->Lspatial)
        uOld = odtl->uvel;

    if (odtP->Iparticles) part->setFracC(); // for tracers only

    //-------------  advance the ODEs one step
#ifndef IMPLICIT // explicit solving
    for (int k = 0; k < neqns; k++) {
        for (int i = 0; i < odtl->ngrd; i++) {
            (*(odtl->props[k]))[i] = (*(odtl->props[k]))[i] + dtStep * k1rhs[k][i];
        }
    }
#else // implicit solving
    //vector<vector<double> > A, B, C, rhsImp;
    //A      = vector<vector<double> >(neqns, vector<double>(odtl->ngrd));
    //B      = vector<vector<double> >(neqns, vector<double>(odtl->ngrd));
    //C      = vector<vector<double> >(neqns, vector<double>(odtl->ngrd));
    //rhsImp = vector<vector<double> >(neqns, vector<double>(odtl->ngrd));
    for (int k = 0; k < neqns; k++){
        A[k].resize(odtl->ngrd);
        B[k].resize(odtl->ngrd);
        C[k].resize(odtl->ngrd);
        rhsImp[k].resize(odtl->ngrd);
    }
    computeMatrix(A,B,C,rhsImp,dtStep);
    for (int k = 0; k < neqns; k++){
        for (int i = 0; i < odtl->ngrd; i++){
            k1rhs[k][i] = (*(odtl->props[k]))[i] + rhsImp[k][i] + dtStep * k1rhs[k][i];
        }
    
        //Thomas-algorithm
        /* Modify the coefficients. */
        double id = 0.0;
        C[k][0]     /= B[k][0]; /* Division by zero risk. */
        k1rhs[k][0] /= B[k][0]; /* Division by zero would imply a singular matrix. */
        for (int i = 1; i < odtl->ngrd-1; i++) {
            id = (B[k][i] - C[k][i - 1] * A[k][i]); /* Division by zero risk. */
            C[k][i]    /= id; /* Last value calculated is redundant. */
            k1rhs[k][i] = (k1rhs[k][i] - k1rhs[k][i - 1] * A[k][i]) / id;
        }
        id = (B[k][odtl->ngrd-1] - C[k][odtl->ngrd-2] * A[k][odtl->ngrd-1]);
        k1rhs[k][odtl->ngrd-1] = (k1rhs[k][odtl->ngrd-1] - k1rhs[k][odtl->ngrd-2] * A[k][odtl->ngrd-1]) / id;
        
        /* Now back substitute. */
        (*(odtl->props[k]))[odtl->ngrd - 1] = k1rhs[k][odtl->ngrd - 1];
        for (int i = odtl->ngrd - 2; i >= 0; i--){
            (*(odtl->props[k]))[i] = k1rhs[k][i] - C[k][i] * (*(odtl->props[k]))[i + 1];
        }
    }
    return;
#endif
        
    if (odtP->Lrxn) {

        //------------ enforce species sum to 1

        for (int i = 0; i < odtl->ngrd; i++) {
            double sum = 0.0;
            for (int k = 0; k < odtl->nspc; k++)
                if (k != iN2) sum += odtl->yspc[k][i];
            odtl->yspc[iN2][i] = 1.0 - sum;
        }
        odtl->enforceYsBounds(odtl0->yspc); // neg spec set to 0, renorm.

        //---------- pressure

        odtl->pres = odtl->pres + dtStep*dpdt;
    } // endif Lrxn

    //------------ particles

    if (odtP->Iparticles) {
        part->diffuseTracerParticle(Gvel); // nothing done if !Ltracer
        computeXZpartPos(dtStep);
        for (int i = 0; i < part->nPart; i++) {
            if (!part->pActive[i]) continue;
            if (!part->Ltracer) {
                part->yPos[i] = part->yPos[i] + dtStep * P1rhs[0][i];
                part->uvel[i] = part->uvel[i] + dtStep * P1rhs[1][i];
                part->vvel[i] = part->vvel[i] + dtStep * P1rhs[2][i];
                part->wvel[i] = part->wvel[i] + dtStep * P1rhs[3][i];

                for (int j = 0; j < part->radial_ngrd; j++) {
                    //------If particles are present, their enthalpy and temperature change
                    if (odtP->Lrxn) {
                        double partTotEnth = part->cellEnth[i][j]*part->cellMass[i][j];
                        double partTotEnth_new = partTotEnth + dtStep * part->cellEnthSource[i][j];
// cout << endl << "part->cellEnthSource = " << part->cellEnthSource[i][j] << endl;        

                        part->cellEnth[i][j] = partTotEnth_new / part->cellMass[i][j];
// cout << endl << "part->cellEnth = " << part->cellEnth[i][j] << endl;        

                        //part->cellEnth[i][j] = part->cellEnth[i][j] + dtStep * part->cellEnthSource[i][j];
                        part->cellTemp[i][j] = part->setPartTemp(i, j);
                    } 
                    //---------If there is particle reaction, mass and density of the particles change.
                    if (odtP->Lprxn) {
                        part->cellMass[i][j] = part->cellMass[i][j] + dtStep * part->cellMassSource[i][j];
                        part->cellMass[i][j] = part->cellMass[i][j] - dtStep * part->cellMoistureSource[i][j];
                        part->moisture_tracker[i][j] = part->moisture_tracker[i][j] + dtStep * part->cellMoistureSource[i][j];
                        part->cellDens[i][j] = part->cellMass[i][j] / part->cellVolume[i][j];
                    } //end if odtP->Lprxn
                    
                } //end particle radial grid points for loop
            } //end if (!part->Ltracer)
        } //end particle for loop
        part->checkBoundsSetInactive(Gvel, time);
    } // end if(odtP->Iparticles loop)

    //---------- Update rho and grid

    if (odtP->Lrxn || odtP->Lprxn || odtP->ItableLookup)
        updateRhoAndGrid(dtStep, uOld);

    ////////////////////////////////////////////////////////////////////////////////////////// 
    ////////////////////////////////////////////////////////////////////////////////////////// 

    if (odtP->LsecondOrder) {

        vector<vector<double> > k2rhs(neqns, vector<double>(odtl->ngrd, 0.0));
        vector<vector<double> > cellMassSource1;
        vector<vector<double> > cellEnthSource1;
        
        if (odtP->Iparticles) {
            cellMassSource1 = vector<vector<double> >(part->nPart, vector<double>(part->radial_ngrd, 0.0));
            cellEnthSource1 = vector<vector<double> >(part->nPart, vector<double>(part->radial_ngrd, 0.0));
        }
        vector<vector<double> > P2rhs;
        if (odtP->Iparticles)
            P2rhs = vector<vector<double> >(4, vector<double>(part->nPart, 0.0));

        double dpdt1 = dpdt;
        brxr.odtl = odtl0; // implicit: integrate again over dt but with fluxes at t_{n+1}

        if (odtP->Iparticles) {
            cellMassSource1 = part->cellMassSource;
            cellEnthSource1 = part->cellEnthSource;
        }

        rhsf(k2rhs, P2rhs, dtStep, time);

        brxr.odtl = odtl; // reset the pointer

        //------------ save uOld for spatial advancement

        if (odtP->Lspatial)
            uOld = odtl->uvel;

        //-------------  advance the ODEs Again (second order) (ngrid should not have changed)

        for (int k = 0; k < neqns; k++)
            for (int i = 0; i < odtl->ngrd; i++)
                (*(odtl->props[k]))[i] = (*(odtl0->props[k]))[i] + 0.5 * dtStep * (k1rhs[k][i] + k2rhs[k][i]);

        if (odtP->Lrxn) {

            //------------ enforce species sum to 1

            for (int i = 0; i < odtl->ngrd; i++) {
                double sum = 0.0;
                for (int k = 0; k < odtl->nspc; k++)
                    if (k != iN2) sum += odtl->yspc[k][i];
                odtl->yspc[iN2][i] = 1.0 - sum;
            }
            odtl->enforceYsBounds(odtl0->yspc); // neg spec set to 0, renorm.

            //---------- pressure

            odtl->pres = odtl0->pres + 0.5 * dtStep * (dpdt1 + dpdt);

        } // endif Lrxn

        //------------ particles

        if (odtP->Iparticles) {
            part->diffuseTracerParticle(Gvel);
            computeXZpartPos(dtStep);
            for (int i = 0; i < part->nPart; i++) {
                if (!part->pActive[i]) continue;
                if (!part->Ltracer) {
                    part->yPos[i] = part0->yPos[i] + 0.5 * dtStep * (P1rhs[0][i] + P2rhs[0][i]);
                    part->uvel[i] = part0->uvel[i] + 0.5 * dtStep * (P1rhs[1][i] + P2rhs[1][i]);
                    part->vvel[i] = part0->vvel[i] + 0.5 * dtStep * (P1rhs[2][i] + P2rhs[2][i]);
                    part->wvel[i] = part0->wvel[i] + 0.5 * dtStep * (P1rhs[3][i] + P2rhs[3][i]);

                    for (int j = 0; j < part->radial_ngrd; j++) {
                        
                        if (odtP->Lrxn) {
                            part->cellEnth[i][j] = part0->cellEnth[i][j] + 0.5 * dtStep * (cellEnthSource1[i][j] + part->cellEnthSource[i][j]);
                            part->cellTemp[i][j] = part->setPartTemp(i, j);
                        } 
                        if (odtP->Lprxn) {
                            part->cellMass[i][j] = part0->cellMass[i][j] + 0.5 * dtStep * (cellMassSource1[i][j] + part->cellMassSource[i][j]);
                            part->cellDens[i][j] = part->cellMass[i][j] / part->cellVolume[i][j];
                        }
                    }
                }
            }
            part->checkBoundsSetInactive(Gvel, time);
        }

        //---------- Update rho and grid

        if (odtP->Lprxn || odtP->Lrxn || odtP->ItableLookup)
            updateRhoAndGrid(dtStep, uOld);

    } // endif LsecondOrder

    ////////////////////////////////////////////////////////////////////////////////////////// 
    ////////////////////////////////////////////////////////////////////////////////////////// 

    ////////////////////// Do grid inflow

    gridInflow(dtStep);

    ////////////////////// delete temporary lines

    if (odtl0) delete odtl0;
    if (part0) delete part0;

}

///////////////////////////////////////////////////////////////////////////////

/**Defines the partial derivative in time of the integration variables.
 * Finite volume discretization (conservative), with central differences
 * used for all derivatives.  Second order on uniform grids.
 *
 *  @param krhs \output vector of rates for equations solved [neq][ngrd]
 *  @param Prhs \output vector of rates for particle equations [eq][ipart]
 *  @param dtStep \input single explicit step size as part of integration.
 */

void diffuser::rhsf(vector<vector<double> > &krhs, vector<vector<double> > &Prhs, double dtStep, double time) {

    setDiffusivities_HSP_RR();
    dpdt = setGasVelocity_return_dpdt();

    //----------Compute RHS for particle and compute gas sources before calling computeFluxes(),
    //setRhsTrn() and setRhsSrc() functions because these functions use gasSource vector.
    if (odtP->Iparticles) {
        for (int kk = 0; kk < part->iPtMass + 1; kk++)
            (*gSource)[kk] = vector<double>(odtl->ngrd, 0);
        part->computeRHSFAndSetGasSource(Gvel, Prhs, dtStep, lambda_f, dxML, gSource, time);

        if(odtP->Ltwoway && !part->Ltracer && !part->Lballistic) // two-way coupling of momentum exchange of paritcles and gas phase
            setRhsPartSrc(Prhs);                                 // there is no need for tracer and ballistic particles    
    }

    computeFluxes();

    setRhsTrn();

    setRhsSrc(dtStep);

    for (int k = 0; k < neqns; k++){ // loop over equations
        for (int i = 0, ip = 1; i < odtl->ngrd; i++, ip++){ // loop over grid points
#ifndef IMPLICIT
            if(odtP->Ltwoway && !part->Ltracer && !part->Lballistic) // two-way coupling of momentum exchange of paritcles and gas phase
                krhs[k][i] = rhsTrn[k][i] + rhsSrc[k][i] + rhsPartSrc[k][i];
            else    
                krhs[k][i] = rhsTrn[k][i] + rhsSrc[k][i]; // full implicit
#else
            if(odtP->Ltwoway && !part->Ltracer && !part->Lballistic) // two-way coupling of momentum exchange of paritcles and gas phase
                krhs[k][i] = rhsSrc[k][i] + rhsPartSrc[k][i];
            else    
            //krhs[k][i] = 0.5*rhsTrn[k][i] + rhsSrc[k][i]; // Crank-Nicolson
                krhs[k][i] = rhsSrc[k][i]; // full implicit
#endif
        }
    }
}


///////////////////////////////////////////////////////////////////////////////

/** Compute the transport part of the RHSF.
 *  Essentially the divergence of the property flux
 */

void diffuser::setRhsTrn() {

    int nlimit;
    if (odtP->Lstrang) nlimit = iptEnth + 1; //for LEM no velocities; iptEnth=0
    else nlimit = neqns;

    for (int k = 0; k < nlimit; k++) { // loop over equations
        for (int i = 0, ip = 1; i < odtl->ngrd; i++, ip++) { // loop over grid points

            rhsTrn[k][i] = -(flxProp[k][ip] - flxProp[k][i]) / (dxML[i] * odtl->rho[i]);

            // ------------------------------- Spatial Modification

            if (odtP->Lspatial)
                rhsTrn[k][i] /= odtl->uvel[i];

                // -------------------------------add extra sources for particle reaction

            else if (odtP->Iparticles && odtP->Lprxn) {
                rhsTrn[k][i] /= odtl->voidFrac[i];
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

/** Compute momentum source term exchanged between the particles and gas phase
 *  two-way coupling
 *  
 *  Guangyuan Sun
 */

void diffuser::setRhsPartSrc(vector<vector<double> > &Prhs) {

    rhsPartSrc = vector<vector<double> >(neqns, vector<double>(odtl->ngrd, 0.0));

//     rhsPartSrc[0].resize(odtl->ngrd, 0.0);
//     rhsPartSrc[1].resize(odtl->ngrd, 0.0);
//     rhsPartSrc[2].resize(odtl->ngrd, 0.0);

    for (int i = 0; i < part->nPart; i++) {

        double pMass = part->nInPseudoPart[i]*part->pDens0[i]*4./3.*3.14159*part->pRadi[i]*part->pRadi[i]*part->pRadi[i];

        part->iyPos[i] = odtl->linePositionToIndex(part->yPos[i], true);

        rhsPartSrc[0][part->iyPos[i]] -= Prhs[1][i] * pMass / (dxML[part->iyPos[i]]*odtl->rho[part->iyPos[i]]);
        rhsPartSrc[1][part->iyPos[i]] -= Prhs[2][i] * pMass / (dxML[part->iyPos[i]]*odtl->rho[part->iyPos[i]]);
        rhsPartSrc[2][part->iyPos[i]] -= Prhs[3][i] * pMass / (dxML[part->iyPos[i]]*odtl->rho[part->iyPos[i]]);

    }
}

///////////////////////////////////////////////////////////////////////////////

/** Compute the source term part of the RHSF 
 *  Spatial fix is at the end
 */

void diffuser::setRhsSrc(double dtStep) {

    //------------ init to zero since all are used in rhsf but not all are explicitly set here.

    rhsSrc = vector<vector<double> >(neqns, vector<double>(odtl->ngrd, 0.0));

    //------------ momentum sources

    if (!odtP->Llem) {
        for (int i = 0; i < odtl->ngrd; i++) {

            rhsSrc[iptUvel][i] -= odtP->dPdx / odtl->rho[i]; // intended for channel flows

            if (odtP->Lbuoyant)
                rhsSrc[iptUvel][i] += (odtl->rho[odtl->ngrd - 1] - odtl->rho[i]) * 9.81 / odtl->rho[i];
            // FALKO debug: hier Auftrieb in v-Richtung einbauen
        }
        if (odtP->Iparticles && odtP->Lprxn) { // add the mass source contribution from eta
            for (int i = 0; i < part->line->ngrd; i++) {
                rhsSrc[iptUvel][i] += (*gSource)[part->iPtUvel][i];
                rhsSrc[iptVvel][i] += (*gSource)[part->iPtVvel][i];
                rhsSrc[iptWvel][i] += (*gSource)[part->iPtWvel][i];
            }
        }
    }

    //------------ eta sources
    if (odtP->Ieta)
        odtl->etaTools->computeSourceTerms();
    

    //------------ mom sources'

    if (odtP->Imom)
        odtl->momTools->computeSourceTerms();

    //------------ Yi, h sources for combustion

    if (odtP->Lrxn && !odtP->Lstrang) {

        for (int i = 0, ip = 1; i < odtl->ngrd; i++, ip++) {

            if (odtP->LimplicitChem) { //------------ implicit rxn integration

                brxr.setFluxSources(i, ip, 1.0 / (dxML[i] * odtl->rho[i]), flxProp, iptYspc, iptEnth,
                        (dpdt + radSource_G[i]) / odtl->rho[i]);

                double stepSize = (!odtP->Lspatial) ? dtStep : dtStep / odtl->uvel[i];
                brxr.integrateCell(i, stepSize);

                rhsSrc[iptEnth][i] = (dpdt + radSource_G[i]) / odtl->rho[i];

                if (odtP->Iparticles) { //added by abinash
                    rhsSrc[iptEnth][i] += (*(gSource))[part->iPtEnth][i];
                }

                for (int k = 0; k < odtl->nspc; k++)
                    rhsSrc[iptYspc + k][i] = brxr.meanRates[k];

            } else { //------------- explicit reaction

                rhsSrc[iptEnth][i] = (dpdt + radSource_G[i]) / odtl->rho[i];

                if (odtP->Iparticles) { //added by abinash
                    rhsSrc[iptEnth][i] += (*(gSource))[part->iPtEnth][i];
                }

                for (int k = 0; k < odtl->nspc; k++)
                    rhsSrc[iptYspc + k][i] = rrSpc[k][i] / odtl->rho[i];

            } // end test implicit or explicit rxn
        } // end loop over cells
    } // end if rxn

    //------------ Fix for spatial
    if (odtP->LhasTemp && odtP->LheatedChannel && timeFS > 0.0){
        if (!(uBulk != uBulk)){ // test if uBulk is not NaN
            for (int i = 0; i < odtl->ngrd; i++){
                rhsSrc[iptTemp][i] += odtl->uvel[i] / uBulk;
                //rhsSrc[iptScalar][i] += 1.0;
            }
        }
    }
    else if (odtP->LhasTemp){
        for (int i = 0; i < odtl->ngrd; i++)
            rhsSrc[iptTemp][i] = 0.0;
    }
    
    if (odtP->Lspatial)
        for (int k = 0; k < neqns; k++) // loop over equations
            for (int i = 0, ip = 1; i < odtl->ngrd; i++, ip++) // loop over grid points
                rhsSrc[k][i] /= odtl->uvel[i];
     
    
}

/*
 * calculate: convecting velocity do to incoming jets and dilatation
 * output: new cell sizes dxML[i]
 */

void diffuser::opposedJets(double timeStep) { //input = time-t0


    // copy of gloabal variables
    int ngrd = odtl->ngrd;

    // create local constants    
    double domainlength = odtl->posf.at(ngrd) - odtl->posf.at(0);
    double pos, pos_1, posf;
    double temp = 0.0;

    vector<double> yi(odtl->nspc); // species at point
    vector<double> vf(ngrd + 1, 0.0); // velocity at face
    vector<double> u_D(ngrd, 0.0); // velocity at cell center from dilatation
    vector<double> u_R(ngrd, 0.0);
    vector<double> rhoOld = odtl->rho; // denstiy at cell center



    if (odtl->LhasRxn) {

        // update properties to new time step
        for (int i = 0; i < ngrd; i++) //why -1
        {
            odtl->getYspVecAtPt(i, yi);
            odtl->gas->setState_PY(odtl->pres, &yi[0]); //set P and Y, hold T constant, change density

            if (odtl->strm->T1 != odtl->strm->T0) //for cold flow H doesn't change --> T doesn't change. only causes noise if left in
                odtl->gas->setState_HP(odtl->enth.at(i), odtl->pres, 1.E-12); //set H and P, incerment T (newton-iteration)
            odtl->temp.at(i) = odtl->gas->temperature();
            odtl->rho.at(i) = odtl->gas->density();
        }

        // calculate dilatation velocity       
        for (int i = 0; i < ngrd; i++) {
            double dRho_dt = (odtl->rho.at(i) - rhoOld.at(i)) / timeStep; // change of density in cell center
            temp += (1 / 3.0) / odtl->rho.at(i) * dRho_dt;
            u_D.at(i) = temp;
        }

        for (int i = 0; i < ngrd; i++) //constraint of BC is modeled as linear distribution of dilatation. u_D_f[ngrd]=u_D[ngrd-1]; u_D_f[0]=0
            u_D.at(i) += (odtl->pos.at(i) - odtl->posf.at(0)) / (domainlength)*(u_D.at(ngrd - 1) - 0.0); //u_D.at(i)=0.0;
    }

    // calculate u_R
    for (int i = 0; i < ngrd; i++) {
        if (odtl->unifCompress) //model 1: uniform compression
            u_R.at(i) = (odtl->pos.at(i) - odtl->posf.at(0))*(odtP->uBChi - odtP->uBClo) / domainlength + odtP->uBClo + u_D.at(i);
        else //model 2: compression based on velocity
            u_R.at(i) = odtl->uvel.at(i) + u_D.at(i);
    }

    // calculate velocity at faces. Used for updating cell sizes
    vf.at(0) = odtP->uBClo;
    for (int i = 1; i < ngrd; i++) {
        posf = odtl->posf.at(i);
        pos_1 = odtl->pos.at(i - 1);
        pos = odtl->pos.at(i);
        vf.at(i) = u_R.at(i - 1) + (posf - pos_1) / (pos - pos_1) * (u_R.at(i) - u_R.at(i - 1));
    }
    vf.at(ngrd) = odtP->uBChi;

    // update cell sizes & calc velocity flowing of odtline
    for (int i = 0; i < ngrd; i++) {
        double deltaCellLength = (vf.at(i + 1) - vf.at(i)) * timeStep;
        dxML.at(i) += deltaCellLength;
        if (!odtP->Llem) odtl->vvel[i] = -deltaCellLength / timeStep; //[m/s]
    }

    // create new cells at left and right boudary: 
    dxML.insert(dxML.begin(), abs(odtP->uBClo) * timeStep); // dxML[0]= |uBClo| * timeStep
    dxML.push_back(abs(odtP->uBChi) * timeStep); // dxML[ngrd+1]= |uBChi| * timeStep

    odtl->ngrd += 2;
    ngrd = odtl->ngrd;
    odtl->ngrdf = ngrd + 1;

    // CHECK: length of domain with new dxML
    double domainLengthPrime = 0.0;
    for (int i = 0; i < ngrd; i++)
        domainLengthPrime += dxML.at(i);

    if (domainLengthPrime > domainlength * (1.0 + 1.0e-5) || domainLengthPrime < domainlength * (1.0 - 1.0e-5)) {
        cout << endl << scientific << "Error in In diffuser::opposedJets. New domain length = " << domainLengthPrime << ", original domain length = " << domainlength << endl;
        exit(0);
    }

    //CHECK: for any negative cell sizes
    for (int i = 0; i < (int)dxML.size(); i++) {
        if (dxML.at(i) < 0.0) {
            ofstream myfile;
            myfile.open("OutputCellSizes.txt", ios::trunc); //ios::app
            myfile << "cell, size " << endl;
            for (int j = 0; j < (int)dxML.size(); j++) {
                myfile << scientific << j << "     " << dxML[j];
                if (dxML[j] <= 0.0) myfile << "<---" << endl; //create arrow at negative position 
                else myfile << endl;
            }
            myfile << "timeStep= " << timeStep << endl;
            myfile.close();
            cout << "Negative cell size created in diffuser::opposedJets. Cell sizes printed to OutputCellSizes.txt. End program. " << endl;
            exit(0);
        }
    }


    // resize pos & posf. Values will be set in function updateGrids
    odtl->pos.resize(ngrd, 0.0);
    odtl->posf.resize(ngrd + 1, 0.0);

    // create states for new cells: Yi, enthalpy, rho, phase, molec (diffusion coeff), u,v,w; fuel->...........<-air

    if (!odtP->Llem) {
        odtl->uvel.insert(odtl->uvel.begin(), odtP->uBClo);
        odtl->vvel.insert(odtl->vvel.begin(), 0.0);
        odtl->wvel.insert(odtl->wvel.begin(), 0.0);
        odtl->uvel.push_back(odtP->uBChi);
        odtl->vvel.push_back(0.0);
        odtl->wvel.push_back(0.0);
    }

    if (odtP->Lrxn) {
        for (int k = 0; k < odtl->nspc; k++) {
            odtl->yspc.at(k).insert(odtl->yspc.at(k).begin(), odtl->strm->y1.at(k));
            odtl->yspc.at(k).push_back(odtl->strm->y0.at(k));
        }

        odtl->gas->setState_TPY(odtl->strm->T1, odtl->pres, &odtl->strm->y1[0]); //sets state, does not updat it 
        odtl->enth.insert(odtl->enth.begin(), odtl->gas->enthalpy_mass());
        odtl->rho.insert(odtl->rho.begin(), odtl->gas->density());
        odtl->phase.insert(odtl->phase.begin(), 0.0);
        odtl->molec.insert(odtl->molec.begin(), odtl->molec[1]); //diffusion coeff= air


        odtl->gas->setState_TPY(odtl->strm->T0, odtl->pres, &odtl->strm->y0[0]);
        odtl->enth.push_back(odtl->gas->enthalpy_mass());
        odtl->rho.push_back(odtl->gas->density());
        odtl->phase.push_back(0.0);
        odtl->molec.push_back(odtl->molec[1]); //diffusion coeff= air
    }

    // clear vectors created
    yi.clear();
    vf.clear();
    u_D.clear();
    u_R.clear();
    rhoOld.clear();

}

///////////////////////////////////////////////////////////////////////////////

/**Compute fluxes at cell faces for uvw velocity, species, enthalpy.
 * Make sure auxiliary variables are up to date.\n
 * Fluxes:                                                         \cond
 * uvw velocities: tau = -mu*grad(u)                             \n
 * species: j_k = -rho*D_k*grad(Y_k) - rho*D_k*Y_k/M_k * grad(M) \n
 * enthalpy: q = -lamda*grad(T) + sum_k( h_k*j_k )               \n\endcond
 *
 */

void diffuser::computeFluxes() {

    int i, im;

    double rho_f; // dummy face values
    double M_f;
    double Y_f;
    double h_f;
    double hjsum;
    
    //---------------- set the grid factor dd used below 

    vector<double> dd(odtl->ngrdf, 0.0); // set the grid factor
    for (int i = 1, im = 0; i < odtl->ngrd; i++, im++) // interior
        dd[i] = 2.0 / (dxML[im] + dxML[i]);
    dd[0] = 2.0 / dxML[0];
    dd[odtl->ngrd] = 2.0 / dxML[odtl->ngrd - 1];
    if (odtP->bcType == 1) { // periodic
        dd[0] = 2.0 / (dxML[0] + dxML[odtl->ngrd - 1]);
        dd[odtl->ngrd] = dd[0];
    }

    //----------------- mmw for reaction cases

    if (odtP->Lrxn)
        mmw = odtl->getMMW();

    //---------- eta

    if ((odtP->IetaType == odtP->IETATYPE_TABLELOOKUP || odtP->IetaType == odtP->IETATYPE_DEFAULT)
            && odtP->Ieta) {
        odtl->etaTools->computeFluxes(dd);
    }

    //---------- mom

    if (odtP->Imom)
        odtl->momTools->computeFluxes(dd);
    /////////////////// Do interior faces ////////////////////////

    for (i = 1, im = 0; i < odtl->ngrd; i++, im++) {
        
        //---------- uvw velocities
        if (!odtP->Llem) { //LEM has no velocities
            if (!odtP->LmultiPhase){
                // all the same phase
                flxProp[iptUvel][i] = -dd[i] * visc_f[i]*(odtl->uvel[i] - odtl->uvel[im]);
                flxProp[iptVvel][i] = -dd[i] * visc_f[i]*(odtl->vvel[i] - odtl->vvel[im]);
                flxProp[iptWvel][i] = -dd[i] * visc_f[i]*(odtl->wvel[i] - odtl->wvel[im]);
                if (odtP->LhasTemp){
                    flxProp[iptTemp][i] = -dd[i] * lambda_f[i]*(odtl->temp[i] - odtl->temp[im]);
                }
            }
            else{ // multi-phase
                if (odtl->phase[i] == odtl->phase[im]){
                    // same phase
                    flxProp[iptUvel][i] = -dd[i] * visc_f[i]*(odtl->uvel[i] - odtl->uvel[im]);
                    flxProp[iptVvel][i] = -dd[i] * visc_f[i]*(odtl->vvel[i] - odtl->vvel[im]);
                    flxProp[iptWvel][i] = -dd[i] * visc_f[i]*(odtl->wvel[i] - odtl->wvel[im]);
                    if (odtP->LhasTemp){
                        flxProp[iptTemp][i] = -dd[i] * lambda_f[i]*(odtl->temp[i] - odtl->temp[im]);
                    }
                }
                else{ // different phases
                    // first:   the velocity on either sides of the surface has to be
                    //          the same
                    // second:  the flux nu1 * d(rho1 * u1) = nu2 * d(rho2 * u2) has to
                    //          be se same on either sides
                    //
                    //  uf is the velocity at the surface. The denseties are 
                    //  considered to be constant in both phases near the face.
                    //  uf = (u1 + A * B * u2) / (1 + A * B)
                    //  with
                    //  A = eta2 / eta1  and  B = dx1 / dx2  and  eta = nu * rho
                    //  where dx1 is the cell width next to the face of phase 1
                    //  and dx2 the cell width next to the face of phase 2.
                    //  B is nearly 1 due to the 2.5-rule for cell sizes.
                    double eta1 = odtl->molec[im];
                    double eta2 = odtl->molec[i];
                    // add part for Temp
                    double A   = 0.0;
                    double B   = 0.0;
                    double uf  = 0.0;
                    double vf  = 0.0;
                    double wf  = 0.0;
                    
                    // velocities
                    if (eta1 > eta2){
                        A  = eta2 / eta1;  // A < 1
                        B  = dxML[im] / dxML[i]; // B ~ 1
                        uf = (odtl->uvel[im] + A * B * odtl->uvel[i]) / (1 + A * B);
                        vf = (odtl->vvel[im] + A * B * odtl->vvel[i]) / (1 + A * B);
                        wf = (odtl->wvel[im] + A * B * odtl->wvel[i]) / (1 + A * B);
                        flxProp[iptUvel][i] = -2/dxML[im] * eta1 * (uf - odtl->uvel[im]);
                        flxProp[iptVvel][i] = -2/dxML[im] * eta1 * (vf - odtl->vvel[im]);
                        flxProp[iptWvel][i] = -2/dxML[im] * eta1 * (wf - odtl->wvel[im]);
                    }
                    else{
                        A  = eta1 / eta2;  // A <= 1
                        B  = dxML[i] / dxML[im]; // B ~ 1
                        uf = (odtl->uvel[i] + A * B * odtl->uvel[im]) / (1 + A * B);
                        vf = (odtl->vvel[i] + A * B * odtl->vvel[im]) / (1 + A * B);
                        wf = (odtl->wvel[i] + A * B * odtl->wvel[im]) / (1 + A * B);
                        flxProp[iptUvel][i] = -2/dxML[i] * eta2 * (odtl->uvel[i] - uf);
                        flxProp[iptVvel][i] = -2/dxML[i] * eta2 * (odtl->vvel[i] - vf);
                        flxProp[iptWvel][i] = -2/dxML[i] * eta2 * (odtl->wvel[i] - wf);
                    }
                    
                    // temperature
                    if (odtP->LhasTemp){
                        eta1 = odtl->lambda[im];
                        eta2 = odtl->lambda[i];
                        if (eta1 > eta2){
                            A  = eta2 / eta1;  // A < 1
                            B  = dxML[im] / dxML[i]; // B ~ 1
                            uf = (odtl->temp[im] + A * B * odtl->temp[i]) / (1 + A * B);
                            flxProp[iptTemp][i] = -2/dxML[im] * eta1 * (uf - odtl->temp[im]);
                        }
                        else{
                            A  = eta1 / eta2;  // A <= 1
                            B  = dxML[i] / dxML[im]; // B ~ 1
                            uf = (odtl->temp[i] + A * B * odtl->temp[im]) / (1 + A * B);
                            flxProp[iptTemp][i] = -2/dxML[i] * eta2 * (odtl->temp[i] - uf);
                        }
                    }
                }
            }
        }
        //--------- reaction

        if (odtP->Lrxn) {

            rho_f = linearInterpToFace(i, *odtl, odtl->rho);
            M_f = linearInterpToFace(i, *odtl, mmw);
            hjsum = 0.0;

            for (int k = 0; k < odtl->nspc; k++) {
                Y_f = linearInterpToFace(i, *odtl, odtl->yspc[k]);
                flxProp[iptYspc + k][i] = -dd[i] * rho_f * DmixYs_f[k][i]*((odtl->yspc[k][i] - odtl->yspc[k][im]) +
                        Y_f / M_f * (mmw[i] - mmw[im]));
                h_f = linearInterpToFace(i, *odtl, hsp[k]);
                hjsum += flxProp[iptYspc + k][i] * h_f;
            }
            flxProp[iptEnth][i] = -dd[i] * lambda_f[i]*(odtl->temp[i] - odtl->temp[im]) + hjsum;
        }
    }

    //////////////////// boundary faces ////////////////////

    
    if(odtP->Lsubdomain){
      
      for (int k = 0; k < neqns; k++) {
            flxProp[k][0] = 0.0;
            flxProp[k][odtl->ngrd] = 0.0;
      }
      
	
      if(subLeft_bc==1 && subRight_bc==0 ){

	
        if (odtP->Lrxn) {
            //Dirichlet for temperature
            flxProp[iptEnth][0] = -dd[0] * lambda_f[0] * (odtl->temp[0] - odtP->tempBClo);
        }
        
      }
      
	
      else if(subLeft_bc==0 && subRight_bc==1 ){
	  if (odtP->Lrxn) {
            //Dirichlet for temperature
            flxProp[iptEnth][odtl->ngrd] = -dd[odtl->ngrd] * lambda_f[odtl->ngrd] * (odtP->tempBChi - odtl->temp[odtl->ngrd - 1]);
	  }
       }
     return;   
     }
      
      
      
    
    


    //////////////////// periodic

    if (odtP->bcType == 1 || odtP->odtPprobType == 6) {

        int i = 0;
        int im = odtl->ngrd - 1;

        //---------- uvw velocities

        if (!odtP->Llem) {
            flxProp[iptUvel][0] = -dd[0] * visc_f[0]*(odtl->uvel[i]-(odtl->uvel[im] - odtP->pJump[0]));
            flxProp[iptVvel][0] = -dd[0] * visc_f[0]*(odtl->vvel[i]-(odtl->vvel[im] - odtP->pJump[1]));
            flxProp[iptWvel][0] = -dd[0] * visc_f[0]*(odtl->wvel[i]-(odtl->wvel[im] - odtP->pJump[2]));
            flxProp[iptUvel][odtl->ngrd] = flxProp[iptUvel][0];
            flxProp[iptVvel][odtl->ngrd] = flxProp[iptVvel][0];
            flxProp[iptWvel][odtl->ngrd] = flxProp[iptWvel][0];
            if (odtP->LhasTemp){
                cout << "\n\nERROR:\nPeriodic boundaries are not tested for passive/aktive scalar!\n\n";
                exit(0);
                flxProp[iptTemp][0] = -dd[0] * lambda_f[0]*(odtl->temp[i]-(odtl->temp[im] - odtP->pJump[3]));
                flxProp[iptTemp][odtl->ngrd] = flxProp[iptTemp][0];
            }
        }

        if (odtP->Lrxn && odtP->odtPprobType != 6) { //for premixed odt scalars are set below (outflow condition)

            //---------- species fluxes

            rho_f = linearInterpToFace(0, *odtl, odtl->rho);
            M_f = linearInterpToFace(0, *odtl, mmw);

            for (int k = 0; k < odtl->nspc; k++) {
                Y_f = linearInterpToFace(i, *odtl, odtl->yspc[k]);
                flxProp[iptYspc + k][0] = -dd[0] * rho_f * DmixYs_f[k][0]*((odtl->yspc[k][i] - odtl->yspc[k][im]) +
                        Y_f / M_f * (mmw[i] - mmw[im]));
                flxProp[iptYspc + k][odtl->ngrd] = flxProp[iptYspc + k][0];
            }

            //---------- heat flux

            hjsum = 0.0;
            for (int k = 0; k < odtl->nspc; k++) {
                h_f = linearInterpToFace(0, *odtl, hsp[k]);
                hjsum += flxProp[iptYspc + k][0] * h_f;
            }

            flxProp[iptEnth][0] = -dd[0] * lambda_f[0]*(odtl->temp[i] - odtl->temp[im]) + hjsum;
            flxProp[iptEnth][odtl->ngrd] = flxProp[iptEnth][0];
        }

    }
        ////////////////////// wall

    else if (odtP->bcType == 2) {

        //---------- uvw velocities : wall bc value

        for (int k = 0; k < neqns; k++) {
            flxProp[k][0] = 0.0;
            flxProp[k][odtl->ngrd] = 0.0;
        }

        if (!odtP->Llem) {
            flxProp[iptUvel][0] = -dd[0] * visc_f[0]* (odtl->uvel[0] - odtP->uBClo);
            flxProp[iptVvel][0] = -dd[0] * visc_f[0]* (odtl->vvel[0] - odtP->vBClo);
            flxProp[iptWvel][0] = -dd[0] * visc_f[0]* (odtl->wvel[0] - odtP->wBClo);

            flxProp[iptUvel][odtl->ngrd] = -dd[odtl->ngrd] * visc_f[odtl->ngrd]*(odtP->uBChi - odtl->uvel[odtl->ngrd - 1]);
            flxProp[iptVvel][odtl->ngrd] = -dd[odtl->ngrd] * visc_f[odtl->ngrd]*(odtP->vBChi - odtl->vvel[odtl->ngrd - 1]);
            flxProp[iptWvel][odtl->ngrd] = -dd[odtl->ngrd] * visc_f[odtl->ngrd]*(odtP->wBChi - odtl->wvel[odtl->ngrd - 1]);
        }
        //HEWSON
        if (odtP->Lrxn) {
            //no flux for species
            for (int k = 0; k < odtl->nspc; k++) {
                flxProp[iptYspc + k][0] = 0.0;
                flxProp[iptYspc + k][odtl->ngrd] = 0.0;
            }

            //Dirichlet for temperature
            flxProp[iptEnth][0] = -dd[0] * lambda_f[0] * (odtl->temp[0] - odtP->tempBClo);
            flxProp[iptEnth][odtl->ngrd] = -dd[odtl->ngrd] * lambda_f[odtl->ngrd] * (odtP->tempBChi - odtl->temp[odtl->ngrd - 1]);
        }

    }
        ////////////////// outflow

    else if (odtP->bcType == 3) {
        int k_min = 0;
        if (odtP->odtPprobType == 6) k_min = 3; //for premixed odt velocity is jump periodic, set above

        for (int k = k_min; k < neqns; k++) {
            flxProp[k][0] = 0.0;
            flxProp[k][odtl->ngrd] = 0.0;

        }

    }
        ////////////////////// wall on left, outflow right

    else if (odtP->bcType == 4 || odtP->bcType == 5) {

        /////// Wall BC on left

        for (int k = 0; k < neqns; k++)
            flxProp[k][0] = 0.0;

        //---------- uvw velocities : wall bc value
        if (!odtP->Llem) {
            flxProp[iptUvel][0] = -dd[0] * visc_f[0]*(odtl->uvel[0] - odtP->uBClo);
            flxProp[iptVvel][0] = 0.0; // -dd[0] * visc_f[0]*( odtl->vvel[0]-odtP->vBClo );
            flxProp[iptWvel][0] = 0.0; // -dd[0] * visc_f[0]*( odtl->wvel[0]-odtP->wBClo );
        }

        if (odtP->Lrxn) {

            //---------- species fluxes : zero flux condition

            rho_f = linearInterpToFace(0, *odtl, odtl->rho);

            odtl->gas->setState_PY(odtl->pres, &odtP->inletYsp[0]);
            //double inletMmw = odtl->gas->meanMolecularWeight();

            for (int k = 0; k < odtl->nspc; k++) {
                // flxProp[iptYspc+k][0] = -dd[0]*rho_f*DmixYs_f[k][0]*( (odtl->yspc[k][0] - odtP->inletYsp[k]) +
                //  odtP->inletYsp[k]/inletMmw*(mmw[0]-inletMmw) );
                flxProp[iptYspc + k][0] = 0.0;
            }

            if (odtP->bcType == 4) {

                //---------- heat flux 

                // flxProp[iptEnth][0] = 0.0; eimdb check with just heat flux

                hjsum = 0.0;
                for (int k = 0; k < odtl->nspc; k++) {
                    h_f = linearInterpToFace(0, *odtl, hsp[k]);
                    hjsum += flxProp[iptYspc + k][0] * h_f; //This is all 0 with no mass flux
                }

                flxProp[iptEnth][0] = -dd[0] * lambda_f[0]*(odtl->temp[0] - odtP->inletTemp) + hjsum;
            } else {

                //---------- heat flux 

                //flxProp[iptEnth][0] = 0.0; 

                hjsum = 0.0;
                for (int k = 0; k < odtl->nspc; k++) {
                    h_f = linearInterpToFace(0, *odtl, hsp[k]);
                    hjsum += flxProp[iptYspc + k][0] * h_f; //This is all 0 with no mass flux
                }

                flxProp[iptEnth][0] = -dd[0] * lambda_f[0]*(odtl->temp[0] - odtP->inletTemp) + hjsum;
            }
        }

        /////// Outflow BC on right

        for (int k = 0; k < neqns; k++)
            flxProp[k][odtl->ngrd] = 0.0;

    }
        ///////////////////// inflow on left and right side

    else if (odtP->bcType == 6) {
        // no flux at inlets
        for (int k = 0; k < neqns; k++) {
            flxProp[k][0] = 0.0;
            flxProp[k][odtl->ngrd] = 0.0;
        }
    }
    else if (odtP->bcType == 70) {
        // Falko's test boundary condition
        flxProp[iptUvel][0] = -dd[0] * visc_f[0]* (odtl->uvel[0] - odtP->uBClo);
        flxProp[iptVvel][0] = -dd[0] * visc_f[0]* (odtl->vvel[0] - odtP->vBClo);
        flxProp[iptWvel][0] = -dd[0] * visc_f[0]* (odtl->wvel[0] - odtP->wBClo);

        flxProp[iptUvel][odtl->ngrd] = -dd[odtl->ngrd] * visc_f[odtl->ngrd]*(odtP->uBChi - odtl->uvel[odtl->ngrd - 1]);
        flxProp[iptVvel][odtl->ngrd] = -dd[odtl->ngrd] * visc_f[odtl->ngrd]*(odtP->vBChi - odtl->vvel[odtl->ngrd - 1]);
        flxProp[iptWvel][odtl->ngrd] = -dd[odtl->ngrd] * visc_f[odtl->ngrd]*(odtP->wBChi - odtl->wvel[odtl->ngrd - 1]);

        if(odtP->LhasTemp){
            // wall BC
            flxProp[iptTemp][0] = -dd[0] * lambda_f[0]* (odtl->temp[0] - 0.0);
            flxProp[iptTemp][odtl->ngrd] = -dd[odtl->ngrd] * lambda_f[odtl->ngrd]*(0.0 - odtl->temp[odtl->ngrd - 1]);
            // slip BC
            //flxProp[iptScalar][0] = 1.0;
            //flxProp[iptScalar][odtl->ngrd] = 1.0;
        }
    } else if (odtP->bcType == -1) {
        // Proxy for individual calculation of BCs for each property depending on an input file.
        // currently not programmed
        *proc.ostrm << endl << "ERROR:\nCurrently the individual BCs are not programmed.";
        exit(0);
    } else {

        *proc.ostrm << endl << "****** ERROR IN DIFFUSER, BCTYPE = "
                << odtP->bcType << " unknown" << endl;
        exit(0);

    }

}

///////////////////////////////////////////////////////////////////////////////

/** Set cell positions in odtl
 */
void diffuser::updateGrids() {

    if (!odtP->Lrxn && !odtP->ItableLookup) {
        *proc.ostrm << "\nERROR diffuser::updateGrids only for combustion" << endl;
        exit(0);
    }

    //----------

    int i, im, ip;

    double sum = 0.0;
    for (int i = 0; i < odtl->ngrd; i++)
        sum += dxML[i];

    if (odtP->Lsubdomain)
      odtl->posf[0] = odtl->posf[0];
    else if (odtP->bcType == 3) {
        odtl->posf[0] = -0.5 * (sum - odtP->domainLength);
    } else
        odtl->posf[0] = 0.0;

    for (i = 1, im = 0; i < odtl->ngrdf; i++, im++)
        odtl->posf[i] = odtl->posf[im] + dxML[im];

    for (i = 0, ip = 1; i < odtl->ngrd; i++, ip++)
        odtl->pos[i] = 0.5 * (odtl->posf[i] + odtl->posf[ip]);

    setGridSize(); // redundant

    //---------- reset the domain size dolcheck iterate for const Ldomain

    odtl->Ldomain = (odtl->posf[odtl->ngrd] - odtl->posf[0]);

    return;
}

///////////////////////////////////////////////////////////////////////////////

/** Note, this does not include periodic jumps (dolcheck) 
 *
 * @param iface \input face index to interpolate to.
 * @param anyl \input line to interpolate.
 * @param vec \input variable being interpolated.
 * @return return interpolated variable at desired face.
 */
double diffuser::linearInterpToFace(int iface, const anyline &anyl,
        const vector<double> &vec) {

    double x1, x2, y1, y2;

    if (iface == 0) {
        if (odtP->Lperiodic) {
            x1 = anyl.pos[anyl.ngrd - 1] - anyl.Ldomain;
            x2 = anyl.pos[0];
            y1 = vec[anyl.ngrd - 1];
            y2 = vec[0];
        } else {
            return vec[0];
        }
    } else if (iface == anyl.ngrd) {
        if (odtP->Lperiodic) {
            x1 = anyl.pos[anyl.ngrd - 1];
            x2 = anyl.pos[0] + anyl.Ldomain;
            y1 = vec[anyl.ngrd - 1];
            y2 = vec[0];
        } else {
            return vec[anyl.ngrd - 1];
        }
    } else {
        x1 = anyl.pos[iface - 1];
        x2 = anyl.pos[iface];
        y1 = vec[iface - 1];
        y2 = vec[iface];
    }

    return y1 + (y2 - y1) / (x2 - x1) * (anyl.posf[iface] - x1);


}

///////////////////////////////////////////////////////////////////////////////

/**Outflow domains expand due to combustion.  Chop the domains
 * symmetrically.  That is, find the middle, and go Ld/2 on each side.
 * Spatial domains often contract due to flux conservation.
 * Issue a warning for any contracting temporal domains.
 * Ldomain is the current (expanded) size of the lines
 * domainLength is the desired size of the lines, so chop. 
 * NOTE, THIS MAY MERGE CELLS, ROUTINE INTENDED TO BE CALLED AFTER
 * DIFFUSION IS DONE, OTHERWISE, THE MEMBERS OF THIS CLASS (LIKE mmw)
 * WILL HAVE THE WRONG SIZE
 */

void diffuser::chopOutflowGrid() {

    if (odtl->probType == 4) //if opposedJets. do nothing
        return;

    if (!(odtP->bcType == 3 || odtP->bcType == 4 || odtP->bcType == 5) || (!odtP->Lrxn && !odtP->ItableLookup)) {
    // if (!(odtP->bcType == 3 || odtP->bcType == 4 || odtP->bcType == 5)) {
        *proc.ostrm << "\nERROR bcType != 3, 4, or 5, or not combustion in diffuser::chopOutflowGrid" << endl;
        exit(0);
    }

    if (odtl->posf[odtl->ngrd] - odtl->posf[0] == odtP->domainLength) {
        odtl->posf[0] = 0.0;
        odtl->posf[odtl->ngrd] = odtP->domainLength;
        odtl->pos[0] = 0.5 * odtl->posf[1];
        odtl->pos[odtl->ngrd - 1] = 0.5 * (odtl->posf[odtl->ngrd - 1] + odtl->posf[odtl->ngrd]);
        return;
    }

    //---------- Domain Contracted:
    // Domain:
    //                |----------------------------|             :
    // Extend the boundaries:
    // |~~~~~~~~~~~ * ----------------------------- * ~~~~~~~~~~~|
    // Then center the edge grid points
    // |~~~~~~*~~~~~-------------------------------~~~~~~~~*~~~~~|
    //
    // Note: rather than extend the cell, we can draw in BC streams

    if (odtl->posf[odtl->ngrd] - odtl->posf[0] < odtP->domainLength) {

        if (!odtP->Lspatial)
            *proc.ostrm << "\nWARNING: grid contracted in temporal outflow simulation "
                << (odtP->domainLength - odtl->posf[odtl->ngrd] + odtl->posf[0]) / odtP->domainLength * 100.0
                << "%" << endl;

        odtl->posf[0] = 0.0;
        odtl->posf[odtl->ngrd] = odtP->domainLength;

        odtl->pos[0] = 0.5 * odtl->posf[1];
        odtl->pos[odtl->ngrd - 1] = 0.5 * (odtl->posf[odtl->ngrd - 1] + odtl->posf[odtl->ngrd]);
    }
        //---------- Domain Expanded: two-sided or one-sided outflow
        //---------- two-sided: chop both sides; one-sided: chop right side

    else {

        //---------- odtl: split cell (if xhi/xlo not on face) 

        double xmid, xhi, xlo;
        int ip1, ip2;

        if (odtP->bcType == 3) { // outflow on both sides
            xmid = 0.5 * (odtl->posf[odtl->ngrd] + odtl->posf[0]);
            xhi = xmid + odtP->domainLength / 2.0;
        } else { // right side outflow
            xhi = odtP->domainLength;
            ip1 = 0;
        }

        //----------- chop right side

        ip2 = odtl->linePositionToIndex(xhi, false); // index of xhi
        if (odtl->posf[ip2 + 1] != xhi) {
            vector<double> splitPos(1, xhi - odtl->posf[ip2]);
            odtl->splitCell(ip2, 1, splitPos, false); // don't interpolate
            // cells above ip2 should be discarded
        }

        //----------- chop left side if two-sided outflow

        if (odtP->bcType == 3) {
            xlo = xhi - odtP->domainLength;
            ip1 = odtl->linePositionToIndex(xlo, true); // index of xlo
            if (odtl->posf[ip1] != xlo) {
                vector<double> splitPos(1, xlo - odtl->posf[ip1]);
                odtl->splitCell(ip1, 1, splitPos, false); // don't interpolate
                ip1++; // cells below ip1 should be discarded 
                ip2++; // splt lower cell, so inc upper cell indx
            }
        }

        //---------- odtl: delete the out of bounds by copying into itself

        *odtl = odtline(odtl, ip1, ip2, false);

        if (odtl->posf[0] != 0.0 || odtl->posf[odtl->ngrd] != odtP->domainLength) {
            odtl->posf[0] = 0.0;
            odtl->posf[odtl->ngrd] = odtP->domainLength;
            odtl->pos[0] = 0.5 * odtl->posf[1];
            odtl->pos[odtl->ngrd - 1] = 0.5 * (odtl->posf[odtl->ngrd - 1] + odtl->posf[odtl->ngrd]);
        }

        //---------- merge any small cells on the boundary due to split/castaway above

        if (odtl->posf[odtl->ngrd] - odtl->posf[odtl->ngrd - 1] < odtP->dxmin)
            odtl->merge2cells(odtl->ngrd - 2);
        if (odtP->bcType == 3)
            if (odtl->posf[1] - odtl->posf[0] < odtP->dxmin)
                odtl->merge2cells(0);

        if (odtP->ItableLookup) {
            odtl->etaTools->updateOdtLineVecs();
        }
        if (odtP->Iparticles) {
            odtl->setVoidFrac(this->part);
        }
        // ------------------------------------

    } // endif expanded

    odtl->Ldomain = odtl->posf[odtl->ngrd] - odtl->posf[0];

    if (odtP->Iparticles) part->checkBoundsSetInactive(Gvel);

    resetVarSizes();
    setGridSize();
    
}

///////////////////////////////////////////////////////////////////////////////

/** Interpolate from cell centers to cell faces.
 * Note, diffc starts as cell centers and ends as cell faces.
 * @param anyl \input line upon which to interpolate from.
 * @param diffc \input output vector of face values interpolated from cell centers.
 */

void diffuser::interpDiffusionCoeffToFaces(const anyline &anyl, vector<double> &diffc) {

    double dx1, dx2, k1, k2;
    int i, im;

    double dfirst; // store the first face diff till end
    double dlast; // store the last face diff till end
    double dmb1, dmb2; // dummies for overwriting of diffc

    if (odtP->Lperiodic) {
        i = 0;
        im = anyl.ngrd - 1;

        dx1 = anyl.posf[anyl.ngrd] - anyl.pos[im];
        dx2 = anyl.pos[0] - anyl.posf[0];
        k1 = diffc[im];
        k2 = diffc[0];
        dfirst = k1 * k2 * (dx1 + dx2) / (k1 * dx2 + k2 * dx1); // first face
        dlast = dfirst; // last face
    } else {
        dfirst = diffc[0]; // first face
        dlast = diffc[anyl.ngrd - 1]; // last face
    }

    dmb2 = diffc[0];
    for (i = 1, im = 0; i < anyl.ngrdf - 1; i++, im++) { // central faces
        dx1 = anyl.pos[im] - anyl.posf[im]; //fixdx
        dx2 = anyl.pos[i] - anyl.posf[i]; //fixdx
        dmb1 = dmb2; // e.g. face 3 dep on cells 2,3
        dmb2 = diffc[i]; // but cell 2 was changed on last call
        k1 = dmb1; // and before change 3, store it for next
        k2 = diffc[i]; // so use two dummies for storage
        if(abs(k1)<1e-9 && abs(k2)<1e-9)            // C. Sch. zero diffusivitites if non-premixed and stream consists of one species only
            diffc[i]=0.0;
        else
            diffc[i] = k1 * k2 * (dx1 + dx2) / (k1 * dx2 + k2 * dx1);
    }

    diffc[0] = dfirst; // insert the first and last
    diffc[anyl.ngrd] = dlast;

}

///////////////////////////////////////////////////////////////////////////////

/** Adapt during diffusion for spatial cases for which grid contraction
 *  results in small grid cells
 */

void diffuser::adaptGridsIfNeeded() {

    if (odtl->LhasRxn) odtl->setTempVec(); //may be redundant

    if (odtP->ItableLookup) odtl->etaTools->updateOdtLineVecs();
    //-------------------------

    if (odtP->Lspatial && *min_element(dxML.begin(), dxML.end()) <= odtP->dxmin) {

        *proc.ostrm << endl << "#------- ADAPTING DURING DIFFUSION";

        odtl->meshAdapter.adaptGrid(0, odtl->ngrd - 1);

        resetVarSizes();

        setGridSize();

        if (odtP->ItableLookup) {
            odtl->etaTools->updateOdtLineVecs();
        }
        if (odtP->Iparticles) {
            odtl->setVoidFrac(this->part);
        }

    }

    if (odtl->probType == 4) {
        *proc.ostrm << endl << "#------- ADAPTING DURING DIFFUSION";

        int girdPntsToAdapt = 5; //odtl->ngrd/2;

        odtl->meshAdapter.adaptGrid(0, 0 + girdPntsToAdapt);

        odtl->meshAdapter.adaptGrid(odtl->ngrd - 1 - girdPntsToAdapt, odtl->ngrd - 1);

        //	odtl->meshAdapter.adaptGrid(0,odtl->ngrd-1);

        resetVarSizes();

        setGridSize();

        if (odtP->ItableLookup) {
            odtl->etaTools->updateOdtLineVecs();
        }
        if (odtP->Iparticles) {
            odtl->setVoidFrac(this->part);
        }
    }



}

///////////////////////////////////////////////////////////////////////////////

/** Implement inlet boundary condition.  Initially implemented for porous ethylene
 *  wall configuration.
 *  @param dtStep \input current step size used to compute amount of inlet
 */

void diffuser::gridInflow(double dtStep) {

    if (odtP->bcType != 4) {
        return;
    }
    if (!odtP->Lspatial) {
        *proc.ostrm << endl << "ERROR diffuser::gridInflow for spatial only" << endl;
        exit(0);
    }
    double dxIn = 0.0;
    double uNew = 0.0;
    //---------- move the grid over to make room for inflow mass
    if (!odtP->Llem) {
        uNew = (-odtP->vBClo * dtStep + sqrt(odtP->vBClo * odtP->vBClo * dtStep * dtStep +
                (odtl->posf[1] - odtl->posf[0]) * odtl->uvel[0] * odtP->vBClo * dtStep)) /
                (odtl->posf[1] - odtl->posf[0]);

        uNew = odtl->uvel[0]; //doldb
        dxIn = dtStep * odtP->vBClo / uNew; // this "volume" mass in
    }

    for (int i = 0; i < odtl->ngrd; i++) {
        odtl->pos[i] += dxIn;
        odtl->posf[i] += dxIn;
    }
    odtl->posf[odtl->ngrd] += dxIn; // do the last face too

    //---------- split the first cell without interpolation

    vector<double> splitPos(1, odtl->pos[0] - odtl->posf[0]); // loc to splt: from left edge of cell
    odtl->splitCell(0, 1, splitPos, false);

    //---------- fix the first and second cell positions

    odtl->posf[0] = 0; // reset the 1st and 2nd cell positions
    odtl->posf[1] = dxIn;
    odtl->pos[0] = 0.5 * (odtl->posf[1] + odtl->posf[0]);
    odtl->pos[1] = 0.5 * (odtl->posf[2] + odtl->posf[1]);

    //---------- fill the contents of first cell with inflow properties
    if (!odtP->Llem) {
        odtl->uvel[0] = uNew;
        odtl->vvel[0] = 0.0; // line vel accounted for with shift
        odtl->wvel[0] = 0.0;
    }
    if (odtP->Lrxn) {
        for (int k = 0; k < odtl->nspc; k++)
            odtl->yspc[k][0] = odtP->inletYsp[k];
        odtl->gas->setState_TPY(odtP->inletTemp, odtl->pres, &odtP->inletYsp[0]);
        odtl->enth[0] = odtl->gas->enthalpy_mass();
        odtl->rho[0] = odtl->gas->density();
    }
    else if (odtP->ItableLookup) {
        odtl->eta[0][0] = 1.0;
        odtl->eta[1][0] = 1.8714e+06;
        odtl->rho[0] = odtl->etaTools->luTable->getValOf(1.0, 1.8714e+06, odtl->etaTools->luTable->index.rho_index);
    }

    //---------- fix up the grid 
    // If merge the small inlet with neighbor, a large cell could be created --> todo fix this

    if (dxIn < odtP->dxmin) { // merge the small inlet cell with neighbor
        odtl->merge2cells(0);
        dxML[0] += dxIn;
    } else { // account for new split cell
        dxML.insert(dxML.begin(), dxIn);
        resetVarSizes();
        setGridSize();
    }

    if (odtP->ItableLookup) {
        odtl->etaTools->updateOdtLineVecs();
    }


    //---------- reset the domain size dolcheck iterate for const Ldomain

    odtl->Ldomain = (odtl->posf[odtl->ngrd] - odtl->posf[0]);

}

///////////////////////////////////////////////////////////////////////////////

/** Set gas velocity at cell faces for particle drag law.
 * \cond
 * Sets Gvel member
 * 
 * div(v) = -1/(gamma*P) * dp/dt + curlyU
 * Integrate this over each cell:
 * delta_Gvel = ve-vw = -1/P*dp/dt * delta_x/gamma + integral(curlyU dx)
 * To evaluate:
 * First, find integral(curlyU dx) in each cell and store it in Gvel.
 * For open domains, this is delta_Gvel and dp/dt = 0
 * For closed domains, evaluate dp/dt by integrating the first eqn over the whole domain:
 * vE - vW = -1/P*dp/dt * sum(delta_x/gamma) + sum(integral(curlyU dx)), where sums are over 
 *    all cells, and vE, vW are currently zero (fixed domain boundaries).
 * Then dp/dt = [ sum(integral(curlyU dx)) ] / [ 1/P * sum(dx/gamma) ]
 * Now, given dp/dt, can evaluate delta_Gvel in each cell be appending the -1/P*dp/dt *delta_x/gamma term.
 *
 * Then, given delta_Gvel, stored in Gvel, we can compute Gvel on cell faces, relative to the 
 *     zero velocity expansion center.
 *  \endcond

 */

double diffuser::setGasVelocity_return_dpdt() {

    //--------- should we be in here?

    if (!odtP->Lrxn) { // no rxn then dpdt = 0 and Gvel = 0
        if (odtP->Iparticles)
            for (int i = 0; i < (int)Gvel.size(); i++)
                Gvel[i] = 0.0;
        return 0.0;
    } else { // reacting cases     
        if (!odtP->Iparticles && odtP->bcType != 1 && odtP->bcType != 2) // no particles in open domains
            return 0.0;
    }

    //--------- All below is for reacting cases that have closed domains or have particles or both:
    //--------- Particles,    open domain   --> compute Gvel, but dpdt = 0
    //--------- Particles,    closed domain --> compute Gvel, and dpdt
    //--------- No particles, closed domain --> compute first part of Gvel (holding int(curlyU dx)) to get dpdt
    //--------- In all these cases we do the first part:

    //--------- first part: compute int(curlyU dx) in each cell (used for particles and to compute dpdt alone)

    double dd, dsum1, dsum2;
    for (int i = 0, ip = 1; i < odtl->ngrd; i++, ip++) {

        dsum1 = 0.0;
        dsum2 = 0.0;
        for (int k = 0; k < odtl->nspc; k++) {
            dd = flxProp[iptYspc + k][ip] - flxProp[iptYspc + k][i] - rrSpc[k][i] * dxML[i];
            dsum1 += hsp[k][i] * dd;
            dsum2 += dd / odtl->gas->molecularWeight(k);
        }

        Gvel[i] = (1.0 - invGamma[i]) / odtl->pres * // Intermediate = int(curlyU dx)
                (-(flxProp[iptEnth][ip] - flxProp[iptEnth][i]) + dsum1) -
                mmw[i] / odtl->rho[i] * dsum2;
    }

    //-------- compute dpdt, open domain

    if (odtP->bcType != 1 && odtP->bcType != 2) // dpdt = 0 for no rxn or open domains
        dpdt = 0.0;

        //-------- compute dpdt, closed domain

    else {
        double intInvGamma = 0;
        dpdt = 0.0;
        for (int i = 0; i < odtl->ngrd; i++) {
            intInvGamma += invGamma[i] * dxML[i];
            dpdt += Gvel[i]; // intermediate
        }
        dpdt *= odtl->pres / intInvGamma;
    }

    if (odtP->Iparticles) {

        //-------- now update Gvel in each cell by adding the pressure term

        if (dpdt != 0.0)
            for (int i = 0; i < odtl->ngrd; i++)
                Gvel[i] += -dpdt / odtl->pres * dxML[i] * invGamma[i];

        //--------- Compute Gvel relative to the expansion center...

        //--------- find the cumulative velocity

        double vBCright = 0.0; // vel of right bc relative to left bc (zero for closed domains)
        for (int i = 0; i < odtl->ngrd; i++)
            vBCright += Gvel[i];

        //--------- Convert (Gvel holding deltaGvel) to (Gvel holding Gvel relative to left BC)

        Gvel[0] = 0.0;
        for (int i = 1, im = 0; i < odtl->ngrdf; i++, im++)
            Gvel[i] += Gvel[im];

        //--------- Now subtract vBCright/2 from all cells to get the gas vel relative to the expansion center

        if (odtP->bcType == 3) {
            double velShift = vBCright / 2.0;
            for (int i = 0; i < odtl->ngrdf; i++)
                Gvel[i] -= velShift;
        }
    }

    return dpdt;
}

///////////////////////////////////////////////////////////////////////////////

/**
 * Enforcing continuity equation.
 *
 * As the combustion happens, density increases or decreases. The method updates the new dxML
 * according to the continuity equation.
 *
 * @param dtStep \input time step
 * @param uOld   \input old u values
 */


void diffuser::updateRhoAndGrid(double& dtStep, vector<double>& uOld) {

    double sumDxML = 0.0;	//subdomain decomposition
    if (odtl->probType == 4)
        opposedJets(dtStep); //calculate velocities: advec. + dilatation. Up-date cell sizes vector dxML
    else {
        vector<double> voidFracOld(odtl->voidFrac.size(), 0.0);
        vector<double> rhoOld(odtl->rho.size(), 0.0);
        rhoOld = odtl->rho;
        voidFracOld = odtl->voidFrac;
        odtl->setVoidFrac(this->part);

        if (odtl->odtP->ItableLookup) {
            odtl->etaTools->updateOdtLineVecs(true);
        }

        for (int i = 0; i < odtl->ngrd; i++) {
            double rhoGammaDeltaX_np1 = 0;
            if (!odtl->odtP->ItableLookup) {
                vector<double> yi(odtl->nspc);
                odtl->getYspVecAtPt(i, yi);
                odtl->gas->setState_PY(odtl->pres, &yi[0]);
                odtl->gas->setState_HP(odtl->enth[i], odtl->pres, 1.E-10);
                odtl->temp[i] = odtl->gas->temperature();
                odtl->rho[i] = odtl->gas->density();
            }
            //----------------Continuity (Calculation of new dxML[i])
            if (odtP->Iparticles && odtP->Lprxn) {
                rhoGammaDeltaX_np1 = rhoOld[i] * voidFracOld[i] * dxML[i] +
                        dtStep * (*gSource)[part->iPtMass][i] * dxML[i] * voidFracOld[i];

                if (odtP->LsecondOrder) { //second order integration in time
                    double rhoGammaDeltaX;
                    double rate_n;
                    double rate_np1;
                    rate_n = (*gSource)[part->iPtMass][i] * voidFracOld[i] * dxML[i];
                    rhoGammaDeltaX = rhoOld[i] * voidFracOld[i] * dxML[i];
                    rate_np1 = (*gSource)[part->iPtMass][i + 1] * dxML[i + 1] * odtl->voidFrac[i + 1];
                    rhoGammaDeltaX_np1 = rhoGammaDeltaX + dtStep + 0.5 * (rate_np1 + rate_n);
               }
                double dxMLOld = dxML[i];
                dxML[i] = rhoGammaDeltaX_np1 / (odtl->rho[i] * odtl->voidFrac[i]);
                Gvel[i] = (dxML[i]-dxMLOld)/dtStep;
            } else {
                if (odtP->Lspatial)
                    dxML[i] = rhoOld[i] * uOld[i] * dxML[i] / odtl->rho[i] / odtl->uvel[i];
                else
                    dxML[i] = rhoOld[i] * dxML[i] / odtl->rho[i];
            }
            if(odtP->Lsubdomain) sumDxML+= dxML[i];
        }
    }
    
     // subdomain decomposition; dealing with expansion of subdomains
    if(odtP->Lsubdomain){
        double sumDxML_sub=0.0;
        for(int k = 0; k< indNextSub_start; k++)
            sumDxML_sub += dxML[k];

        if((sumDxML_sub  + odtl->posf[0] - odtl->posf[indNextSub_start]) > 1.e-30){
            dx_exp += sumDxML_sub - odtl->posf[indNextSub_start] + odtl->posf[0];
            posNextSub_start = sumDxML_sub + posNextSub_start - odtl->posf[indNextSub_start] + odtl->posf[0];
        }


        double oldLength = odtl->posf[odtl->ngrd] - odtl->posf[0];
        //shrink begin
        if(oldLength - sumDxML > 0) {
            double addDx = 0.5 * (oldLength - sumDxML);
            dxML[0] = dxML[0] + addDx;
            dxML[odtl->ngrd-1] = dxML[odtl->ngrd-1] + addDx;
        }
        //shrink end
    }    
    // subdomain decomposition end  
    
        
    updateGrids();
    if (odtP->Iparticles) {
        odtl->setVoidFrac(this->part);
    }
    
    // FALKO debug: hier das Update der Dichte reinschreiben; wie haengt rho mit T zusammen?
    
    
    
}

////////////////////////////////////////////////////////////////////////////

/**
 * Computes the radiate heat source for gas and for particles
 */
void diffuser::setRadHeatSources() {

    if (odtP->Iradiation > 0) {
        vector<vector<double> > xMoleSp =
                vector<vector<double> >(odtl->ngrd, vector<double>(odtl->nspc));

        if (odtP->Lrxn) {
            for (int i = 0; i < odtl->ngrd; i++) {
                vector<double> yi(odtl->nspc); // working array
                odtl->getYspVecAtPt(i, yi);
                odtl->gas->setState_PY(odtl->pres, &yi[0]);
                odtl->gas->getMoleFractions(&xMoleSp[i][0]);
            }
        }

        //----------- set gas absorption coefficient ka if lookup table, 
        //-----------else pass zero size and compute in radiation

        vector<double> ka;
        if (odtP->ItableLookup) {
            ka.resize(odtl->ngrd);
            for (int i = 0; i < odtl->ngrd; i++)
                ka[i] = odtl->etaTools->luTable->
                    getValAtGridPoint(i, odtl->etaTools->luTable->index.kabs_index);
        } else
            ka = vector<double>();

        //------------ compute radiation sources: gas and particles

        if (odtP->ImomType != odtP->IMOMTYPE_SOOT) {
            rad.getRadHeatSource(xMoleSp, odtl->temp, odtl->pres, odtl->posf, radSource_G, radSource_P, ka, part);
        }
        else { // call with soot volume fraction

            vector<double> fvSoot(odtl->ngrd);
            odtl->momTools->setVolumeFraction(fvSoot);
            rad.getRadHeatSource(xMoleSp, odtl->temp, odtl->pres, odtl->posf, radSource_G, radSource_P, ka, part, fvSoot);
        }


    } else {
        if (part)
            radSource_P = vector<double>(part->nPart, 0.0);
        radSource_G = vector<double>(odtl->ngrd, 0.0);
    }

}

///////////////////////////////////////////////////////////////////////

/** helper function for sort function in particles::updateEddyInfoArray
 *
 *  @Guangyuan Sun 07/2012
 */

bool compare(diffuser::eddyCombination a, diffuser::eddyCombination b) {
    return a.eddyDeadTime < b.eddyDeadTime;
}

///////////////////////////////////////////////////////////////////////
/**Store the effect of active eddies in the process of diffustion.                 
 * called in odtSolver::sampleEddy                                                
 *
 * @author Guangyuan Sun 10/2012
 */

void diffuser::updateActiveEddy(double leftEdge, double rightEdge, double eddySize, 
                            double eddyLife, double uEddy, double vEddy, double wEddy,
                            double eddyEndTime, double Pvel, vector<int> iPartInEddy,
                            vector<int> iPartInEddyIC, vector<double> xPartPos, vector<double> zPartPos
                            ) {

    int nActiveEddy = ActiveEddy.size()+1;
    ActiveEddy.resize(nActiveEddy);
    ActiveEddy[nActiveEddy-1].eddyLeftEdge = leftEdge;
    ActiveEddy[nActiveEddy-1].eddyRightEdge = rightEdge;
    ActiveEddy[nActiveEddy-1].eddyLength   = eddySize;
    ActiveEddy[nActiveEddy-1].eddyLife     = eddyLife;
    ActiveEddy[nActiveEddy-1].eddyDeadTime = eddyEndTime;
    ActiveEddy[nActiveEddy-1].eddyUvel     = uEddy;
    ActiveEddy[nActiveEddy-1].eddyVvel     = vEddy;
    ActiveEddy[nActiveEddy-1].eddyWvel     = wEddy;
    ActiveEddy[nActiveEddy-1].eddyPvel     = Pvel; 
    ActiveEddy[nActiveEddy-1].eddyiPartIC  = iPartInEddyIC;
    ActiveEddy[nActiveEddy-1].eddyiPart    = iPartInEddy;
    ActiveEddy[nActiveEddy-1].eddyxPartPos = xPartPos;
    ActiveEddy[nActiveEddy-1].eddyzPartPos = zPartPos;

    sort (ActiveEddy.begin(),ActiveEddy.end(),compare);
}

///////////////////////////////////////////////////////////////////////

/** check and update the array of interacting particles during type-C or -IC interaction at the certain time step
 *
 *  Guangyuan Sun 09/11
 */

void diffuser::updatePartIxn(double time) {

    if(part->PeddyType == 2 || part->PeddyType == 3) {

        if(ActiveEddy.empty()) return;    

        std::vector<eddyCombination>::iterator it; 
        for ( it=ActiveEddy.begin() ; it != ActiveEddy.end(); ++it) {       
            for (int i=0; i<part->nPart; i++) {
                if (!part->pActive[i]) continue;                                     // skip inactive particles (e.g. wall collisions/outflow)


                bool inList = false;
                for (int j=0; j<it->eddyiPart.size(); j++) {
                    if (part->PeddyType == 3)
                    {
                        if (i == it->eddyiPartIC[j])
                            goto nextParticle;
                    }

                    if (i==it->eddyiPart[j]) {
                        // check and update the status of current particles in the array
                        inList = true;
                        if ( part->yPos[i]>it->eddyRightEdge || part->yPos[i]<it->eddyLeftEdge || it->eddyxPartPos[j]<-0.5*it->eddyLength || it->eddyxPartPos[j]>0.5*it->eddyLength || it->eddyzPartPos[j]<-0.5*it->eddyLength || it->eddyzPartPos[j]>0.5*it->eddyLength ) {
                            it->eddyiPart.erase(it->eddyiPart.begin()+j);
                            it->eddyxPartPos.erase(it->eddyxPartPos.begin()+j);
                            it->eddyzPartPos.erase(it->eddyzPartPos.begin()+j);

                            // sun increment crossing frequency of particle during type-C interaction for Snyder&Lumley

//                             part->partCross[i]++;
//                             if (time >= 0.2637 && time < 0.3637) part->partCross1[i]++;
//                             if (time >= 0.3637 && time < 0.4637) part->partCross2[i]++;
//                             if (time >= 0.4637 && time < 0.5637) part->partCross3[i]++;
//                             if (time >= 0.5637 && time < 0.6637) part->partCross4[i]++;
//                             if (time >= 0.6637 && time <= 0.7637) part->partCross5[i]++;
//                             if (time >= 0.2637 && time <= 0.7637) part->partCrossFirstFiveSec[i]++;
//                             if (time >= 0 && time <= 0.7637) part->partCrossFirstFiveSecAll[i]++;

                            goto nextParticle;
                        }
                    }
                }
                if (!inList) {
                    if (part->yPos[i] >= it->eddyLeftEdge && part->yPos[i] <= it->eddyRightEdge) {
                        if( ((it->eddyLeftEdge == part->yPos[i])  && (part->vvel[i] < 0.0)) ||   // skip if part on edge and moving out of eddy box
                                ((it->eddyRightEdge == part->yPos[i]) && (part->vvel[i] > 0.0)) ) 
                            continue;
                        it->eddyiPart.push_back(i);
                        it->eddyxPartPos.push_back(0.);
                        it->eddyzPartPos.push_back(0.);

                        // sun increment interaction frequency of particle during type-C interaction for Snyder&Lumley
//                         part->partIxn[i]++;
//                         if (time >= 0.2637 && time < 0.3637) part->partIxn1[i]++;
//                         if (time >= 0.3637 && time < 0.4637) part->partIxn2[i]++;
//                         if (time >= 0.4637 && time < 0.5637) part->partIxn3[i]++;
//                         if (time >= 0.5637 && time < 0.6637) part->partIxn4[i]++;
//                         if (time >= 0.6637 && time <= 0.7637) part->partIxn5[i]++;
//                         if (time >= 0.2637 && time <= 0.7637) part->partIxnFirstFiveSec[i]++;
//                         if (time >= 0 && time <= 0.7637) part->partIxnFirstFiveSecAll[i]++;
                    }
                }
nextParticle:
                continue;
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////

/** compute x (eddyxPartPos) and z (eddyzPartPos) position of particles for every active eddy 
 *
 *  Guangyuan Sun 09/11
 */

void diffuser::computeXZpartPos(double dtStep) {
    if(part->PeddyType == 2 || part->PeddyType == 3) {
        std::vector<eddyCombination>::iterator it;
        for ( it=ActiveEddy.begin() ; it != ActiveEddy.end(); ++it) {
            for (int j=0; j<it->eddyiPart.size(); j++) {
                int iPart = it->eddyiPart[j];
                part->iyPos[iPart] = odtl->linePositionToIndex(part->yPos[iPart], true); //should be redundant: delete

                //            it->eddyxPartPos[j] += (odtl->uvel[part->iyPos[iPart]] + (part->TauP[iPart] / part->f[iPart] * part->AGx)) * dtStep
                //                + part->TauP[iPart] / part->f[iPart] * (part->uvel[iPart] - odtl->uvel[part->iyPos[iPart]] - part->TauP[iPart]/part->f[iPart]*part->AGx)
                //                * (1 - exp(-part->f[iPart]/part->TauP[iPart]*dtStep));
                //            it->eddyzPartPos[j] += (odtl->wvel[part->iyPos[iPart]] + (part->TauP[iPart] / part->f[iPart] * part->AGz)) * dtStep
                //                + part->TauP[iPart] / part->f[iPart] * (part->wvel[iPart] - odtl->wvel[part->iyPos[iPart]] - part->TauP[iPart]/part->f[iPart]*part->AGz)
                //                * (1 - exp(-part->f[iPart]/part->TauP[iPart]*dtStep));

                //            it->eddyxPartPos[j] += - odtl->uvel[part->iyPos[iPart]] * dtStep + (odtl->uvel[part->iyPos[iPart]] + (part->TauP[iPart] / part->f[iPart] * part->AGx)) * dtStep
                //                + part->TauP[iPart] / part->f[iPart] * (part->uvel[iPart] - odtl->uvel[part->iyPos[iPart]] - part->TauP[iPart]/part->f[iPart]*part->AGx)
                //                * (1 - exp(-part->f[iPart]/part->TauP[iPart]*dtStep));
                //            it->eddyzPartPos[j] += - odtl->wvel[part->iyPos[iPart]] * dtStep + (odtl->wvel[part->iyPos[iPart]] + (part->TauP[iPart] / part->f[iPart] * part->AGz)) * dtStep
                //                + part->TauP[iPart] / part->f[iPart] * (part->wvel[iPart] - odtl->wvel[part->iyPos[iPart]] - part->TauP[iPart]/part->f[iPart]*part->AGz)
                //                * (1 - exp(-part->f[iPart]/part->TauP[iPart]*dtStep));

                it->eddyxPartPos[j] += - it->eddyUvel * dtStep + (odtl->uvel[part->iyPos[iPart]] + (part->TauP[iPart] / part->f[iPart] * part->AGx)) * dtStep
                    + part->TauP[iPart] / part->f[iPart] * (part->uvel[iPart] - odtl->uvel[part->iyPos[iPart]] - part->TauP[iPart]/part->f[iPart]*part->AGx)
                    * (1 - exp(-part->f[iPart]/part->TauP[iPart]*dtStep));
                it->eddyzPartPos[j] += - it->eddyWvel * dtStep + (odtl->wvel[part->iyPos[iPart]] + (part->TauP[iPart] / part->f[iPart] * part->AGz)) * dtStep
                    + part->TauP[iPart] / part->f[iPart] * (part->wvel[iPart] - odtl->wvel[part->iyPos[iPart]] - part->TauP[iPart]/part->f[iPart]*part->AGz)
                    * (1 - exp(-part->f[iPart]/part->TauP[iPart]*dtStep));

            }
        }
    }
}

///////////////////////////////////////////////////////////////////////

void diffuser::getEddyUWvel(double time) {
    std::vector<eddyCombination>::iterator it; 
    for ( it=ActiveEddy.begin() ; it != ActiveEddy.end(); ++it) {
        int iEdLeft  = odtl->linePositionToIndex(it->eddyLeftEdge,true);
        int iEdRight = odtl->linePositionToIndex(it->eddyRightEdge,true);
        it->eddyUvel = odtl->uvel[iEdLeft]*(odtl->posf[iEdLeft+1]-it->eddyLeftEdge) + odtl->uvel[iEdRight]*(it->eddyRightEdge-odtl->posf[iEdRight]);
        it->eddyVvel = odtl->vvel[iEdLeft]*(odtl->posf[iEdLeft+1]-it->eddyLeftEdge) + odtl->vvel[iEdRight]*(it->eddyRightEdge-odtl->posf[iEdRight]);
        it->eddyWvel = odtl->wvel[iEdLeft]*(odtl->posf[iEdLeft+1]-it->eddyLeftEdge) + odtl->wvel[iEdRight]*(it->eddyRightEdge-odtl->posf[iEdRight]);
        for(int i=iEdLeft+1; i<=iEdRight-1; i++) {
            it->eddyUvel += odtl->uvel[i]*(odtl->posf[i+1]-odtl->posf[i]);
            it->eddyVvel += odtl->vvel[i]*(odtl->posf[i+1]-odtl->posf[i]);
            it->eddyWvel += odtl->wvel[i]*(odtl->posf[i+1]-odtl->posf[i]);
        }
        it->eddyUvel = it->eddyUvel/it->eddyLength;
        it->eddyVvel = it->eddyVvel/it->eddyLength;
        it->eddyWvel = it->eddyWvel/it->eddyLength;
    }
}

///////////////////////////////////////////////////////////////////////
/** Check active eddies at certain time step
 * called in diffuser::diffuseodtline  
 * @time: current diffustion time step
 *
 * @Guangyuan Sun 10/2012
 */

void diffuser::checkActiveEddyeffect(double time) {

    if(part->PeddyType == 2 || part->PeddyType == 3) {

        if(ActiveEddy.empty()) return;    

        // --------- find out how many eddies are dead at current time
        std::vector<eddyCombination>::iterator it = ActiveEddy.begin();
        std::vector<eddyCombination>::iterator end = ActiveEddy.begin();
        for ( it=ActiveEddy.begin() ; it < ActiveEddy.end(); it++) {

            if (it->eddyDeadTime <= time && ActiveEddy.size() == 1) {
                    ActiveEddy.resize(0);
                    break;
            }
            if (it->eddyDeadTime <= time) {
                continue;
            }
            else {
                ActiveEddy.erase(ActiveEddy.begin(), it);
                break;
            }
        }
    }
// cout << " ActiveEddy.size() after erase " << ActiveEddy.size() << endl;    
}

///////////////////////////////////////////////////////////////////////
/** calculate gas velocity due to sum effects of active eddies in the process of diffusion (for each particle)
 * called in particles::computeRHSFandSetGasSource
 * @iPart: particle index
 *
 * @Guangyuan Sun 10/2012
 */

void diffuser::getPartEddyVel(double time) {

    std::fill(part->PeddyUvel.begin(), part->PeddyUvel.end(), 0.);
    std::fill(part->PeddyVvel.begin(), part->PeddyVvel.end(), 0.);
    std::fill(part->PeddyWvel.begin(), part->PeddyWvel.end(), 0.);
    
    // ------- TypeC (2) or Type IC (3)
    if(part->PeddyType == 2 || part->PeddyType == 3) {
        std::vector<eddyCombination>::iterator it = ActiveEddy.begin();
        std::vector<eddyCombination>::iterator end = ActiveEddy.begin();

        for ( it=ActiveEddy.begin() ; it < ActiveEddy.end(); it++) {

            for (int j=0; j<it->eddyiPart.size(); j++) {
                int iPart = it->eddyiPart[j];

                bool L = true;

                // ------- if TypeIC (3), get rid of active eddies that has typeI interaction with particle (particle hits bottomline of eddy box)
                // ------- if TypeC  (2), do not need this step
                if(part->PeddyType == 3 && it->eddyiPartIC.size() != 0) {
                    for(int i = 0; i<it->eddyiPartIC.size(); i++){
                        if(iPart == it->eddyiPartIC[i]) {
                            bool L = false;   
                            break;
                        }
                    }
                }

                // ------- sum up the eddy velocity of active eddies 

                if(L == true) {
                    for (int j=0; j<it->eddyiPart.size(); j++) {
                        if (iPart==it->eddyiPart[j]) {
                            part->PeddyVvel[iPart] += it->eddyPvel; //GYSun 
                        }
                    }
                }
                part->PeddyUvel[iPart] = odtl->uvel[part->iyPos[iPart]];
                part->PeddyWvel[iPart] = odtl->wvel[part->iyPos[iPart]];
            }
        }
    }
    // ------- TypeI (1)

    if(part->PeddyType == 1) {
        double fC, gv;
        int    ii;

        for (int iPart=0; iPart<part->nPart; iPart++) {

            if(odtP->Lrxn) {
                ii = part->iyPos[iPart];                                                         // line cell index of particle
                fC         = (part->yPos[iPart] - odtl->posf[ii]) / (odtl->posf[ii + 1] - odtl->posf[ii]);    // fractional location in cell
                gv = Gvel[ii] + fC * (Gvel[ii+1] - Gvel[ii]);                    // gas velocity at particle relative to expansion center
            }
            else {
                gv  = 0;
            }
            part->PeddyVvel[iPart] = gv;
            //        part->PeddyVvel[iPart] = odtl->vvel[part->iyPos[iPart]];
            part->PeddyUvel[iPart] = odtl->uvel[part->iyPos[iPart]];
            part->PeddyWvel[iPart] = odtl->wvel[part->iyPos[iPart]];
        }
    }
}

////////////////////////////////////////////////////////////////////////////
/**
 *  strang splitting of diffusion and chemistry 
 *  in species equation with implicit diffusion and chemistry 
 *  (should be used for premixed combustion to run stable)
 **/
void diffuser::strangSplitting(double dtStep, double time) {

    vector<double> uOld;
    if (odtP->Lspatial)
        uOld = odtl->uvel;

    vector< vector<double> > Y_old(odtl->nspc);
    for (int k = 0; k < odtl->nspc; k++) {
        Y_old[k].resize(odtl->ngrd, 0.0);

        for (int i = 0; i < odtl->ngrd; i++)
            Y_old[k][i] = odtl->yspc[k][i];
    }

    double dtStep2 = dtStep / 2;

    rhsf(k1rhs, P1rhs, dtStep, time);


    //---------- species: do splitting: dt/2 diffusion; dt reaction rate; dt/2 diffusion ; do nspc-1 enforcing sum=1 with yN2
    //---------- dt/2 diffusion--------------------------------------------------
    for (int k = 0; k < iN2; k++) {

        computeImplDiff(k, dtStep2, flxProp);

    }
    for (int k = iN2 + 1; k < odtl->nspc; k++) {

        computeImplDiff(k, dtStep2, flxProp);
    }

    for (int i = 0; i < odtl->ngrd; i++) {


        double sum = 0.0;
        for (int k = 0; k < iN2; k++) {
            sum += odtl->yspc[k][i];
        }
        for (int k = iN2 + 1; k < odtl->nspc; k++) {
            sum += odtl->yspc[k][i];
        }
        odtl->yspc[iN2][i] = 1.0 - sum;



    }
    odtl->enforceYsBounds(Y_old); // neg spec set to 0, renorm., for bad state go back to old state

    //---------- dt reaction rate--------------------------------------------

    //double radSource = 0.0; //  !!!!!  unused variable

    for (int i = 0; i < odtl->ngrd; i++) {

        for (int k = 0; k < odtl->nspc; k++) Y_old[k][i] = odtl->yspc[k][i];

        //double dd = 1.0 / (dxML[i] * odtl->rho[i]); //  !!!!! unused variable


        // 	    if(odtP->Iradiation==1) { // optically thin radiation model
        // 	      vector<double> xMoleSp(odtl->ngrd, 0.0);
        // 	      vector<double> yi(odtl->nspc); // working array
        // 	      odtl->getYspVecAtPt(i, yi);
        // 	      odtl->gas->setState_PY(odtl->pres, &yi[0]);
        // 	      odtl->gas->getMoleFractions(&xMoleSp[0]);
        // 	      rad.opthinRadHeatSourcePt(xMoleSp, odtl->temp[i], odtl->pres, radSource);
        // 	    }


        int ip = i + 1;
#pragma omp critical (cantera)
        {
            brxr.setFluxSources(i, ip, 1.0 / (dxML[i] * odtl->rho[i]), flxProp, iptYspc, iptEnth, (dpdt + radSource_G[i]) / odtl->rho[i]);
            brxr.integrateCell(i, dtStep);

            for (int k = 0; k < odtl->nspc; k++)
                k1rhs[iptYspc + k][i] = brxr.meanRates[k];
        }
        double sum = 0.0;
        for (int k = 0; k < iN2; k++) {
            odtl->yspc[k][i] = odtl->yspc[k][i] + dtStep * k1rhs[iptYspc + k][i];
            sum += odtl->yspc[k][i];
        }
        for (int k = iN2 + 1; k < odtl->nspc; k++) {
            odtl->yspc[k][i] = odtl->yspc[k][i] + dtStep * k1rhs[iptYspc + k][i];
            sum += odtl->yspc[k][i];
        }
        odtl->yspc[iN2][i] = 1.0 - sum;


    }
    odtl->enforceYsBounds(Y_old); // neg spec set to 0, renorm., for bad state go back to old state


    //---------- dt/2 diffusion--------------------------------------------------
    for (int i = 0; i < odtl->ngrd; i++)
        for (int k = 0; k < odtl->nspc; k++) Y_old[k][i] = odtl->yspc[k][i];

    for (int k = 0; k < iN2; k++) {
        computeImplDiff(k, dtStep2, flxProp);

    }
    for (int k = iN2 + 1; k < odtl->nspc; k++) {
        computeImplDiff(k, dtStep2, flxProp);
    }

    for (int i = 0; i < odtl->ngrd; i++) {

        double sum = 0.0;
        for (int k = 0; k < iN2; k++) {
            if (odtl->yspc[k][i] < 0.0) odtl->yspc[k][i] = 0.0;
            sum += odtl->yspc[k][i];
        }
        for (int k = iN2 + 1; k < odtl->nspc; k++) {
            if (odtl->yspc[k][i] < 0.0) odtl->yspc[k][i] = 0.0;
            sum += odtl->yspc[k][i];
        }
        odtl->yspc[iN2][i] = 1.0 - sum;


    }
    odtl->enforceYsBounds(Y_old); // neg spec set to 0, renorm., for bad state go back to old state


    int nlimit;
    nlimit = iptEnth + 1;
    for (int k = 0; k < nlimit; k++) {
        for (int i = 0; i < odtl->ngrd; i++) {
            (*(odtl->props[k]))[i] = (*(odtl->props[k]))[i] + dtStep * k1rhs[k][i];
        }
    }
    //---------- enthalpy 

    //         for(int i=0, ip=1; i < odtl->ngrd; i++, ip++) {
    // 
    //            int dd = 1.0/(dxML[i]*odtl->rho[i]);
    // 	  
    //            k1rhs[iptEnth][i] = -dd*(flxProp[iptEnth][ip] - flxProp[iptEnth][i]) + (dpdt+radSource) / odtl->rho[i];
    //         
    //            odtl->enth[i] = odtl->enth[i] + dtStep*k1rhs[iptEnth][i];
    // 	   
    //         }

    // momentum equation !!



    //---------- pressure

    odtl->pres = odtl->pres + dtStep*dpdt;

    if (odtP->Lrxn || odtP->Lprxn || odtP->ItableLookup)
        updateRhoAndGrid(dtStep, uOld);
}

////////////////////////////////////////////////////////////////////////////
/** C.Schroedinger
 * computation of implicit diffusion terms with Thomas-algorithm
 * @param k input: species index.
 * @param dtStep input: diffusion time step.
 */

void diffuser::computeImplDiff(int k, double dtStep, std::vector<std::vector<double> > &flxProp) {

    int i, im, ip;
    double dd, dd1, dd2;

    double rho_f1; // dummy face values
    double M_f1;
    double rho_f2; // dummy face values
    double M_f2;


    vector<double> A(odtl->ngrd);
    vector<double> B(odtl->ngrd);
    vector<double> C(odtl->ngrd);
    vector<double> D(odtl->ngrd);
    vector<double> Y_old(odtl->ngrd);



    mmw = odtl->getMMW();

    for (i = 0, ip = 1; i < odtl->ngrd; i++, ip++) {

        Y_old[i] = odtl->yspc[k][i];

        dd = 1.0 / (dxML[i] * odtl->rho[i]);

        D[i] = -0.5 * dd * (flxProp[iptYspc + k][ip] - flxProp[iptYspc + k][i]) * dtStep + odtl->yspc[k][i];


        rho_f1 = linearInterpToFace(i, *odtl, odtl->rho);
        M_f1 = linearInterpToFace(i, *odtl, mmw);
        rho_f2 = linearInterpToFace(ip, *odtl, odtl->rho);
        M_f2 = linearInterpToFace(ip, *odtl, mmw);



        if (i == 0) {

            dd2 = 2.0 * rho_f2 / (dxML[i] + dxML[ip]);

            A[i] = 0.0;

            B[i] = (-dd2 * DmixYs_f[k][ip]-(1 - 1 / (odtl->pos[ip] - odtl->pos[i])*(odtl->posf[ip] - odtl->pos[i]))*(mmw[ip] - mmw[i]) / M_f2) / (dxML[i] * odtl->rho[i]);

            C[i] = (dd2 * DmixYs_f[k][ip]-(1 / (odtl->pos[ip] - odtl->pos[i])*(odtl->posf[ip] - odtl->pos[i]))*(mmw[ip] - mmw[i]) / M_f2) / (dxML[i] * odtl->rho[i]);
            // 
            B[i] = -0.5 * dtStep * B[i] + 1;
            C[i] *= -0.5 * dtStep;
        } else if (i == odtl->ngrd - 1) {

            im = i - 1;
            dd1 = 2.0 * rho_f1 / (dxML[im] + dxML[i]);

            C[i] = 0.0;
            A[i] = (dd1 * DmixYs_f[k][i]+(1 - 1 / (odtl->pos[i] - odtl->pos[im])*(odtl->posf[i] - odtl->pos[im]))*(mmw[i] - mmw[im]) / M_f1) / (dxML[i] * odtl->rho[i]);
            B[i] = (-dd1 * DmixYs_f[k][i]+(1 / (odtl->pos[i] - odtl->pos[im])*(odtl->posf[i] - odtl->pos[im]))*(mmw[i] - mmw[im]) / M_f1) / (dxML[i] * odtl->rho[i]);

            A[i] *= -0.5 * dtStep;
            B[i] = -0.5 * dtStep * B[i] + 1;

        } else {

            im = i - 1;

            dd1 = 2.0 * rho_f1 / (dxML[im] + dxML[i]);
            dd2 = 2.0 * rho_f2 / (dxML[i] + dxML[ip]);

            A[i] = (dd1 * DmixYs_f[k][i]+(1 - 1 / (odtl->pos[i] - odtl->pos[im])*(odtl->posf[i] - odtl->pos[im]))*(mmw[i] - mmw[im]) / M_f1) / (dxML[i] * odtl->rho[i]);

            B[i] = (-dd2 * DmixYs_f[k][ip]-(1 - 1 / (odtl->pos[ip] - odtl->pos[i])*(odtl->posf[ip] - odtl->pos[i]))*(mmw[ip] - mmw[i]) / M_f2) / (dxML[i] * odtl->rho[i])+ (-dd1 * DmixYs_f[k][i]+(1 / (odtl->pos[i] - odtl->pos[im])*(odtl->posf[i] - odtl->pos[im]))*(mmw[i] - mmw[im]) / M_f1) / (dxML[i] * odtl->rho[i]);

            C[i] = (dd2 * DmixYs_f[k][ip]-(1 / (odtl->pos[ip] - odtl->pos[i])*(odtl->posf[ip] - odtl->pos[i]))*(mmw[ip] - mmw[i]) / M_f2) / (dxML[i] * odtl->rho[i]);

            A[i] = A[i]*(-0.5) * dtStep;
            B[i] = -0.5 * dtStep * B[i] + 1;
            C[i] = C[i]*(-0.5) * dtStep;
        }
    }



    //Thomas-algorithm

    /* Modify the coefficients. */
    C[0] /= B[0]; /* Division by zero risk. */
    D[0] /= B[0]; /* Division by zero would imply a singular matrix. */
    for (i = 1; i < odtl->ngrd; i++) {
        double id = (B[i] - C[i - 1] * A[i]); /* Division by zero risk. */
        C[i] /= id; /* Last value calculated is redundant. */
        D[i] = (D[i] - D[i - 1] * A[i]) / id;
    }

    /* Now back substitute. */
    odtl->yspc[k][odtl->ngrd - 1] = D[odtl->ngrd - 1];
    for (i = odtl->ngrd - 2; i >= 0; i--)
        odtl->yspc[k][i] = D[i] - C[i] * odtl->yspc[k][i + 1];



}


///////////////////////////////////////////////////////////////////////////////
/** 
 *  
 *  
 */
#ifdef IMPLICIT
void diffuser::computeMatrix(vector<vector<double> > &A, vector<vector<double> > &B,
        vector<vector<double> > &C, vector<vector<double> > &rhs, double &dtStep){
    vector<double> xipMxi; // = x_{i+1} - x_i
    xipMxi = vector<double>(odtl->ngrd-1, 0.0);
    double CN = 1.0;
    //CN = 2.0 // factor for Crank-Nicolson
    double tempRHS  = 0.0;
    double fracMu   = 0.0, fracDx    = 0.0, fracLamb = 0.0;
    double prodFrac = 0.0, prodFracL = 0.0;
    double tempA    = 0.0, tempB     = 0.0, tempC = 0.0;
    double tempAL   = 0.0, tempBL    = 0.0;
    int nlimit = 3;
    if (odtP->LhasTemp) nlimit++;
    
    if (odtP->bcType == 1){
        *proc.ostrm << endl << "WARNING:\n" << "Periodic boundaries are currently not supported by the implicit solver!\n";
        exit(0);
    }
    if (odtP->Lrxn){
        *proc.ostrm << endl << "WARNING:\n" << "Reactions are currently not supported by the implicit solver!\n";
        exit(0);
    }
    if (odtP->Iparticles){
        *proc.ostrm << endl << "WARNING:\n" << "Particles are currently not supported by the implicit solver!\n";
        exit(0);
    }
    if (odtP->Lspatial){
        *proc.ostrm << endl << "WARNING:\n" << "The spatial formulation is currently not supported by the implicit solver!\n";
        exit(0);
    }
    
    
    for (int i = 0, ip = 1; i < odtl->ngrd; i++, ip++){
        
        // u velocity,  default BC is Dirichlet
        if (i == 0){
            xipMxi[i] = odtl->pos[ip] - odtl->pos[i];
            A[0][i] = 0.0;
            C[0][i] = -dtStep * visc_f[ip] / ( CN * odtl->rho[i] * dxML[i] * xipMxi[i] );
            B[0][i] = -C[0][i] + dtStep * visc_f[i] / ( CN * odtl->rho[i] * dxML[i] * (odtl->pos[i]-odtl->posf[i]) );
            tempRHS = dtStep * visc_f[i] / ( CN * odtl->rho[i] * dxML[i] * (odtl->pos[i]-odtl->posf[i]) );
            rhs[0][i] = odtl->bcprops[0][1] * tempRHS;
            rhs[1][i] = odtl->bcprops[1][1] * tempRHS;
            rhs[2][i] = odtl->bcprops[2][1] * tempRHS;
            if (odtP->LhasTemp){
                A[3][i] = 0.0;
                C[3][i] = -dtStep * lambda_f[ip] / ( CN * odtl->rho[i] * dxML[i] * xipMxi[i] );
                B[3][i] = -C[3][i] + dtStep * lambda_f[i] / ( CN * odtl->rho[i] * dxML[i] * (odtl->pos[i]-odtl->posf[i]) );
                rhs[3][i] = odtl->bcprops[3][1] * tempRHS * lambda_f[i] / visc_f[i];
            }
        }
        else if (i == odtl->ngrd -1) {
            A[0][i] = -dtStep * visc_f[i]  / ( CN * odtl->rho[i] * dxML[i] * xipMxi[i-1] );
            C[0][i] = 0.0;
            B[0][i] = -A[0][i] + dtStep * visc_f[ip] / ( CN * odtl->rho[i] * dxML[i] * (odtl->posf[ip]-odtl->pos[i]) );
            tempRHS = dtStep * visc_f[ip] / ( CN * odtl->rho[i] * dxML[i] * (odtl->posf[ip]-odtl->pos[i]) );
            rhs[0][i] = odtl->bcprops[0][5] * tempRHS;
            rhs[1][i] = odtl->bcprops[1][5] * tempRHS;
            rhs[2][i] = odtl->bcprops[2][5] * tempRHS;
            if (odtP->LhasTemp){
                A[3][i] = -dtStep * lambda_f[i]  / ( CN * odtl->rho[i] * dxML[i] * xipMxi[i-1] );
                C[3][i] = 0.0;
                B[3][i] = -A[3][i] + dtStep * lambda_f[ip] / ( CN * odtl->rho[i] * dxML[i] * (odtl->posf[ip]-odtl->pos[i]) );
                rhs[3][i] = odtl->bcprops[3][5] * tempRHS * lambda_f[i] / visc_f[i];
            }
        }
        else{
            xipMxi[i] = odtl->pos[ip] - odtl->pos[i];
            A[0][i] = -dtStep * visc_f[i]  / ( CN * odtl->rho[i] * dxML[i] * xipMxi[i-1] );
            C[0][i] = -dtStep * visc_f[ip] / ( CN * odtl->rho[i] * dxML[i] * xipMxi[i] );
            B[0][i] = -A[0][i] -C[0][i];
            if (odtP->LhasTemp){
                A[3][i] = -dtStep * lambda_f[i]  / ( CN * odtl->rho[i] * dxML[i] * xipMxi[i-1] );
                C[3][i] = -dtStep * lambda_f[ip] / ( CN * odtl->rho[i] * dxML[i] * xipMxi[i] );
                B[3][i] = -A[3][i] -C[3][i];
            }
        }
        
        // multi phase
        if (i != 0 && odtl->phase[i-1] != odtl->phase[i]){ // previous cell has different phase
            // the first cell has no previous cell. For periodic boundaries use
            // the next part for changes caused by boundaries.
            
            fracMu = odtl->molec[i-1] / odtl->molec[i];
            tempA = ( 3*dxML[i] +dxML[ip] ) * ( dxML[i] +dxML[ip] ) / ( (2 *dxML[i] +dxML[ip]) * (2 *dxML[i] +dxML[ip]) );
                 // +fracMu * dxML[i] * (dxML[i] +dxML[ip] ) / ( dxML[i-1] * (2 *dxML[i] +dxML[ip]) );
            tempB =  fracMu * dxML[i] * (dxML[i] +dxML[ip] ) / ( dxML[i-1] * (2 *dxML[i] +dxML[ip]) );
            tempA += tempB;
            tempC = dxML[i] * dxML[i] / ( (2 *dxML[i] +dxML[ip]) * (2 *dxML[i] +dxML[ip]) );
            B[0][i] = -C[0][i] - 2/CN * odtl->molec[i] * (1-1/tempA) / ( odtl->rho[i] * dxML[i] * dxML[i]);
            C[0][i] -= 2/CN * odtl->molec[i] * tempC / ( tempA * odtl->rho[i] * dxML[i] * dxML[i] );
            A[0][i] = 2/CN * odtl->molec[i] * tempB / ( tempA * odtl->rho[i] * dxML[i] * dxML[i] );
            if (odtP->LhasTemp){
                fracMu = odtl->lambda[i-1] / odtl->lambda[i];
                tempAL = ( 3*dxML[i] +dxML[ip] ) * ( dxML[i] +dxML[ip] ) / ( (2 *dxML[i] +dxML[ip]) * (2 *dxML[i] +dxML[ip]) );
                     //  +fracMu * dxML[i] * (dxML[i] +dxML[ip] ) / ( dxML[i-1] * (2 *dxML[i] +dxML[ip]) );
                tempBL =  fracMu * dxML[i] * (dxML[i] +dxML[ip] ) / ( dxML[i-1] * (2 *dxML[i] +dxML[ip]) );
                tempAL += tempBL;
                B[3][i] = -C[3][i] - 2/CN * odtl->lambda[i] * (1-1/tempA) / ( odtl->rho[i] * dxML[i] * dxML[i]);
                C[3][i] -= 2/CN * odtl->lambda[i] * tempC / ( tempAL * odtl->rho[i] * dxML[i] * dxML[i] );
                A[3][i] = 2/CN * odtl->lambda[i] * tempBL / ( tempAL * odtl->rho[i] * dxML[i] * dxML[i] );
            }
            
            ////A[0][i] = -2/CN * odtl->molec[i-1] / (odtl->rho[i] * dxML[i-1] * dxML[i] * prodFrac);
            ////B[0][i] = A[0][i] *fracMu *fracDx - C[0][i];
            ////A[0][i] += 2/CN * odtl->molec[i-1] / (odtl->rho[i-1] * dxML[i-1] * dxML[i]);
            //A[0][i] = 2/CN * odtl->molec[i] / (odtl->rho[i] * dxML[i] * dxML[i] * prodFrac);
            //B[0][i] = -C[0][i] +A[0][i] *fracMu *fracDx - 2/CN *odtl->molec[i] /(odtl->rho[i] *dxML[i] *dxML[i]);
            //if (odtP->LhasTemp){
            //    //A[3][i] = - 2/CN * odtl->lambda[i-1] / (odtl->rho[i] * dxML[i-1] * dxML[i] * prodFracL);
            //    //B[3][i] = A[3][i] *fracLamb *fracDx - C[3][i];
            //    //A[3][i] += 2/CN * odtl->lambda[i-1] / (odtl->rho[i-1] * dxML[i-1] * dxML[i]);
            //    A[3][i] = 2/CN * odtl->lambda[i] / (odtl->rho[i] * dxML[i] * dxML[i] * prodFracL);
            //    B[3][i] = -C[3][i] +A[3][i] *fracLamb *fracDx - 2/CN *odtl->lambda[i] /(odtl->rho[i] *dxML[i] *dxML[i]);
            //}
        }
        if (i != odtl->ngrd-1 && odtl->phase[i] != odtl->phase[ip]){ // next cell has different phase
            // the last cell has no next cell
            
            fracMu = odtl->molec[ip] / odtl->molec[i];
            tempA  = ( 3*dxML[i] +dxML[i-1] ) * ( dxML[i] +dxML[i-1] ) / ( (2 *dxML[i] +dxML[i-1]) * (2 *dxML[i] +dxML[i-1]) );
                  // +fracMu * dxML[i] * (dxML[i] +dxML[i-1] ) / ( dxML[ip] * (2 *dxML[i] +dxML[i-1]) );
            tempB  =  fracMu * dxML[i] * (dxML[i] +dxML[i-1] ) / ( dxML[ip] * (2 *dxML[i] +dxML[i-1]) );
            tempA += tempB;
            tempC  = dxML[i] * dxML[i] / ( (2 *dxML[i] +dxML[i-1]) * (2 *dxML[i] +dxML[i-1]) );
            B[0][i] = -A[0][i] - 2/CN * odtl->molec[i] * (1-1/tempA) / ( odtl->rho[i] * dxML[i] * dxML[i]);
            A[0][i] -= 2/CN * odtl->molec[i] * tempC / ( tempA * odtl->rho[i] * dxML[i] * dxML[i] );
            C[0][i] = 2/CN * odtl->molec[i] * tempB / ( tempA * odtl->rho[i] * dxML[i] * dxML[i] );
            if (odtP->LhasTemp){
                fracMu = odtl->lambda[ip] / odtl->lambda[i];
                tempAL = ( 3*dxML[i] +dxML[i-1] ) * ( dxML[i] +dxML[i-1] ) / ( (2 *dxML[i] +dxML[i-1]) * (2 *dxML[i] +dxML[i-1]) );
                      // +fracMu * dxML[i] * (dxML[i] +dxML[i-1] ) / ( dxML[ip] * (2 *dxML[i] +dxML[i-1]) );
                tempBL =  fracMu * dxML[i] * (dxML[i] +dxML[i-1] ) / ( dxML[ip] * (2 *dxML[i] +dxML[i-1]) );
                tempAL += tempBL;
                B[3][i] = -A[3][i] - 2/CN * odtl->lambda[i] * (1-1/tempAL) / ( odtl->rho[i] * dxML[i] * dxML[i]);
                A[3][i] -= 2/CN * odtl->lambda[i] * tempC / ( tempAL * odtl->rho[i] * dxML[i] * dxML[i] );
                C[3][i] = 2/CN * odtl->lambda[i] * tempBL / ( tempAL * odtl->rho[i] * dxML[i] * dxML[i] );
            }
            
            //// old things; doesn't work
            //fracMu   = odtl->molec[ip] / odtl->molec[i];
            //fracDx   = dxML[i] / dxML[ip];
            //prodFrac = 1 + fracMu * fracDx;
            //C[0][i] = 2/CN * odtl->molec[i] / (odtl->rho[i] * dxML[i] * dxML[i] * prodFrac);
            //B[0][i] = -A[0][i] +C[0][i] -2/CN *odtl->molec[i] / (odtl->rho[i] *dxML[i] *dxML[i]);
            //C[0][i] *= fracMu * fracDx;
            //if (odtP->LhasTemp){
            //    fracLamb  = odtl->lambda[ip] / odtl->lambda[i];
            //    prodFracL = 1 + fracLamb * fracDx;
            //    C[3][i] = 2/CN * odtl->lambda[i] / (odtl->rho[i] * dxML[i] * dxML[i] * prodFracL);
            //    B[3][i] = -A[3][i] +C[3][i] -2/CN *odtl->lambda[i] / (odtl->rho[i] *dxML[i] *dxML[i]);
            //    C[3][i] *= fracLamb * fracDx;
            //}
        }
        
        B[0][i] += 1.0;
        if (odtP->LhasTemp)
            B[3][i] += 1.0;
        
        // v and w velocity
        A[1][i] = A[0][i]; A[2][i] = A[0][i];
        B[1][i] = B[0][i]; B[2][i] = B[0][i];
        C[1][i] = C[0][i]; C[2][i] = C[0][i];
    }
    
    // reactions
    if (odtP->Lrxn){
        // TBS
    }
    
    //----- adjustment for other boundary conditions
    double err = 0.000001;
    int    N   = odtl->ngrd-1;
    for (int k = 0; k < nlimit; k++){
        // lower boundary
        if (abs(odtl->bcprops[k][0] - 1.0) > err){
            // adjust lower boundary
            if (abs(odtl->bcprops[k][0] - 0.0) < err){
                // periodic bc
                *proc.ostrm << endl << "ERROR:\nPeriodic boundaries are currently not supported by the implicit solver.";
                *proc.ostrm << endl << "odtl->bcprops[" << k << "][0] = " << odtl->bcprops[k][0] << endl;
                exit(0);
            }
            else if (abs(odtl->bcprops[k][0] - 2.0) < err){
                // Neuman bc
                B[k][0]   = -C[k][0] + 1.0;
                rhs[k][0] = -dtStep * odtl->bcprops[k][1] * visc_f[0] / (CN * odtl->rho[0] * dxML[0]);
                if (odtP->LhasTemp && k==3) rhs[k][0] = rhs[k][0] * lambda_f[0] / visc_f[0];
            }
            else if (abs(odtl->bcprops[k][0] - 3.0) < err){
                // Cauchy bc
                // TBS
            }
            else if (abs(odtl->bcprops[k][0] - 4.0) < err){
                // Robin bc
                // TBS
            }
            else{
                *proc.ostrm << endl << "ERROR:\n" << "boundary condition could not be interpreted";
                *proc.ostrm << endl << "odtl->bcprops[" << k << "][0] = " << odtl->bcprops[k][0] << endl;
                exit(0);
            }
        }
        
        // upper boundary
        if (abs(odtl->bcprops[k][4] - 1.0) > err){
            // adjust upper boundary
            if (abs(odtl->bcprops[k][0] - 0.0) < err){
                // periodic bc
                *proc.ostrm << endl << "ERROR:\nPeriodic boundaries are currently not supported by the implicit solver.";
                *proc.ostrm << endl << "odtl->bcprops[" << k << "][4] = " << odtl->bcprops[k][4] << endl;
                exit(0);
            }
            else if (abs(odtl->bcprops[k][4] - 2.0) < err){
                // Neuman bc
                B[k][N]   = -A[k][N] + 1.0;
                rhs[k][N] =  dtStep *odtl->bcprops[k][5] * visc_f[N+1] / (CN * odtl->rho[N] * dxML[N]);
                if (odtP->LhasTemp && k==3) rhs[k][N] = rhs[k][N] * lambda_f[N+1] / visc_f[N+1];
            }
            else if (abs(odtl->bcprops[k][4] - 3.0) < err){
                // Cauchy bc
                // TBS
            }
            else if (abs(odtl->bcprops[k][4] - 4.0) < err){
                // Robin bc
                // TBS
            }
            else{
                *proc.ostrm << endl << "ERROR:\n" << "boundary condition could not be interpreted";
                *proc.ostrm << endl << "odtl->bcprops[" << k << "][4] = " << odtl->bcprops[k][4] << endl;
                exit(0);
            }
        }
    }
}
#endif


/////////////////////////////DOXYGEN DOCUMENTATION//////////////////////////////


/*! \fn void diffuser::computeFluxes()
 * uvw velocities: \fun{\tau = -\mu \nabla  u} \n
 * species: \fun{ j_k = -\rho*D_k*\nabla (Y_k) - \rho*D_k*\frac{Y_k}{M_k} * \nabla (M)} \n
 * enthalpy: \fun{ q = -\lambda*\nabla T  + \displaystyle\sum\limits_k( h_k*j_k )}               \n
 */

/*! \fn double diffuser::setGasVelocity_return_dpdt()
 *
 * Sets \fun{G_{vel}} member
 * \f[
 *       div (v) = \frac{-1}{\gamma P} \frac{ dp}{dt } + \mho
 * \f]
 * Integrate this over each cell:
 *
 * \f[
 *       \Delta G_{vel} = v_e-v_w = \frac{-1}{P} \frac{dp}{dt} \frac{\Delta x}{\gamma} + \int{\mho dx}
 * \f]
 *
 * To evaluate:
 * First, find \fun{\int{\mho dx}} in each cell and store it in \fun{G_{vel}}.\n
 * For open domains, this is \fun{\Delta G_{vel}} and \fun{\frac{dp}{dt} = 0}\n
 * For closed domains, evaluate dp/dt by integrating the first eqn over the whole domain:
 *
 * \f[
 *       v_E - v_W = -\frac{1}{P} \frac{dp}{dt} \sum{\frac{\Delta x}{\gamma}} + \sum{\int{\mho dx}}
 * \f]
 *
 * where sums are over all cells, and \fun{v_E}, \fun{v_W} are currently zero (fixed
 * domain boundaries).
 * Then
 *   \f[
 *       \frac{dp}{dt} = \frac{\sum{\int{\mho dx}} } { \frac{1}{P} \sum{\frac{dx}{\gamma} }}
 *   \f]
 * Now, given \fun{\frac{dp}{dt}}, we can evaluate \fun{\Delta G_{vel}} in each cell be appending the \fun{-\frac{1}{P} \frac{dp}{dt} \frac{\Delta x}{\gamma}} term.
 *
 * Then, given \fun{\Delta G_{vel}}, stored in \fun{G_{vel}}, we can compute \fun{G_{vel}} on cell faces, relative to the
 *     zero velocity expansion center.
 */


////////////////////////////END DOXYGEN DOCUMENTATION///////////////////////////////////////
