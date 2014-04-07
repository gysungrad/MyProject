/**
 * @file batchRXR.cc
 * Source file for class batchRXR
 */

#include <cstdlib>

#include "batchRXR.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
// Declare the global function prototypes so they can be used in this source file

static int f_CONP(double t, N_Vector y, N_Vector ydot, void* f_data);
void       getProblemSpecificRR(double rho, double temp, double pres, double *yi, double *rr); 

///////////////////////////////////////////////////////////////////////////////

/** Constructor function.
 *
 *  @param odtlp \input combustion line object, set pointer with.
 *  @param odtpp \input odt line object, set pointer with.
 */
batchRXR::batchRXR(odtline *odtlp, odtParam *odtpp) {

    //---------- set pointers

    odtl = odtlp;
    odtP = odtpp;

    if(!odtP->Lrxn)
        return;

    //---------- set some CVode params

    atol = 1.0E-10;
    rtol = 1.0E-4;
    neq = odtl->nspc-1;           // nspc species (h trivial, not included)

    //---------- initialize the dependent variable

    dVar = N_VNew_Serial(neq);

    //---------- set the CVode object

    int flag;

    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);

    if(!cvode_mem) {
        cout << endl << "ERROR INITIALIZING CVODE MEMORY" << endl;
        exit(0);
    }
    
    flag = CVodeMalloc(cvode_mem, f_CONP, 0.0, dVar, CV_SS, rtol, &atol);
            testCVflag(flag, "CVodeMalloc");

    flag = CVDense(cvode_mem, neq);
            testCVflag(flag, "CVDense");

    flag = CVodeSetFdata(cvode_mem, this);
            testCVflag(flag, "CVodeSetFdata");

    flag = CVodeSetMaxNumSteps(cvode_mem, 2000);
            testCVflag(flag, "CVodeSetMaxNumSteps");

    //---------- other

    upToDate = false;
    meanRates       = vector<double>(odtl->nspc, 0.0);
    massFluxSources = vector<double>(odtl->nspc, 0.0);
    yd              = vector<double>(odtl->nspc, 0.0);
    rd              = vector<double>(odtl->nspc, 0.0);

}

///////////////////////////////////////////////////////////////////////////////

/** destructor */
batchRXR::~batchRXR() { 

    if(odtP->Lrxn) {
        N_VDestroy_Serial(dVar);
        CVodeFree(&cvode_mem);
    }
    odtP=0; 
    odtl=0; 
}  

///////////////////////////////////////////////////////////////////////////////

/** Public interface, user sets flux source terms before each call to integrate.
 *
 *  @param i            \input which combustion line cell.
 *  @param ip           \input next combustion line cell (usually i+1).
 *  @param inverseRhoDx \input \fun{1/(\rho*dx)}.
 *  @param flxProp      \input heat and mass fluxes
 *  @param iptYspc      \input pointer to species fluxes in flxProp
 *  @param iptEnth      \input pointer to heat fluxe in flxProp
 *  @param dpdtTerm     \input pressure term in energy equation
 *
 *  @sideeffect sets upToDate to true \n
 *              updates batchRXR::heatFluxSource \n
 *              updates batchRXR::massFluxSources
 */
void batchRXR::setFluxSources(int i, int ip, double inverseRhoDx, std::vector<std::vector<double> > &flxProp, 
                    int iptYspc, int iptEnth, double dpdtTerm) {

    upToDate = true;

    heatFluxSource = -inverseRhoDx * (flxProp[iptEnth][ip] - flxProp[iptEnth][i]) + dpdtTerm;
    
    for(int k=0; k < odtl->nspc; k++) 
        massFluxSources[k] = -inverseRhoDx * (flxProp[iptYspc+k][ip] - flxProp[iptYspc+k][i]);

}

///////////////////////////////////////////////////////////////////////////////

/** Public interface: Integrates a given cell (iC) in combustion line for time tres.  
 *  The output is the meanRates vector.  The integration is done with 
 *  CVODE.
 * 
 *  @param iC \input integrate this cell in the combustion line.
 *  @param tres \input time to integrate for.
 */
void batchRXR::integrateCell(int iC, double tres) {


    iCell = iC;

    //---------- Check that setFluxSources has been called first

    if(!upToDate) {
        cout << "\n\n********* ERROR in batchRXR, fluxes not set" << endl << endl;
        exit(0);
    }
    
    hCellStart = odtl->enth[iCell];

    //---------- Initialize the gas object
    odtl->getYspVecAtPt(iCell, yd);
    odtl->gas->setState_PY(odtl->pres, &yd[0]);
    odtl->gas->setState_HP(odtl->enth[iCell], odtl->pres, 1.E-10);

    //---------- Initialize the Dependent Variable

    for(int i=0; i < neq; i++) 
        NV_Ith_S(dVar, i) = odtl->gas->massFraction(i);

    //---------- Reset CVode for the new cell

    int flag = CVodeReInit(cvode_mem, f_CONP, 0.0, dVar, CV_SS, rtol, &atol); 
               testCVflag(flag, "CVodeReInit");

    //---------- Integrate the solution
    
    double t;
    flag = CVodeSetStopTime(cvode_mem, tres);            testCVflag(flag, "CVodeSetStopTime");
    flag = CVode(cvode_mem, tres, dVar, &t, CV_NORMAL);  testCVflag(flag, "CVode");

    //----------- Recover the solution
    // composition should be normalized already as part of solution
    // dolcheck (seems inefficient)

    yd[neq] = 1.0;
    for(int k=0; k < neq; k++) {
        yd[k]    = NV_Ith_S(dVar,k);
        yd[neq] -= yd[k];
    }
    odtl->gas->setMassFractions_NoNorm(&yd[0]);
    odtl->gas->setState_HP(getEnthAtTime(tres), odtl->pres, 1.E-10);   

    //----------- Compute mean species reaction rate
    // dY/dt = Flux + Rxn
    // [(Y2-Y1]/dt)_mean = (Flux + Rxn)_mean = Flux_mean + Rxn_mean
    // but Flux term is constant here, and solve for Rxn_mean

    for(int k=0; k < odtl->nspc; k++)
        meanRates[k] = ( yd[k] - odtl->yspc[k][iCell] )/tres - 
                              massFluxSources[k];

    //----------- reset flag

    upToDate = false;

}

///////////////////////////////////////////////////////////////////////////////

/** Compute enthalpy in RXR at time t from start of integration.  Enthalpy
 *  variation is by constant flux with trivial ODE that is solved analytically
 *  to give the enthalpy.
 *
 *  @param t \input time from start of integration to compute enthalpy.
 *  @return enthalpy at given time.
 */
double batchRXR::getEnthAtTime(double t) {

 return hCellStart + t * heatFluxSource;

}

///////////////////////////////////////////////////////////////////////////////

/** RHS function for batch reactor integration.  This is a constant pressure problem. 
 *  f_data is cast to the batchRXR that calls cvode (so cvode is passed pointer "this").
 *  "static" means only defined in this file (just following cvode standard use).
 *
 *  @param t input, current time (from start of integration).
 *  @param y input, current dependent variable state (mass fractions at t during integration).
 *  @param ydot output, rhs of dy/dt = ydot.  This is the rate for the time integration.
 *  @param f_data input, pointer to be cast to user data for computing ydot from y.
 */
static int f_CONP(double t, N_Vector y, N_Vector ydot, void* f_data) {


    batchRXR *brxr;
    brxr = static_cast<batchRXR*>(f_data);

    int N = brxr->neq;

    brxr->yd[N] = 1.0;
    for(int k=0; k < N; k++)  {
        brxr->yd[k] = NV_Ith_S(y,k); 
        brxr->yd[N] -= brxr->yd[k];
    }
    brxr->odtl->gas->setState_PY(brxr->odtl->gas->pressure(), &brxr->yd[0]);
    brxr->odtl->gas->setState_HP(brxr->getEnthAtTime(t), brxr->odtl->gas->pressure(),1.E-10);
#ifdef PROBLEMSPECIFICRR
    getProblemSpecificRR(brxr->odtl->gas->density(), brxr->odtl->gas->temperature(), 
                         brxr->odtl->pres, &brxr->yd[0], &brxr->rd[0]);
#else
    brxr->odtl->gas->getNetProductionRates(&brxr->rd[0]);
#endif
    double rhoInv = 1.0/brxr->odtl->gas->density();
    for(int k=0; k < N; k++)
        NV_Ith_S(ydot,k) = brxr->rd[k] * rhoInv * brxr->odtl->gas->molecularWeight(k) + 
                           brxr->massFluxSources[k];

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
/** Message for bad cvode flag 
 * @param flag \input value flag
 * @param func \input function with error
 */

void batchRXR::testCVflag(int flag, std::string func) {

    if(flag != CV_SUCCESS) {
        cout << endl << "ERROR in " << func << endl;
        exit(0);
    }

}





