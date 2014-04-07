/**
 * @file batchRXR.h
 * Header file for class batchRXR
 */

#ifndef BATCHRXR_H
#define BATCHRXR_H

#include "odtParam.h"
#include "odtline.h"

#include "cvode/cvode.h"
#include "cvode/nvector_serial.h"
#include "cvode/cvode_dense.h"
#include "cvode/sundials_dense.h"
#include "cvode/sundials_types.h"

////////////////////////////////////////////////////////////////////////////////////////

/** A class to integrate cell chemistry implicitly.  Allows explicit time advancement
 *  with stiff gas chemistry.  Integrate cells as for a batch reactor (but with
 *  additional flux terms.  Compute the average rate over the time advancement).
 *  Solve \cond
 *        dY/dt = -1/(rho*dx)*(j_e - j_w) + RR/rho,
 *        dh/dt = -1/(rho*dx)*(q_e - q_w)
 *        \endcond
 *  
 *  \f[   
 *        \frac{dY}{dt} = - \frac{1}{\rho \times dx} \left( j_{e} - j_{w} \right) + \frac{RR}{\rho} \\
 *        \frac{dh}{dt} = - \frac{1}{\rho \times dx} \left( q_{e} - q_{w} \right) 
 *  \f]
 *
 *  The first term on RHS is the flux source and is assumed constant.  Hence
 *  the enthalpy solve is trivial with h = h_init + t * fluxTerm.
 *  The mean (RR/rho) is computed as (Y2-Y1)/dt - fluxTerm.
 *  Note, the mean rate (dY/dt)_mean could be computed as (Y2-Y1)/dt also.
 *  Here, Y1 and Y2 denote mass fractions at times 1 and 2, respectively.
 *
 *  USER calls setFluxSources, then integrateCell, then use member meanRates
 *  
 *  @author David O. Lignell
 */

class batchRXR {

    public :

        ////////////////////// DATA MEMBERS /////////////////////

        std::vector<double> meanRates;       ///< mean spc rate

        int                 neq;             ///< number of eqns solved
        odtParam            *odtP;           ///< pointer to odt param object
        odtline             *odtl;           ///< pointer to combustion line

        double              atol;            ///< CVODE tol
        double              rtol;            ///< CVODE tol

        void                *cvode_mem;      ///< CVode memory

        N_Vector            dVar;            ///< composition vector for CVode

        double              heatFluxSource;  ///< \fun{-1/(\rho*dx) * (q_e - q_w)}
        std::vector<double> massFluxSources; ///< \fun{-1/(\rho*dx) * (j_e - j_w)}

        double              hCellStart;      ///< initial enthalpy in cell 

        bool                upToDate;        ///< flag for fluxSources

        int                 iCell;           ///< which cell are we integrating

        std::vector<double>  yd;             ///< species array dummy
        std::vector<double>  rd;             ///< species rxn rate array dummy

        ////////////////////// MEMBER FUNCTIONS  /////////////////////

        void integrateCell(int iCell, double tres);
        void setFluxSources(int i, int ip, double inverseRhoDx, std::vector<std::vector<double> > &flxProp, 
                            int iptYspc, int iptEnth, double dpdtTerm); 

        double getEnthAtTime(double t);
        void testCVflag(int flag, std::string func);

        ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

        batchRXR(odtline *odtlp, odtParam *odtpp);   // constructor
        ~batchRXR();                                 // destructor

};


#endif
