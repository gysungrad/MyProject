/**
 * @file stats.cc
 * Header file for class stats
 */

#include "stats.h"
#include "anyline.h"
#include "odtline.h"
#include "odtParam.h"
#include <iomanip>
#include <fstream>
#include <cmath>          // fabs
#include <stdlib.h>

#ifdef NEWSTATS
#include <stdio.h>
#include <iostream>
#include "processor.h"

//#ifdef NETCDF // tested on Linux for NetCDF 3.6.3
//#include "netcdfcpp.h"
//#endif

extern processor proc;
#endif

using namespace std;


///////////////////////////////////////////////////////////////////////////////
/**Constructor function (The one to use).
 * Stats has its own uniform grid and does minimal statistics gathering.
 * It computes the mean and mean square values of velocity and mixture fraction.
 * The stats grid goes from 0 to Ldomain, even for periodic domains, so there
 * is some wrapping to do when the domain is periodic.
 * vTrans is just a dummy vector for transfering from the odt grid to the
 * stats grid.
 *
 * @param odtpp \input parameters object to set pointer.
 * @param Ld    \input domain length.
 * @param npts  \input number of evenly spaced stats points along domain.
 */
stats::stats(odtParam &odtpp, double Ld, int npts)
//: edstat ( vector<vector<vector<vector<double> > > >
//            (8, vector<vector<vector<double> > >
//            (4, vector<vector<double> >
//            (odtP->nStat, vector<double>(ngrd, 0.0) ) ) ) )
{


    odtP = &odtpp;
    
    Ldomain = Ld;
    ngrd    = npts;
    ngrdf   = ngrd+1;
    pos     = vector<double>(ngrd,  0.0);
    posf    = vector<double>(ngrdf, 0.0);

    uMean   = vector<double>(ngrd,  0.0);
    vMean   = vector<double>(ngrd,  0.0);
    wMean   = vector<double>(ngrd,  0.0);

    uMsqr   = vector<double>(ngrd, 0.0);
    vMsqr   = vector<double>(ngrd, 0.0);
    wMsqr   = vector<double>(ngrd, 0.0);

    vTrans  = vector<double>(ngrd, 0.0);

    double dp = Ldomain/ngrd;
    for(int i=1; i<ngrdf; i++)
        posf[i] = posf[i-1] + dp;
    pos[0] = dp/2;
    for(int i=1; i<ngrd; i++)
        pos[i] = pos[i-1] + dp;
    
#ifdef NEWSTATS
    // edstat(j,k,l,i): i = coordinate
    //                  j = del_U_eddy, del_U_eddy^2, del_U_flow, del_U_flow^2, del_U_all, del_U_all^2, del_U_adpt, del_U_adpt^2
    //                  k = u, v, w, T
    //                  l = average period
    edstat      = vector<vector<vector<vector<double> > > > 
                    (8, vector<vector<vector<double> > >
                    (4, vector<vector<double> >
                    (odtP->nStat, vector<double>(ngrd, 0.0) ) ) );
    // cstat(m,k,l,i):  i = coordinate
    //                  m = u, u^2, du   (eddy, eddy near wall, eddy influence region)
    //                  k = u, v, w, T, eddy
    //                  l = average period
    // ctime(l):        l = average period
    cstat       = vector<vector<vector<vector<double> > > >
                    (3, vector<vector<vector<double> > >
                    (5, vector<vector<double> >
                    (odtP->nStat, vector<double>(ngrd, 0.0) ) ) );
    ctime       = vector<double>(odtP->nStat, 0.0);
    // oldVars(k,i)     i = coordinate
    //                  k = u, v, w
    oldVars     = vector<vector<double> > (4, vector<double>(ngrd, 0.0) );
    
    // currently unused, needed for further work on 2-phase flows (multi-phase flows)
    // phstat(n,l,i)    i = coordinate
    //                  n = phase corresponding to phases(n)
    //                  l = average period
    phstat      = vector<vector<vector<double> > >
                    (1, vector<vector<double> >
                    (odtP->nStat, vector<double>(ngrd, odtP->phase) ) );
    // phases(n)        n = phase
    phases      = vector<double>(1, odtP->phase);
    
    // tau_wstat(o,p)   o = (tau_w, min bound, max bound)
    //                  p = (-4:0.1:6) array for containers of tau_w+
    //                      initial calculation for tau_w+ = 1 is based on domain
    //                      size l and pressure gradient gradP. 
    //                      <tau_w> = l/2 * gradP * tau_w+
    // Derivation: Re_tau^2 = h^3 * gradP / nu^2        h = l/2
    //             Re_tau   = h * u_tau / nu
    //             => u_tau = sqrt( h * gradP )
    //                u_tau = sqrt( <tau_w> )
    tau_wstat   = vector<vector<double> > (4, vector<double>(1000, 0.0));
    mTauWlo = 0.0; mTauWlo2 = 0.0;
    mTauWup = 0.0; mTauWup2 = 0.0;
    double minTauWp = -4; double maxTauWp = 6;
    double delTau = (maxTauWp-minTauWp) / tau_wstat[0].size();
    meanTauW = odtP->domainLength / 2 * abs(odtP->dPdx);
    for (int ii = 0; ii < (int)tau_wstat[0].size(); ii++)
    {
        tau_wstat[2][ii] = meanTauW * (delTau * ii      + minTauWp); // lower bound of container
        tau_wstat[3][ii] = meanTauW * (delTau * (ii +1) + minTauWp); // upper bound of container
    }

    // currently unused, needed for further steps where the Domain can be increased
    // setting default values for ngrd_av, Ldomain_av, and their max values
    ngrd_av = vector<int>(1, ngrd);
    max_ngrd = ngrd;
    Ldomain_av = vector<double>(1, Ldomain);
    max_Ldomain = Ldomain;
#endif

}
///////////////////////////////////////////////////////////////////////////////
/** This constructor is used when the domain size is changed during the
 *  simulation. The new stats are adjusted to the new domain size, while keeping
 *  the old information. After this change, each following period uses the new
 *  stats grid. The unused space of edstat and cstat for previous periods are
 *  filled with zeros and can be ignorred.
 *  Each following realisation now always setts ngrd, ngrdf, Ldomain, pos, posf,
 *  (u,v,w)Mean, (u,v,w)Msqr to the right ngrd for current period. The needed
 *  data for this behavior is saved in ngrd_av, max_ngrd and Ldomain.
 *
 *  @param oldStat \input the old stats object.
 *  @param iStat \input number of the current period
 *  @param newLdomain \input the new domain length after extension
 */
#ifdef NEWSTATS
stats::stats(const stats *oldStats, const int iStat, double newLdomain) {
    
    *proc.ostrm << endl << "#--- calling stats-constructor to adjust to a changed domain size";
    odtP       = oldStats->odtP;
    
    ngrd_av    = vector<int>(odtP->nStat, oldStats->ngrd);
    Ldomain_av = vector<double>(odtP->nStat, oldStats->Ldomain);
    
    
    // FALKO: Muss getestet werden !!
    ngrd       = static_cast<int>(oldStats->ngrd * newLdomain / oldStats->Ldomain) +1;
    //double dely    = oldStats-Ldomain / oldStats->ngrd;
    //int    Nnew    = static_cast<int>(newLdomain / dely);
    //double delynew = newLdomain / Nnew;
    //if (dely == delynew){
    //    ngrd       = Nnew;
    //}
    //else{
    //    ngrd       = Nnew + 1;
    //    *proc.ostrm << endl << "#-- the new domain size is not a multiple of "
    //        << "the old grid distance.\n#-- Therefore the old static grid is "
    //        << "changed for the rest of the simulation.";
    //}
    
    max_ngrd   = max(ngrd,oldStats->max_ngrd);
    max_Ldomain= max(newLdomain,oldStats->max_Ldomain);
    ngrdf      = ngrd+1;
    Ldomain    = newLdomain;
    
    pos     = vector<double>(ngrd,  0.0);
    posf    = vector<double>(ngrdf, 0.0);
    
    double dp = Ldomain/ngrd;
    for(int i=1; i<ngrdf; i++)
        posf[i] = posf[i-1] + dp;
    pos[0] = dp/2;
    for(int i=1; i<ngrd; i++)
        pos[i] = pos[i-1] + dp;
    
    // filling the vectors ngrd_av and Ldomain_av
    if (oldStats->ngrd_av.size() == 1) {
        // ngrd_av of oldStats does not exist
        for(int i=iStat; i<=odtP->nStat; i++){
            // setting all following average period ngrds to this ngrd
            ngrd_av[i-1]    = ngrd;
            Ldomain_av[i-1] = Ldomain; 
        }
    }
    else{
        // ngrd_av of odtStats does exist
        ngrd_av = oldStats->ngrd_av;
        for(int i=iStat; i<=odtP->nStat; i++){
            // setting all following average period ngrds to this ngrd
            ngrd_av[i-1]    = ngrd;
            Ldomain_av[i-1] = Ldomain;
        }
        if (oldStats->ngrd > max_ngrd) max_ngrd = oldStats->ngrd;
    }
    
    uMean   = vector<double>(ngrd,  0.0);
    vMean   = vector<double>(ngrd,  0.0);
    wMean   = vector<double>(ngrd,  0.0);
    
    uMsqr   = vector<double>(ngrd,  0.0);
    vMsqr   = vector<double>(ngrd,  0.0);
    wMsqr   = vector<double>(ngrd,  0.0);
    
    vTrans  = vector<double>(ngrd,  0.0);
    
    edstat  = vector<vector<vector<vector<double> > > >
                (8, vector<vector<vector<double> > >
                (4, vector<vector<double> >
                (odtP->nStat, vector<double>
                (max_ngrd ,0.0) ) ) );
    
    cstat   = vector<vector<vector<vector<double> > > >
                (3, vector<vector<vector<double> > >
                (5, vector<vector<double> >
                (odtP->nStat, vector<double>
                (max_ngrd ,0.0) ) ) );
    
    phases  = oldStats->phases;
    phstat  = vector<vector<vector<double> > >
                (phases.size(), vector<vector<double> >
                (odtP->nStat, vector<double>
                (max_ngrd, 0.0) ) );
    
    // copy curren content of edstat and cstat from oldStat to the new edstat
    // and cstat
    
    int N = 0;
    for(int i=0; i<odtP->nStat; i++) {
        if(oldStats->ngrd_av.size() == 1)
            N = oldStats->ngrd;
        else
            N = oldStats->ngrd_av[i];
        
        //for(int j=0; j<N; j++){
        for(int k=0; k<(int)edstat.size(); k++){ // flow, flow^2, eddy, eddy^2 ...
            for(int l=0; l<(int)edstat[k].size(); l++){ // u, v, w, T
                for(int j=0; j<N; j++){
                    edstat[k][l][i][j] = oldStats->edstat[k][l][i][j];
                }
            }
        }
        for(int k=0; k<(int)cstat.size(); k++){ // u, u^2, du
            for(int l=0; l<(int)cstat[k].size(); l++) { // u, v, w, T, eddy
                for(int j=0; j<N; j++){
                     cstat[k][l][i][j] = oldStats->cstat[k][l][i][j];
                }
            }
        }
        for(int k=0; k<(int)phases.size(); k++){ // each phase
            for(int j=0; j<N; j++){
                phstat[k][i][j] = oldStats->phstat[k][i][j];
            }
        }
    }
    
    ctime   = oldStats->ctime;
    
    oldVars = vector<vector<double> > (4, vector<double>(ngrd ,0.0));
    
}


///////////////////////////////////////////////////////////////////////////////
/** This is a copy-constuctor used during reading a change file. See 
 *  odtSolver::ReadChangeFile for more details.
 *  
 *  @Param newStats The stats-object to copy
 */
void stats::operator=(const stats &newStats){
    ngrd    = newStats.ngrd;
    ngrdf   = ngrd+1;
    Ldomain = newStats.Ldomain;
    odtP    = newStats.odtP;
    pos     = newStats.pos;
    posf    = newStats.posf;
    uMean   = newStats.uMean;
    vMean   = newStats.vMean;
    wMean   = newStats.wMean;
    uMsqr   = newStats.uMsqr;
    vMsqr   = newStats.vMsqr;
    wMsqr   = newStats.wMsqr;
    vTrans  = newStats.vTrans;
    edstat  = newStats.edstat;
    cstat   = newStats.cstat;
    ctime   = newStats.ctime;
    oldVars = newStats.oldVars;
    phstat  = newStats.phstat;
    phases  = newStats.phases;
    ngrd_av = newStats.ngrd_av;
    Ldomain_av  = newStats.Ldomain_av;
    max_ngrd    = newStats.max_ngrd;
    max_Ldomain = newStats.max_Ldomain;
}
#endif


///////////////////////////////////////////////////////////////////////////////
/**Call this once the odtline has been initialized. This is done so that the 
 * stats are not initially biased by the starting value in the constructor.
 * This could probably be moved into the constructor too, if desired.
 *
 * @param odtl \input odtline object to get stats for.
 */
void stats::initStats(odtline *odtl) {


    //odtGrd2statGrd(odtl->posf, odtl->uvel, odtP->pJump[0]);
    odtGrd2statGrd_c(odtl->posf, odtl->pos, odtl->uvel, odtP->pJump[0],
            (int)(odtl->bcprops[0][0]), odtl->bcprops[0][1], odtl->bcprops[0][2], odtl->bcprops[0][3],
            (int)(odtl->bcprops[0][4]), odtl->bcprops[0][5], odtl->bcprops[0][6], odtl->bcprops[0][7]);
    uMean = vTrans;

    //odtGrd2statGrd(odtl->posf, odtl->vvel, odtP->pJump[1]);
    odtGrd2statGrd_c(odtl->posf, odtl->pos, odtl->vvel, odtP->pJump[1],
            (int)(odtl->bcprops[1][0]), odtl->bcprops[1][1], odtl->bcprops[1][2], odtl->bcprops[1][3],
            (int)(odtl->bcprops[1][4]), odtl->bcprops[1][5], odtl->bcprops[1][6], odtl->bcprops[1][7]);
    vMean = vTrans;

    //odtGrd2statGrd(odtl->posf, odtl->wvel, odtP->pJump[2]);
    odtGrd2statGrd_c(odtl->posf, odtl->pos, odtl->wvel, odtP->pJump[2],
            (int)(odtl->bcprops[2][0]), odtl->bcprops[2][1], odtl->bcprops[2][2], odtl->bcprops[2][3],
            (int)(odtl->bcprops[2][4]), odtl->bcprops[2][5], odtl->bcprops[2][6], odtl->bcprops[2][7]);
    wMean = vTrans;

    Ldomain = odtl->Ldomain;

    for(int i=0; i<ngrd; i++) {
        uMsqr[i] = uMean[i]*uMean[i];
        vMsqr[i] = vMean[i]*vMean[i];
        wMsqr[i] = wMean[i]*wMean[i];
    }

}

///////////////////////////////////////////////////////////////////////////////
#ifdef NEWSTATS
/** This function returns the bulk velocity in i direction
 *
 *  @param i \input the coordinate direction (0,1,2) = (u,v,w)
 *  @param uBulk \output the bulk velocity
 */
double stats::getUBulk(int i, int istat){
    double uBulk = 0.0;

    if (istat > 1){
        if (ctime[istat-1] < ctime[istat-2]*0.25){
            for (int j=0; j<ngrd; j++){
                uBulk += cstat[0][i][istat-1][j] + cstat[0][i][istat-2][j];
            }
            uBulk /= ngrd;
            return uBulk/(ctime[istat-1]+ctime[istat-2]);
        }else{
            for (int j=0; j<ngrd; j++){
                uBulk += cstat[0][i][istat-1][j];
            }
        }
    }
    else{
        for(int j=0; j<ngrd; j++){
            uBulk += cstat[0][i][istat-1][j];
        }
    }
    uBulk /= ngrd;
    return uBulk/ctime[istat-1];
}
#endif


///////////////////////////////////////////////////////////////////////////////
/**Compute current values of mean variables.
 * You could also pass in an anyline, but then need to get the right position
 * in the props array to get uvel, etc.
 * Currently this will give a Reynolds average, not a Favre Average
 *
 * @param odtl \input odtline for computing stats.
 * @param t \input current realization time.
 * @param dt \input time step size.
 */
void stats::updateMeans(odtline &odtl, const double &t, const double &dt) {


    //odtGrd2statGrd(odtl.posf, odtl.uvel, odtP->pJump[0]);
    odtGrd2statGrd_c(odtl.posf, odtl.pos, odtl.uvel, odtP->pJump[0],
            (int)(odtl.bcprops[0][0]), odtl.bcprops[0][1], odtl.bcprops[0][2], odtl.bcprops[0][3],
            (int)(odtl.bcprops[0][4]), odtl.bcprops[0][5], odtl.bcprops[0][6], odtl.bcprops[0][7]);
    addToVec( uMean, vTrans, t, dt);
    addToVec2(uMsqr, vTrans, t, dt);

    //odtGrd2statGrd(odtl.posf, odtl.vvel, odtP->pJump[1]);
    odtGrd2statGrd_c(odtl.posf, odtl.pos, odtl.vvel, odtP->pJump[1],
            (int)(odtl.bcprops[1][0]), odtl.bcprops[1][1], odtl.bcprops[1][2], odtl.bcprops[1][3],
            (int)(odtl.bcprops[1][4]), odtl.bcprops[1][5], odtl.bcprops[1][6], odtl.bcprops[1][7]);
    addToVec( vMean, vTrans, t, dt);
    addToVec2(vMsqr, vTrans, t, dt);

    //odtGrd2statGrd(odtl.posf, odtl.wvel, odtP->pJump[2]);
    odtGrd2statGrd_c(odtl.posf, odtl.pos, odtl.wvel, odtP->pJump[2],
            (int)(odtl.bcprops[2][0]), odtl.bcprops[2][1], odtl.bcprops[2][2], odtl.bcprops[2][3],
            (int)(odtl.bcprops[2][4]), odtl.bcprops[2][5], odtl.bcprops[2][6], odtl.bcprops[2][7]);
    addToVec( wMean, vTrans, t, dt);
    addToVec2(wMsqr, vTrans, t, dt);

}


///////////////////////////////////////////////////////////////////////////////
/**Adds an instantaneous profile to a mean profile. Weighting each
 * by their respective ages: the current mean has a weight of the
 * simulation age, while the profile to add gets a weight of dt.
 * \cond
 * <x> = ( x1*dt1 + x2*dt2 + x3*dt3 ) / (t=dt1+dt2+dt3)
 * <x> = ( <x>*t + x4*dt4 ) / (t + dt4)
 *
 * NOTE: BELOW IS JUST LATEX STYLE FORMULAS UNTIL PARAMTERS
 * \endcond
 *
 * \f[
 *      \langle x \rangle = \frac{ x_1 d t_1 + x_2 d t_2 + x_3 d t_3 }{(t = d t_1+d t_2+d t_3)}
 * \f]
 * \f[
 *      \langle x \rangle = \frac{ \langle x \rangle t + x_4 d t_4 }{ t + d t_4}
 * \f]
 *
 *
 * @param vecBase  \inout mean variable to update.
 * @param vecToAdd \input variable to augment the mean.
 * @param t        \input current realization time.
 * @param dt       \input time step size.
 */
void stats::addToVec(std::vector<double> &vecBase,
                     std::vector<double> &vecToAdd,
                     const double &t, const double &dt) {


    for(int i=0; i<(int)vecToAdd.size(); i++)
        vecBase[i] = (vecBase[i] * t + vecToAdd[i]*dt)/(t+dt);

}


///////////////////////////////////////////////////////////////////////////////
/**Like addToVec, but square quantities instead of the value itself.
 * \cond
 * <x> = ( x1*dt1 + x2*dt2 + x3*dt3 ) / (t=dt1+dt2+dt3)
 * <x> = ( <x>*t + x4*dt4 ) / (t + dt4)
 *
 * NOTE: BELOW IS JUST LATEX STYLE FORMULAS UNTIL PARAMTERS
 * \endcond
 *
 * \f[
 *      \langle x \rangle = \frac{ x_1 d t_1 + x_2 d t_2 + x_3 d t_3 }{(t = d t_1+d t_2+d t_3)}
 * \f]
 * \f[
 *      \langle x \rangle = \frac{ \langle x \rangle t + x_4 d t_4 }{ t + d t_4}
 * \f]
 *
 *
 * @param vecBase  \inout mean variable to update.
 * @param vecToAdd \input variable to augment the mean.
 * @param t        \input current realization time.
 * @param dt       \input time step size.
 */
void stats::addToVec2(std::vector<double> &vecBase,
                     std::vector<double> &vecToAdd,
                     const double &t, const double &dt) {


    for(int i=0; i<(int)vecToAdd.size(); i++)
        vecBase[i] = (vecBase[i] * t + vecToAdd[i]*vecToAdd[i]*dt)/(t+dt);

}


///////////////////////////////////////////////////////////////////////////////
/** This function wrapps around the odtposf vector and the odtvec vector for
 *  later use in odtGrd2statGrd().
 *
 *  @param odtposf \input vector of face positions on the odtline.
 *  @param odtvec  \input vector of an odtline variable.
 *  @param pJump   \input periodic jumps.
 */
void stats::wrapAround(std::vector<double> &odtposf,
                        std::vector<double> &odtvec, double &pJump) {
    int odtngrd = odtposf.size()-1;
    
    //---------- If odtposf[0] != 0 (e.g, periodic) then rebuild the odt vectors
    //---------- by wrapping around the domain

    if(odtP->Lperiodic && odtposf[0] != 0.0) {
        
        vector<double> vd;

        int nDoff = static_cast<int> (odtposf[0] / Ldomain);   // # of full domains offset (usually 0)
        if(nDoff != 0) 
            for(int i=0; i<(int)odtposf.size(); i++)
                odtposf[i] -= nDoff * Ldomain;                 // shift domain 

        for(int i=odtngrd-1; i>=0; i--)
            if(odtposf[i] < Ldomain) {
                int ipos = i;                                 // cell that splits hi boundary

                //---------- split cell ipos at the domain boundary

                if(odtposf[ipos+1] != Ldomain) {
                    odtposf.insert(odtposf.begin()+ipos+1, Ldomain);
                    odtvec.insert( odtvec.begin() +ipos+1, odtvec[ipos]);
                    odtngrd++;        
                    ipos++;
                }
                
                //---------- now wrap the posf and variable vecs

                vd = vector<double>(odtvec.begin()+ipos, odtvec.end());
                if(pJump != 0.0) 
                    for(int j=0; j<(int)vd.size(); j++)
                        vd[j] -= pJump;
                odtvec.erase(odtvec.begin()+ipos, odtvec.end());
                odtvec.insert(odtvec.begin(), vd.begin(), vd.end());


                vd = vector<double>(odtposf.begin()+ipos, odtposf.end()-1);
                for(int j=0; j<(int)vd.size(); j++)
                    vd[j] -= Ldomain;
                odtposf.erase(odtposf.begin()+ipos, odtposf.end()-1);
                odtposf.insert(odtposf.begin(), vd.begin(), vd.end());
                odtposf[odtngrd] = Ldomain;
                
                break;
            }
    }
}


///////////////////////////////////////////////////////////////////////////////
// 1st order constant interpolation (used for old statistics)
/**Used to transfer a vector on the odt grid to the stats grid.
 * Routine integrates odt variable odtVec, filling vTrans data member.
 * Loop over each stat grid cell and for each cell march along the odt grid.
 * Consider the grids: 
 * \vc{
 * stat:  |           |           |           |
 * odt:   |   |   |   :     |     |           |
 *                  pLeft
 * }
 * \c pLeft refers to the overlap.  Note overlap on left of second stat cell.
 * When start second stat cell, go from pLeft to fifth stat face.
 * This will work with nonuniform stat grids.
 *
 * @param odtposf \input vector of face positions on the odtline.
 * @param odtvec  \input vector of an odtline variable.
 * @param pJump   \input periodic jumps.
 */
void stats::odtGrd2statGrd(std::vector<double> odtposf, 
                           std::vector<double> odtvec, double pJump) {

    int    odtngrd = odtposf.size()-1;

    //---------- If odtposf[0] != 0 (e.g, periodic) then rebuild the odt vectors
    //---------- by wrapping around the domain
    
    wrapAround(odtposf, odtvec, pJump);

    //---------- Now transfer grids

    int    i, ip;                       // for stat grid
    int    j, jp;                       // for odt grid
    double pLeft;                       // for the overlap

    //---------- transfer grids

    j     = 0;
    jp    = 1;
    pLeft = posf[0];

    //vTrans = vector<double>(ngrd, 0.0);
    std::fill(vTrans.begin(), vTrans.end(), 0.0);

    for(i=0, ip=1; i<ngrd; i++, ip++) {                     // loop over stat grd
        while( odtposf[jp] < posf[ip] && jp!=odtngrd) {     // loop over odt grd
            vTrans[i] += odtvec[j] * (odtposf[jp]-pLeft);
            j++;
            jp++;
            pLeft = odtposf[j];
        }
        vTrans[i] += odtvec[j] * (posf[ip] - pLeft);        // get the leftovers
        pLeft = posf[ip];

    }

    //---------- normalize

    for(i=0, ip=1; i<ngrd; i++, ip++) 
        vTrans[i] /= (posf[ip]-posf[i]);
}


///////////////////////////////////////////////////////////////////////////////
// interpolation using a cubic spline interpolation
//
// upper (BCu) and lower (BC_l) bc's;  nn. = not needed
// periodic-BC      BC = 0, a = nn.,   b = nn.,   c = nn.
// Dirichtet-BC:    BC = 1, a = f(x),  b = nn.,   c = nn.
// Neumann-BC:      BC = 2, a = f'(x), b = nn.,   c = nn.  (outflow BC)
// Cauchy-BC:       BC = 3, a = f(x),  b = f'(x), c = nn.
// Robin-BC:        BC = 4, c = a*f(x) + b*f'(x)
void stats::odtGrd2statGrd_c(std::vector<double> odtposf, std::vector<double> odtpos,
                             std::vector<double> odtvec, double pJump,
                             int BCl, double al, double bl, double cl,
                             int BCu, double au, double bu, double cu) {
    
    int    odtngrd = odtpos.size();
    
    //---------- If odtposf[0] != 0 (e.g, periodic) then rebuild the odt vectors
    //---------- by wrapping around the domain
    
    wrapAround(odtposf, odtvec, pJump);
    
    // to avoid errors set the last boundey
    if(odtposf[odtngrd] != posf[ngrd])
        posf[ngrd] = odtposf[odtngrd];
    //---------- Now transfer grids

    int    i, ip;                       // for stat grid
    int    j;                           // for odt grid
    double pLeft;                       // for the overlap
    double x0, x1, x2, x3;
    double x0S, x1S, x2S, x3S, xS;      // xS = x_shift
    double xpos, pLeftS;
    double dx0, dx1, dx2;
    double y0, y1, y11, y2, y22, y3;
    // y1 = f(x1);  y11 = f'(x1);  y2 = f(x2);  y22 = f'(x2)
    double a, b, c, d;      // parameters of the kubic curve between x1 and x2

    //---------- transfer grids
    j = 0;
    pLeft = posf[j];
    x0    = odtpos[odtngrd-1]-Ldomain;
    x2    = odtpos[j];
    y2    = odtvec[j];
    x3    = odtpos[j+1];
    y3    = odtvec[j+1];
    dx2   = x3 - x2;
    
    //vTrans = vector<double>(ngrd, 0.0);
    std::fill(vTrans.begin(), vTrans.end(), 0.0);
    
    if(BCl == 0){ // periodic
        x1  = odtpos[odtngrd]-Ldomain;
        dx0 = x0 - x1;
        dx1 = x2 - x1;
        y0  = odtvec[odtngrd-1];
        y1  = odtvec[odtngrd];
        y11 = (y2-y1) * dx0 / ( dx1 * (dx0+dx1) )
              -(y0-y1) * dx1 / ( dx0 * (dx0+dx1) );
    }
    else if(BCl == 1){ // Dirichlet
        x1  = pLeft;
        dx0 = 0.0;
        dx1 = x2 - x1;
        y1  = al;
        y11 = (y2-y1) * (dx1+dx2) / ( dx1 * dx2 )
              -(y3-y1) * dx1 / ( dx2 * (dx1+dx2) );
    }
    else if(BCl == 2){ // Neumann (outflow)
        x1  = pLeft;
        dx0 = 0.0;
        dx1 = x2 - x1;
        y11 = al;
        y1  = ( (dx1+dx2) * (dx1+dx2) * (y2-dx1*y11)
              -( dx1 * dx1 * (y3-(dx1+dx2)*y11) ) )
              / ( (dx1+dx2) * (dx1+dx2) - dx1 * dx1 );
    }
    else if(BCl == 3){ // Cauchy
        x1  = pLeft;
        dx0 = 0.0;
        dx1 = x2 - x1;
        y1  = 0.0;
        y11 = 0.0;
        cout << "\nERROR: in stats::odtGrd2statGrd_c\n";
        cout << "The Cauchy-BC is currently not programmed";
        exit(0);
    }
    else if(BCl == 4){ // Robin
        x1  = pLeft;
        y1  = al;
        y11 = bl;
        cout << "\nERROR: in stats::odtGrd2statGrd_c\n";
        cout << "The Robin-BC is currently not programmed";
        exit(0);
    }
    else{
        cout << "\nERROR: in stats::odtGrd2statGrd_c\n";
        cout << "Unknown BC given. BCl = " << BCl;
        exit(0);
    }
    y22 = (y3-y2) * dx1 / ( dx2 * (dx1+dx2) )
          -(y1-y2) * dx2 / ( dx1 * (dx1+dx2) );
    xS  = x1;
    x0S = x0-xS;
    x1S = 0.0;
    x2S = x2-xS;
    x3S = x3-xS;
    
    //// calculating parameters of cubic curve between x1 and x2
    a = 2*(y1-y2)/(x2S*x2S*x2S) + (y11+y22)/(x2S*x2S);
    b = 3*(y2-y1) / (x2S*x2S) - (2*y11+y22) / x2S;
    c = y11;
    d = y1;
    //*proc.ostrm << endl << scientific << setprecision(16)
    //    << setw(25) << a*x1S*x1S*x1S+b*x1S*x1S+c*x1S+d
    //    << setw(25) << y1
    //    << setw(25) << (a*x1S*x1S*x1S+b*x1S*x1S+c*x1S+d-y1)/y1
    //    << setw(25) << a*x2S*x2S*x2S+b*x2S*x2S+c*x2S+d
    //    << setw(25) << y2
    //    << setw(25) << (a*x2S*x2S*x2S+b*x2S*x2S+c*x2S+d-y2)/y2
    //    << setw(25) << 3*a*x1S*x1S+2*b*x1S+c
    //    << setw(25) << y11
    //    << setw(25) << (3*a*x1S*x1S+2*b*x1S+c-y11)/y11
    //    << setw(25) << 3*a*x2S*x2S+2*b*x2S+c
    //    << setw(25) << y22 
    //    << setw(25) << (3*a*x2S*x2S+2*b*x2S+c-y22)/y22;

    // interpolating the adaptive grid to the static grid
    for(i = 0, ip = 1; i<ngrd; i++, ip++) {
        while( x2 < posf[ip] && j!=odtngrd+1) {
            //cout << endl << "Falko: " << scientific << setprecision(16) << setw(25) << x2 << "< " << setw(25) << posf[ip];
            x2S    = x2 - xS;
            pLeftS = pLeft - xS;
            //*proc.ostrm << endl << scientific << setprecision(16)
            //    << setw(25) << pLeft
            //    << setw(25) << a*pLeftS*pLeftS*pLeftS+b*pLeftS*pLeftS+c*pLeftS+d << endl
            //    << setw(25) << (x2+pLeft)/2
            //    << setw(25) << a/8.0*(x2S+pLeftS)*(x2S+pLeftS)*(x2S+pLeftS)+b/4.0*(x2S+pLeftS)*(x2S+pLeftS)+c/2.0*(x2S+pLeftS)+d << endl
            //    << setw(25) << x2
            //    << setw(25) << a*x2S*x2S*x2S+b*x2S*x2S+c*x2S+d;
            vTrans[i] += a * (x2S *x2S *x2S *x2S - pLeftS *pLeftS *pLeftS *pLeftS) / 2.0 / 2.0
                        +b * (x2S *x2S *x2S      - pLeftS *pLeftS *pLeftS)        / 3.0
                        +c * (x2S *x2S           - pLeftS *pLeftS)               / 2.0
                        +d * (x2S                - pLeftS);
            j++;
            pLeft = x2;
            x0  = x1;  x1  = x2; x2 = x3;
            dx0 = dx1; dx1 = dx2;
            y0  = y1;  y1  = y2; y2 = y3;
            y11 = y22;
            if (j == odtngrd-1) { // second last cell
                if(BCu == 0){ // periodic
                    x3  = odtpos[0] + Ldomain;
                    dx2 = x3 - x2;
                    y3  = odtvec[0];
                    y22 = (y3-y2) * dx1 / ( dx2 * (dx1+dx2) )
                            -(y1-y2) * dx2 / ( dx1 * (dx1+dx2) );
                }else if(BCu == 1){ // Dirichlet
                    x3  = odtposf[j+1];
                    dx2 = x3 - x2;
                    y3  = au;
                    y22 = (y3-y2) * dx1 / ( dx2 * (dx1+dx2) )
                            -(y1-y2) * dx2 / ( dx1 * (dx1+dx2) );
                }else if(BCu == 2){ // Neumann
                    x3  = odtposf[j+1];
                    dx2 = x3 - x2;
                    y3  = ((dx2 +dx1) * (dx2 +dx1) * (y2 +dx2 *au)
                            -dx2 * dx2 * (y1 + (dx2 +dx1) *au) )
                            / ((dx2+dx1)*(dx2+dx1)-dx2*dx2);
                    y22 = (y3-y2) * dx1 / ( dx2 * (dx1+dx2) )
                            -(y1-y2) * dx2 / ( dx1 * (dx1+dx2) );
                }else if(BCu == 3){ // Cauchy
                    cout << "\nERROR: in stats::odtGrd2statGrd_c\n";
                    cout << "The Cauchy-BC is currently not programmed";
                    exit(0);
                }else if(BCu == 4){ // Robin
                    cout << "\nERROR: in stats::odtGrd2statGrd_c\n";
                    cout << "The Robin-BC is currently not programmed";
                    exit(0);
                }else{
                    cout << "\nERROR: in stats::odtGrd2statGrd_c\n";
                    cout << "Unknown BC given. BCu = " << BCu;
                    exit(0);
                }
            }
            else if (j == odtngrd) { // last cell
                if(BCu == 0){ // periodic
                    x3  = odtpos[1] + Ldomain;
                    dx2 = x3 - x2;
                    y3  = odtvec[1];
                    y22 = (y3-y2) * dx1 / ( dx2 * (dx1+dx2) )
                            -(y1-y2) * dx2 / ( dx1 * (dx1+dx2) );
                }else if(BCu == 1){ // Dirichlet
                    x3  = 0.0;
                    y3  = 0.0;
                    dx2 = 0.0;
                    y22 = ((y2 -y1) *(dx1 +dx0) *(dx1 *dx0) +(y0 -y2) *dx1 *dx1) 
                          / ( (dx0 +dx1) *dx0 *dx1);
                }else if(BCu == 2){ // Neumann
                    x3  = 0.0;
                    dx2 = 0.0;
                    y3  = 0.0;
                    y22 = au;
                }else if(BCu == 3){ // Cauchy
                    cout << "\nERROR: in stats::odtGrd2statGrd_c\n";
                    cout << "The Cauchy-BC is currently not programmed";
                    exit(0);
                }else if(BCu == 4){ // Robin
                    cout << "\nERROR: in stats::odtGrd2statGrd_c\n";
                    cout << "The Robin-BC is currently not programmed";
                    exit(0);
                }else{
                    cout << "\nERROR: in stats::odtGrd2statGrd_c\n";
                    cout << "Unknown BC given. BCu = " << BCu;
                    exit(0);
                }
            }
            else{
                x3  = odtpos[j+1];
                dx2 = x3 - x2;
                y3  = odtvec[j+1];
                y22 = (y3-y2) * dx1 / ( dx2 * (dx1+dx2) )
                        -(y1-y2) * dx2 / ( dx1 * (dx1+dx2) );
            }
            xS  = x1;
            x0S = x0-xS;
            x1S = 0.0;
            x2S = x2-xS;
            x3S = x3-xS;
            a   = 2*(y1-y2)/(x2S*x2S*x2S) + (y11+y22)/(x2S*x2S);
            b   = 3*(y2-y1) / (x2S*x2S) - (2*y11+y22) / x2S;
            c   = y11;
            d   = y1;
            //*proc.ostrm << endl << scientific << setprecision(16)
            //    << setw(25) << a*x1S*x1S*x1S+b*x1S*x1S+c*x1S+d
            //    << setw(25) << y1
            //    << setw(25) << (a*x1S*x1S*x1S+b*x1S*x1S+c*x1S+d-y1)/y1
            //    << setw(25) << a*x2S*x2S*x2S+b*x2S*x2S+c*x2S+d
            //    << setw(25) << y2
            //    << setw(25) << (a*x2S*x2S*x2S+b*x2S*x2S+c*x2S+d-y2)/y2
            //    << setw(25) << 3*a*x1S*x1S+2*b*x1S+c
            //    << setw(25) << y11
            //    << setw(25) << (3*a*x1S*x1S+2*b*x1S+c-y11)/y11
            //    << setw(25) << 3*a*x2S*x2S+2*b*x2S+c
            //    << setw(25) << y22 
            //    << setw(25) << (3*a*x2S*x2S+2*b*x2S+c-y22)/y22;
        }
        xpos   = posf[ip] -xS;
        pLeftS = pLeft    -xS;
        //*proc.ostrm << endl << scientific << setprecision(16)
        //    << setw(25) << pLeft
        //    << setw(25) << a*pLeftS*pLeftS*pLeftS+b*pLeftS*pLeftS+c*pLeftS+d << endl
        //    << setw(25) << (posf[ip]+pLeft)/2
        //    << setw(25) << a/8.0*(xpos+pLeftS)*(xpos+pLeftS)*(xpos+pLeftS)+b/4.0*(xpos+pLeftS)*(xpos+pLeftS)+c/2.0*(xpos+pLeftS)+d << endl
        //    << setw(25) << posf[ip]
        //    << setw(25) << a*xpos*xpos*xpos+b*xpos*xpos+c*xpos+d;
        vTrans[i] += a * (xpos *xpos *xpos *xpos - pLeftS *pLeftS *pLeftS *pLeftS) /2.0 /2.0
                    +b * (xpos *xpos *xpos       - pLeftS *pLeftS *pLeftS) / 3.0
                    +c * (xpos *xpos             - pLeftS *pLeftS) / 2.0
                    +d * (xpos                   - pLeftS);
        pLeft = posf[ip];
    }
    
    // normalize
    for(i = 0; i < ngrd; i++)
        vTrans[i] /= (posf[i+1]-posf[i]);
    
    //*proc.ostrm << "\n\n odtposf odtpos odtvec\n";
    //for(j = 0; j < odtngrd; j++)
    //    cout << scientific << setprecision(16) << setw(25) << odtposf[j] << " " << setw(25) << odtpos[j] << " " << setw(25) << odtvec[j] << endl;
    //*proc.ostrm << scientific << setprecision(16) << setw(25) << odtposf[odtngrd] << endl;
    //*proc.ostrm << "\n\n posf pos vTrans\n";
    //for(i=0; i < ngrd; i++)
    //    cout << scientific << setprecision(16) << setw(25) << posf[i] << " " << setw(25) << pos[i] << " " << setw(25) << vTrans[i] << endl;
    //*proc.ostrm << scientific << setprecision(16) << setw(25) << odtposf[odtngrd] << endl;
    //exit(0);
}

///////////////////////////////////////////////////////////////////////////////
/** Outputs the properties of the stats.   
 *
 *  @param fname \input output file name.
 */
void stats::outputProperties(std::string fname, const int iStat) {

   if(odtP->Llem) return;
   
   
   ofstream ofile(fname.c_str()); 
   if(!ofile) 
       cout << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
   
   ofile << "# grid points = "   << ngrd;
   ofile << "\n# Domain Size = " << Ldomain;
   if(fabs(posf[ngrd]-posf[0]- Ldomain) > 1.0E-6)
       ofile << "\n# last posf-first posf != Ldomain, last posf = " << posf[ngrd];
   ofile << "\n# 1_pos             "
         <<     "2_posf            "
         <<     "3_uMean           "
         <<     "4_vMean           " 
         <<     "5_wMean           "
         <<     "6_uMsqr           "
         <<     "7_vMsqr           "
         <<     "8_wMsqr           "
         ;
   ofile << scientific;
   ofile << setprecision(10);
#ifndef NEWSTATS
   for(int i=0; i<ngrd; i++)
       ofile << endl 
           << setw(19) << pos[i] 
           << setw(19) << posf[i]
           << setw(19) << uMean[i]
           << setw(19) << vMean[i]
           << setw(19) << wMean[i]
           << setw(19) << uMsqr[i]
           << setw(19) << vMsqr[i]
           << setw(19) << wMsqr[i]
           ;
#else
   for(int i=0; i<ngrd;i++)
       ofile << endl
           << setw(19) << pos[i]
           << setw(19) << posf[i]
           << setw(19) << cstat[0][0][iStat-1][i] / ctime[iStat-1]
           << setw(19) << cstat[0][1][iStat-1][i] / ctime[iStat-1]
           << setw(19) << cstat[0][2][iStat-1][i] / ctime[iStat-1]
           << setw(19) << cstat[1][0][iStat-1][i] / ctime[iStat-1]
           << setw(19) << cstat[1][1][iStat-1][i] / ctime[iStat-1]
           << setw(19) << cstat[1][2][iStat-1][i] / ctime[iStat-1]
           ;
#endif
   ofile.close();
}

///////////////////////////////////////////////////////////////////////////////
/** This function is only called if at some point during a simulation a change
 *  of the domain size has occurred. If this has happened the ngrd_av vector 
 *  was set by the change-constructor (the second one). Now for each period
 *  the current right ngrd has to be chosen to make sure that the number of 
 *  gridpoints for the static grid is correct.
 *
 *  @param iStat \input number of the current period
 */
#ifdef NEWSTATS
void stats::ChangeStatGrid(const int iStat){
    ngrd    = ngrd_av[iStat-1];
    ngrdf   = ngrd + 1;
    Ldomain = Ldomain_av[iStat-1];
    
    double dp = Ldomain/ngrd;
    posf[0] = 0.0;
    for(int i=1; i<ngrdf; i++)
        posf[i] = posf[i-1] + dp;
    pos[0]  = dp/2;
    for(int i=1; i<ngrd; i++)
        pos[i] = pos[i-1] + dp;
    
    vTrans.resize(ngrd);
    std::fill(vTrans.begin(), vTrans.end(), 0.0);
    
    oldVars = vector<vector<double> > (3, vector<double>(ngrd ,0.0));
    
    uMean   = vector<double>(ngrd, 0.0);
    vMean   = vector<double>(ngrd, 0.0);
    wMean   = vector<double>(ngrd, 0.0);
    uMsqr   = vector<double>(ngrd, 0.0);
    vMsqr   = vector<double>(ngrd, 0.0);
    wMsqr   = vector<double>(ngrd, 0.0);
    
}
#endif

///////////////////////////////////////////////////////////////////////////////
#ifdef NEWSTATS
/**
 *  This function initializes a temporary vector of the same size as the odtline
 *  used for the data gathering during the time steps . Afterwards, the gathered
 *  date is converted to the stats-grid and added to the original data gathering.
 *  This is only done for temporal flows to leave out the conversion to a stats 
 *  grid. For spacial flows this wont work and each time step has to be converted
 *
 *  @param odtl \input odtline from whiche the data is gathered
 */
void stats::initBStats(odtline &odtl){
    
    //// initialization of the temporal cstats vektor to gather data for cstat
    //if(!(odtP->Lspatial))
    //    if(!odtP->LheatedChannel)
    //        cstats = vector<vector<vector<double> > >
    //                (cstat.size(), vector<vector<double> >
    //                (cstat[0].size()-2, vector<double>
    //                (odtl.ngrd, 0.0) ) );
    //    else
    //        cstats = vector<vector<vector<double> > >
    //                (cstat.size(), vector<vector<double> >
    //                (cstat[0].size()-1, vector<double>
    //                (odtl.ngrd, 0.0) ) );
    
    // The rest of this function is for further work with multi-phase flows
    //
    // scanning the current odtline for unknown phases
    int marker     = 0;
    int num_phases = phases.size();
    for(int j=0; j<odtl.ngrd; j++){
        marker = 0;
        for(int i=0; i<(int)phases.size(); i++){
            if(phases[i] == odtl.phase[j]){
                marker = 1;
                break;
            }
        }
        if(marker == 0){
            //marker = 0;
            //*proc.ostrm << endl << "Phases old: ";
            //for (int k = 0; k < phases.size(); k++)
            //    *proc.ostrm << phases[k];
            //*proc.ostrm << endl;
            *proc.ostrm << endl << "# new phase found at j = " << j;
            // the phase is unknown
            phases.push_back(odtl.phase[j]);
            //*proc.ostrm << endl << "Phases new: ";
            //for (int k = 0; k < phases.size(); k++)
            //    *proc.ostrm << phases[k];
            //*proc.ostrm << endl;
        }
    }
    
    //// initialisation of the temporal phstats vector to gather data for phstat
    //if(!(odtP->Lspatial))
    //    phstats     = vector<vector<double> >
    //                (phases.size(), vector<double>(odtl.ngrd, 0.0) );
    
    // if new phases are detacted then extending the phstat
    if (num_phases != (int)phases.size()){
        odtP->LmultiPhase = true;
        if(odtP->LconstProp == true){
            *proc.ostrm << endl << endl << "ERROR:" << endl;
            *proc.ostrm << "LconsProp = 1 and multiple phases" << endl;
            *proc.ostrm << "A multi phase simulation needs LconstProp = 0!\n";
            exit(0);
        }
        vector< vector< vector<double> > > new_phstat;
        new_phstat = vector< vector< vector<double> > >
                    (phases.size(), vector< vector<double> >
                    (odtP->nStat, vector<double>(max_ngrd ,0.0)));
        
        for(int i=0; i<num_phases; i++){
            for(int k=0; k<odtP->nStat; k++){
                for(int jj=0; jj<max_ngrd; jj++){
                    new_phstat[i][k][jj] = phstat[i][k][jj];
                }
            }
        }
        phstat = new_phstat;
    }
}


///////////////////////////////////////////////////////////////////////////////
/**
 *  This function switches between a temporal and a spazial simulation.
 *  Depending on the simulation case the odtline is saved in the cstats and 
 *  phstats arrays (temporal case) or the information of the odtline is 
 *  converted to the stats grid and added directly to the cstat and phstat 
 *  arrays (spazial case).
 * 
 * @param odtl  \input the current odtline to be gathered
 * @param tStep \input the time step for the gathering
 * @param iStat \input the current statistical period
 *
 */

void stats::BStats(odtline &odtl, const double tStep, const int iStat){
    bool fast = false;
    //if(odtP->Lspatial){
    if(!fast){
        this->BStatsSpace(odtl, tStep, iStat);
        //*proc.ostrm << endl << "### Spatial Stats";
    }
    //else{
    //    this->BStatsTime(odtl, tStep, iStat);
    //    //*proc.ostrm << endl << "### Temporal Stats;";
    //}
    return;
}

//-----------------------------------------------------------------------------
void stats::BStatsTime(odtline &odtl, const double tStep, const int iStat){
    
    int N = odtl.ngrd;
    double delx, delxp;
    // time t
    ctime[iStat-1] += tStep;
    
    // u velocity
    for(int i=0; i<N; i++){
        cstats[0][0][i] += tStep *     odtl.uvel[i];
        cstats[1][0][i] += tStep * pow(odtl.uvel[i],2);
        cstats[0][1][i] += tStep *     odtl.vvel[i];
        cstats[1][1][i] += tStep * pow(odtl.vvel[i],2);
        cstats[0][2][i] += tStep *     odtl.wvel[i];
        cstats[1][2][i] += tStep * pow(odtl.wvel[i],2);
        if (odtP->LheatedChannel){
            cstats[0][3][i] += tStep *     (*odtl.props[3])[i];
            cstats[1][3][i] += tStep * pow((*odtl.props[3])[i],2);
        }
    }
    
    // derivative along y of u,v,w
    delx  = (odtl.posf[3]-odtl.posf[1])/2.0; // second distance between cell centers
    delxp = (odtl.posf[2]-odtl.posf[0])/2.0; // first distance between cell centers
    // first cell, forward approximation second order
    cstats[2][0][0] += pow( (odtl.uvel[0]-odtl.uvel[2]) * delxp / (delx*(delxp+delx))
                                    +(odtl.uvel[1]-odtl.uvel[0]) * (delxp+delx) / (delxp*delx)
                                    ,2) * tStep;
    cstats[2][1][0] += pow( (odtl.vvel[0]-odtl.vvel[2]) * delxp / (delx*(delxp+delx))
                                    +(odtl.vvel[1]-odtl.vvel[0]) * (delxp+delx) / (delxp*delx)
                                    ,2) * tStep;
    cstats[2][2][0] += pow( (odtl.wvel[0]-odtl.wvel[2]) * delxp / (delx*(delxp+delx))
                                    +(odtl.wvel[1]-odtl.wvel[0]) * (delxp+delx) / (delxp*delx)
                                    ,2) * tStep;
    if (odtP->LheatedChannel){
        cstats[2][3][0] += pow( (odtl.temp[0]-odtl.temp[2]) * delxp / (delx*(delxp+delx))
                                    +(odtl.temp[1]-odtl.temp[0]) * (delxp+delx) / (delxp*delx)
                                    ,2) * tStep;
    }
    
    for(int i=1; i<N-1; i++){
        // internal field, central approximation second order
        delx  = delxp; // copy following cell center distance to current center
        delxp = (odtl.posf[i+2]-odtl.posf[i])/2.0; // calculate next cell center distance
        cstats[2][0][i] += pow( (odtl.uvel[i+1]-odtl.uvel[i]) * delx / (delxp*(delx+delxp))
                                        +(odtl.uvel[i]-odtl.uvel[i-1]) * delxp / (delx*(delx+delxp))
                                        ,2) * tStep;
        cstats[2][1][i] += pow( (odtl.vvel[i+1]-odtl.vvel[i]) * delx / (delxp*(delx+delxp))
                                        +(odtl.vvel[i]-odtl.vvel[i-1]) * delxp / (delx*(delx+delxp))
                                        ,2) * tStep;
        cstats[2][2][i] += pow( (odtl.wvel[i+1]-odtl.wvel[i]) * delx / (delxp*(delx+delxp))
                                        +(odtl.wvel[i]-odtl.wvel[i-1]) * delxp / (delx*(delx+delxp))
                                        ,2) * tStep;
        if (odtP->LheatedChannel){
            cstats[2][3][i] += pow( (odtl.temp[i+1]-odtl.temp[i]) * delx / (delxp*(delx+delxp))
                                +(odtl.temp[1]-odtl.temp[0]) * (delxp+delx) / (delxp*delx) ,2) * tStep;
        }
    }
    // last cell, backward approximation second order
    delx  = delxp; // copy last distance between cell centers
    delxp = (odtl.posf[N-1]-odtl.posf[N-3])/2.0; // calculate second last distance
    cstats[2][0][N-1] += pow( (odtl.uvel[N-3]-odtl.uvel[N-1]) * delx / (delxp*(delx+delxp))
                                    +(odtl.uvel[N-1]-odtl.uvel[N-2]) * (delx+delxp) / (delx*delxp)
                                    ,2) * tStep;
    cstats[2][1][N-1] += pow( (odtl.vvel[N-3]-odtl.vvel[N-1]) * delx / (delxp*(delx+delxp))
                                    +(odtl.vvel[N-1]-odtl.vvel[N-2]) * (delx+delxp) / (delx*delxp)
                                    ,2) * tStep;
    cstats[2][2][N-1] += pow( (odtl.wvel[N-3]-odtl.wvel[N-1]) * delx / (delxp*(delx+delxp))
                                    +(odtl.wvel[N-1]-odtl.wvel[N-2]) * (delx+delxp) / (delx*delxp)
                                    ,2) * tStep;
    if (odtP->LheatedChannel){
        cstats[2][3][N-1] += pow( (odtl.temp[N-3]-odtl.temp[N-1]) * delx / (delxp*(delx+delxp))
                +(odtl.temp[N-1]-odtl.temp[N-2]) * (delx+delxp) / (delx*delxp) ,2) * tStep;
    }
    
    // phstats
    for(int i=0; i<(int)phases.size(); i++){
        for(int j=0; j<N; j++){
            //if (curr_phase[j] == phases[i]){
            if (abs(odtl.phase[j] -phases[i]) < 0.01){
                phstats[i][j] += tStep;}
        }
    }
    
    return;
}

//-----------------------------------------------------------------------------
void stats::BStatsSpace(odtline &odtl, const double tStep, const int iStat){
    
    int    N = 0;
    double delx, delxp;

    // time t
    ctime[iStat-1] += tStep;
    
    /////   u velocity   /////
    odtGrd2statGrd_c(odtl.posf, odtl.pos, odtl.uvel, odtP->pJump[0],
                    (int)(odtl.bcprops[0][0]), odtl.bcprops[0][1], odtl.bcprops[0][2], odtl.bcprops[0][3],
                    (int)(odtl.bcprops[0][4]), odtl.bcprops[0][5], odtl.bcprops[0][6], odtl.bcprops[0][7]);
    N = vTrans.size();
    for(int i=0; i<N; i++){
        cstat[0][0][iStat-1][i] += tStep *     vTrans[i];
        cstat[1][0][iStat-1][i] += tStep * pow(vTrans[i],2);
    }

    // derivative along y of u
    delx  = (posf[3]-posf[1])/2.0; // second distance between cell centers
    delxp = (posf[2]-posf[0])/2.0; // first distance between cell centers
    cstat[2][0][iStat-1][0] += pow( (vTrans[0]-vTrans[2]) * delxp / (delx*(delxp+delx))
                                         +(vTrans[1]-vTrans[0]) * (delxp+delx) / (delxp*delx)
                                         ,2) * tStep;
    for(int i=1; i<N-1; i++){
        // internal field, central approximation second order
        delx  = delxp; // copy following cell center distance to current center
        delxp = (posf[i+2]-posf[i])/2.0; // calculate next cell center distance
        cstat[2][0][iStat-1][i] += pow( (vTrans[i+1]-vTrans[i]) * delx / (delxp*(delx+delxp))
                                             +(vTrans[i]-vTrans[i-1]) * delxp / (delx*(delx+delxp))
                                             ,2) * tStep;
    }
    delx  = delxp; // copy last distance between cell centers
    delxp = (posf[N-1]-posf[N-3])/2.0; // calculate second last distance
    cstat[2][0][iStat-1][N-1] += pow( (vTrans[N-3]-vTrans[N-1]) * delx / (delxp*(delx+delxp))
                                           +(vTrans[N-1]-vTrans[N-2]) * (delx+delxp) / (delx*delxp)
                                           ,2) * tStep;

    /////   v velocity   /////
    odtGrd2statGrd_c(odtl.posf, odtl.pos, odtl.vvel, odtP->pJump[1],
                    (int)(odtl.bcprops[1][0]), odtl.bcprops[1][1], odtl.bcprops[1][2], odtl.bcprops[1][3],
                    (int)(odtl.bcprops[1][4]), odtl.bcprops[1][5], odtl.bcprops[1][6], odtl.bcprops[1][7]);
    N = vTrans.size();
    for(int i=0; i<N; i++){
        cstat[0][1][iStat-1][i] += tStep *     vTrans[i];
        cstat[1][1][iStat-1][i] += tStep * pow(vTrans[i],2);
    }

    // derivative along y of u
    delx  = (posf[3]-posf[1])/2.0; // second distance between cell centers
    delxp = (posf[2]-posf[0])/2.0; // first distance between cell centers
    cstat[2][1][iStat-1][0] += pow( (vTrans[0]-vTrans[2]) * delxp / (delx*(delxp+delx))
                                         +(vTrans[1]-vTrans[0]) * (delxp+delx) / (delxp*delx)
                                         ,2) * tStep;
    for(int i=1; i<N-1; i++){
        // internal field, central approximation second order
        delx  = delxp; // copy following cell center distance to current center
        delxp = (posf[i+2]-posf[i])/2.0; // calculate next cell center distance
        cstat[2][1][iStat-1][i] += pow( (vTrans[i+1]-vTrans[i]) * delx / (delxp*(delx+delxp))
                                             +(vTrans[i]-vTrans[i-1]) * delxp / (delx*(delx+delxp))
                                             ,2) * tStep;
    }
    delx  = delxp; // copy last distance between cell centers
    delxp = (posf[N-1]-posf[N-3])/2.0; // calculate second last distance
    cstat[2][1][iStat-1][N-1] += pow( (vTrans[N-3]-vTrans[N-1]) * delx / (delxp*(delx+delxp))
                                           +(vTrans[N-1]-vTrans[N-2]) * (delx+delxp) / (delx*delxp)
                                           ,2) * tStep;

    /////   w velocity   /////
    odtGrd2statGrd_c(odtl.posf, odtl.pos, odtl.wvel, odtP->pJump[2],
                    (int)(odtl.bcprops[2][0]), odtl.bcprops[2][1], odtl.bcprops[2][2], odtl.bcprops[2][3],
                    (int)(odtl.bcprops[2][4]), odtl.bcprops[2][5], odtl.bcprops[2][6], odtl.bcprops[2][7]);
    N = vTrans.size();
    for(int i=0; i<N; i++){
        cstat[0][2][iStat-1][i] += tStep *     vTrans[i];
        cstat[1][2][iStat-1][i] += tStep * pow(vTrans[i],2);
    }

    // derivative along y of u
    delx  = (posf[3]-posf[1])/2.0; // second distance between cell centers
    delxp = (posf[2]-posf[0])/2.0; // first distance between cell centers
    cstat[2][2][iStat-1][0] += pow( (vTrans[0]-vTrans[2]) * delxp / (delx*(delxp+delx))
                                         +(vTrans[1]-vTrans[0]) * (delxp+delx) / (delxp*delx)
                                         ,2) * tStep;
    for(int i=1; i<N-1; i++){
        // internal field, central approximation second order
        delx  = delxp; // copy following cell center distance to current center
        delxp = (posf[i+2]-posf[i])/2.0; // calculate next cell center distance
        cstat[2][2][iStat-1][i] += pow( (vTrans[i+1]-vTrans[i]) * delx / (delxp*(delx+delxp))
                                             +(vTrans[i]-vTrans[i-1]) * delxp / (delx*(delx+delxp))
                                             ,2) * tStep;
    }
    delx  = delxp; // copy last distance between cell centers
    delxp = (posf[N-1]-posf[N-3])/2.0; // calculate second last distance
    cstat[2][2][iStat-1][N-1] += pow( (vTrans[N-3]-vTrans[N-1]) * delx / (delxp*(delx+delxp))
                                           +(vTrans[N-1]-vTrans[N-2]) * (delx+delxp) / (delx*delxp)
                                           ,2) * tStep;
    
    /////   temp      /////
    if (odtP->LheatedChannel){
        odtGrd2statGrd_c(odtl.posf, odtl.pos, (*odtl.props[3]), odtP->pJump[3],
                         (int)(odtl.bcprops[3][0]), odtl.bcprops[3][1], odtl.bcprops[3][2], odtl.bcprops[3][3],
                         (int)(odtl.bcprops[3][4]), odtl.bcprops[3][5], odtl.bcprops[3][6], odtl.bcprops[3][7]);
        N = vTrans.size();
        for(int i=0; i<N; i++){
            cstat[0][3][iStat-1][i] += tStep *     vTrans[i];
            cstat[1][3][iStat-1][i] += tStep * pow(vTrans[i],2);
        }
        // derivative along y of temp
        delx  = (posf[3]-posf[1])/2.0; // second distance between cell centers
        delxp = (posf[2]-posf[0])/2.0; // first distance between cell centers
        cstat[2][3][iStat-1][0] += pow( (vTrans[0]-vTrans[2]) * delxp / (delx*(delxp+delx))
                                                   +(vTrans[1]-vTrans[0]) * (delxp+delx) / (delxp*delx)
                                                   ,2) * tStep;
        for(int i=1; i<N-1; i++){
            // internal field, central approximation second order
            delx  = delxp; // copy following cell center distance to current center
            delxp = (posf[i+2]-posf[i])/2.0; // calculate next cell center distance
            cstat[2][3][iStat-1][i] += pow( (vTrans[i+1]-vTrans[i]) * delx / (delxp*(delx+delxp))
                                                       +(vTrans[i]-vTrans[i-1]) * delxp / (delx*(delx+delxp))
                                                       ,2) * tStep;
        }
        delx  = delxp; // copy last distance between cell centers
        delxp = (posf[N-1]-posf[N-3])/2.0; // calculate second last distance
        cstat[2][3][iStat-1][N-1] += pow( (vTrans[N-3]-vTrans[N-1]) * delx / (delxp*(delx+delxp))
                                                     +(vTrans[N-1]-vTrans[N-2]) * (delx+delxp) / (delx*delxp)
                                                     ,2) * tStep;
    }
    
    /////   phstats   /////
    vector<double> curr_phase;
    N = odtl.ngrd;
    for(int i=0; i<(int)phases.size(); i++){
        curr_phase = vector<double>(N,0.0);
        for(int j=0; j<N; j++){
            //if (curr_phase[j] == phases[i]){
            if (abs(odtl.phase[j] -phases[i]) < 0.01)
                curr_phase[j] = 1.0;
        }
        odtGrd2statGrd(odtl.posf, curr_phase, 0.0); // needs a low order interpolation
        for(int j=0; j<(int)vTrans.size(); j++){
            phstat[i][iStat-1][j] += vTrans[j]*tStep;}
    }
    
    /////   tauW statistic   /////
    N = odtl.ngrd-1;
    double tauWlo = odtl.molec[0] * odtl.uvel[0] / (odtl.pos[0] * odtl.rho[0]);
    double tauWup = odtl.molec[N] * odtl.uvel[N] / ((odtl.Ldomain -odtl.pos[N]) * odtl.rho[N]);
    mTauWlo += tauWlo * tStep; mTauWlo2 += pow(tauWlo,2.0) * tStep;
    mTauWup += tauWup * tStep; mTauWup2 += pow(tauWup,2.0) * tStep;
    for(int ii = 0; ii < (int)tau_wstat[0].size(); ii++)
    {
        if (tauWlo >= tau_wstat[2][ii] && tauWlo < tau_wstat[3][ii])
            tau_wstat[0][ii] += tStep;
        if (tauWup >= tau_wstat[2][ii] && tauWup < tau_wstat[3][ii])
            tau_wstat[1][ii] += tStep;
    }
}


///////////////////////////////////////////////////////////////////////////////
/**
 *  This function is only used if the simulation is temporal. In this case the 
 *  gathered data is temporally stored in cstats and phstats. These arrays are
 *  now converted to the stats grid and added to the arrays cstat and phstat.
 *
 *  @param odtl  \input the faces of the odtline are needed for odtGrd2statGrd()
 *  @param iStat \input the current statistical period
 *
 */
void stats::BStats2statGrd(odtline &odtl, const int iStat){
    
    if(odtP->Lspatial) return;
    
    int N = odtl.ngrd;
    std::vector<double> temp; // temporal vector for use in odtGrd2statGrd()
    temp = vector<double>(N,0.0);

    // adding cstats to cstat
    // k = 0 => uvel   k = 1 => vvel   k = 2 => wvel   k = 3 => T
    for(int k=0; k<(int)cstats[0].size(); k++){
        for(int i=0; i<N; i++){
            temp[i] = cstats[0][k][i];
        }
        odtGrd2statGrd_c(odtl.posf, odtl.pos, temp, odtP->pJump[k],
                        1, 0.0, 0.0, 0.0, 1, 0.0, 0.0, 0.0);
        for(int i=0; i<(int)vTrans.size(); i++){
            cstat[0][k][iStat-1][i] += vTrans[i];}
    }
    
    // k = 0 => uvel^2   k = 1 => vvel^2   k = 2 => wvel^2   k = 3 => T^2
    for(int k=0; k<(int)cstats[0].size(); k++){
        for(int i=0; i<N; i++){
            temp[i] = cstats[1][k][i];}
        odtGrd2statGrd_c(odtl.posf, odtl.pos, temp, odtP->pJump[k]*odtP->pJump[k],
                1, 0.0, 0.0, 0.0, 1, 0.0, 0.0, 0.0);
        for(int i=0; i<(int)vTrans.size(); i++){
            cstat[0][k][iStat-1][i] += vTrans[i];}
    }

    // k = 0 => del uvel   k = 1 => del vvel   k = 2 => del wvel   k = 3 => del T
    double bcl, bcu;
    double dx10, dx20, dx30, dx21, dx31, dx32;
    double dx01, dx02, dx03, dx12, dx13, dx23;
    dx10 = odtl.pos[0]-odtl.posf[0]; dx01 = odtl.posf[N]-odtl.pos[N-1];
    dx20 = odtl.pos[1]-odtl.posf[0]; dx02 = odtl.posf[N]-odtl.pos[N-2];
    dx30 = odtl.pos[2]-odtl.posf[0]; dx03 = odtl.posf[N]-odtl.pos[N-3];
    dx21 = odtl.pos[1]-odtl.pos[0];  dx12 = odtl.pos[N-1]-odtl.pos[N-2];
    dx31 = odtl.pos[2]-odtl.pos[0];  dx13 = odtl.pos[N-1]-odtl.pos[N-3];
    dx32 = odtl.pos[2]-odtl.pos[1];  dx23 = odtl.pos[N-2]-odtl.pos[N-3];
    for(int k=0; k<(int)cstats[0].size(); k++){
        for(int i=0; i<N; i++){
            temp[i] = cstats[2][k][i];
        }
        bcl = ( (dx20*temp[0]-dx10*temp[1]) *dx30 *dx31
                +(dx10*temp[2]-dx30*temp[0]) *dx20 *dx21 )
                / (dx21*dx31*dx32);
        bcu = ( (dx02*temp[N-1]-dx01*temp[N-2]) *dx03 *dx13
                -(dx03*temp[N-1]-dx01*temp[N-3]) *dx02 *dx12 )
                / (dx12 *dx13 *dx23);
        odtGrd2statGrd_c(odtl.posf, odtl.pos, temp, 0.0,
                1, bcl, 0.0, 0.0, 1, bcu, 0.0, 0.0);
        for(int i=0; i<(int)vTrans.size(); i++){
            cstat[2][k][iStat-1][i] += vTrans[i];
        }
    }
    

    // adding phases to phstat
    for(int k=0; k<(int)phases.size(); k++){
        // phase k
        for(int i=0; i<N; i++){
            temp[i] = phstats[k][i];
        }
        odtGrd2statGrd_c(odtl.posf, odtl.pos, temp, 0.0,
                        2, 0.0, 0.0, 0.0, 2, 0.0, 0.0, 0.0);
        for(int i=0; i<(int)vTrans.size(); i++){
            phstat[k][iStat-1][i] += vTrans[i];
        }
    }
    
    return;
}


///////////////////////////////////////////////////////////////////////////////
/**
 *  BSetOld saves the current odtline converted to the stats grid to be used in 
 *  BChange to calculate the difference generated through an eddy event or a 
 *  diffusion step. The conversion to the stats grid is needed due to the change 
 *  of the odt grid during the process.
 *
 *  @param odtl \input the current odtline
 */
void stats::BSetOld(odtline &odtl){
    
    //odtGrd2statGrd(odtl.posf, odtl.uvel, odtP->pJump[0]);
    odtGrd2statGrd_c(odtl.posf, odtl.pos, odtl.uvel, odtP->pJump[0],
            (int)(odtl.bcprops[0][0]), odtl.bcprops[0][1], odtl.bcprops[0][2], odtl.bcprops[0][3],
            (int)(odtl.bcprops[0][4]), odtl.bcprops[0][5], odtl.bcprops[0][6], odtl.bcprops[0][7]);
    for (int i=0; i<(int)vTrans.size(); i++)
        oldVars[0][i] = vTrans[i];
    
    //odtGrd2statGrd(odtl.posf, odtl.vvel, odtP->pJump[1]);
    odtGrd2statGrd_c(odtl.posf, odtl.pos, odtl.vvel, odtP->pJump[1],
            (int)(odtl.bcprops[1][0]), odtl.bcprops[1][1], odtl.bcprops[1][2], odtl.bcprops[1][3],
            (int)(odtl.bcprops[1][4]), odtl.bcprops[1][5], odtl.bcprops[1][6], odtl.bcprops[1][7]);
    for (int i=0; i<(int)vTrans.size(); i++)
        oldVars[1][i] = vTrans[i];
    
    //odtGrd2statGrd(odtl.posf, odtl.wvel, odtP->pJump[2]);
    odtGrd2statGrd_c(odtl.posf, odtl.pos, odtl.wvel, odtP->pJump[2],
            (int)(odtl.bcprops[2][0]), odtl.bcprops[2][1], odtl.bcprops[2][2], odtl.bcprops[2][3],
            (int)(odtl.bcprops[2][4]), odtl.bcprops[2][5], odtl.bcprops[2][6], odtl.bcprops[2][7]);
    for (int i=0; i<(int)vTrans.size(); i++)
        oldVars[2][i] = vTrans[i];
    
    if (odtP->LheatedChannel){
        //odtGrd2statGrd(odtl.posf, odtl.temp, odtP->pJump[3]);
        odtGrd2statGrd_c(odtl.posf, odtl.pos, odtl.temp, odtP->pJump[3],
                (int)(odtl.bcprops[3][0]), odtl.bcprops[3][1], odtl.bcprops[3][2], odtl.bcprops[3][3],
                (int)(odtl.bcprops[3][4]), odtl.bcprops[3][5], odtl.bcprops[3][6], odtl.bcprops[3][7]);
        for (int i=0; i<(int)vTrans.size(); i++)
            oldVars[3][i] = vTrans[i];
    }
    return;
}


///////////////////////////////////////////////////////////////////////////////
/**
 *  BChange saves and gathers the changes generated through a process (eddy 
 *  event, diffusion, ...) into edstat.
 *
 *  @param odtl  \input the current odtline
 *  @param jj    \input switch: 0: eddy event  2: diffusion  4: all processes
 *  @param iStat \input the current gathering period
 */
void stats::BChange(odtline &odtl, const int &jj, const int &iStat){
    // jj = 0  ==>  differences caused by eddies
    // jj = 2  ==>  differences caused by diffusion
    // jj = 4  ==>  differences caused by eddies and diffusion (all changes) 
    // jj = 6  ==>  differences caused by adaption after diffusion and eddies
    
    //odtGrd2statGrd(odtl.posf, odtl.uvel, odtP->pJump[0]);
    odtGrd2statGrd_c(odtl.posf, odtl.pos, odtl.uvel, odtP->pJump[0],
            (int)(odtl.bcprops[0][0]), odtl.bcprops[0][1], odtl.bcprops[0][2], odtl.bcprops[0][3],
            (int)(odtl.bcprops[0][4]), odtl.bcprops[0][5], odtl.bcprops[0][6], odtl.bcprops[0][7]);
    for (int i=0; i<(int)vTrans.size(); i++){
        edstat[0+jj][0][iStat-1][i] +=     vTrans[i]    -     oldVars[0][i];
        edstat[1+jj][0][iStat-1][i] += pow(vTrans[i],2) - pow(oldVars[0][i], 2);
        // it would be better first to copy the square of uvel into a new vector
        // and than convert this new vector to the static grid. But it is costly
    }
    
    //odtGrd2statGrd(odtl.posf, odtl.vvel, odtP->pJump[1]);
    odtGrd2statGrd_c(odtl.posf, odtl.pos, odtl.vvel, odtP->pJump[1],
            (int)(odtl.bcprops[1][0]), odtl.bcprops[1][1], odtl.bcprops[1][2], odtl.bcprops[1][3],
            (int)(odtl.bcprops[1][4]), odtl.bcprops[1][5], odtl.bcprops[1][6], odtl.bcprops[1][7]);
    for (int i=0; i<(int)vTrans.size(); i++){
        edstat[0+jj][1][iStat-1][i] +=     vTrans[i]    -     oldVars[1][i];
        edstat[1+jj][1][iStat-1][i] += pow(vTrans[i],2) - pow(oldVars[1][i], 2);
    }
    
    //odtGrd2statGrd(odtl.posf, odtl.wvel, odtP->pJump[2]);
    odtGrd2statGrd_c(odtl.posf, odtl.pos, odtl.wvel, odtP->pJump[2],
            (int)(odtl.bcprops[2][0]), odtl.bcprops[2][1], odtl.bcprops[2][2], odtl.bcprops[2][3],
            (int)(odtl.bcprops[2][4]), odtl.bcprops[2][5], odtl.bcprops[2][6], odtl.bcprops[2][7]);
    for (int i=0; i<(int)vTrans.size(); i++){
        edstat[0+jj][2][iStat-1][i] +=     vTrans[i]    -     oldVars[2][i];
        edstat[1+jj][2][iStat-1][i] += pow(vTrans[i],2) - pow(oldVars[2][i], 2);
    }
    
    if (odtP->LheatedChannel){
        //odtGrd2statGrd(odtl.posf, odtl.temp, odtP->pJump[3]);
        odtGrd2statGrd_c(odtl.posf, odtl.pos, odtl.temp, odtP->pJump[3],
                (int)(odtl.bcprops[3][0]), odtl.bcprops[3][1], odtl.bcprops[3][2], odtl.bcprops[3][3],
                (int)(odtl.bcprops[3][4]), odtl.bcprops[3][5], odtl.bcprops[3][6], odtl.bcprops[3][7]);
        for (int i=0; i<(int)vTrans.size(); i++){
            edstat[0+jj][3][iStat-1][i] +=     vTrans[i]    -     oldVars[3][i];
            edstat[1+jj][3][iStat-1][i] += pow(vTrans[i],2) - pow(oldVars[3][i], 2);
        }
    }

    // This function and BSetOld currently do not track the changes of the
    // phases through this events!
    return;
}


///////////////////////////////////////////////////////////////////////////////
/**
 *  The function BSnap uses the gatherd data to calculat the mean values, 
 *  the variances, the tke, the budget terms of tke and further data using the 
 *  stats arrays ctime, cstat, edstat, phstat. After calculation the information 
 *  is stored in a file located in the data folder.
 *
 *  @input myid \input the id of the current processor (only used for mpi runs)
 */
void stats::BSnap(const int myid){
    
    *proc.ostrm << "Start subroutine BSnap; myID = " << myid << endl;
    
    int    N  = 0;
    double dz = 0.0;
    double fluxfac2 = 0.0;
    if (ngrd_av.size() == 1){
        ngrd_av    = vector<int>(odtP->nStat, ngrd);
        Ldomain_av = vector<double>(odtP->nStat, Ldomain);
    }
    
    // deklaration and initialisation of needed arrays
    std::vector< std::vector< std::vector< std::vector<double> > > > eavg;
    std::vector< std::vector< std::vector< std::vector<double> > > > cavg;
    std::vector< std::vector< std::vector<double> > >                phavg;
    
    // eavg(j,k,l,i):   i = coordinate
    //                  j = del_U_eddy, del_U_eddy^2, del_U_flow, del_U_flow^2,
    //                      del_U_all, del_U_all^2, del_U_adpt, del_U_adpt^2
    //                  k = u, v, w, T
    //                  l = average periode
    eavg = vector<vector<vector<vector<double> > > > (8,
           vector<vector<vector<double> > > (4,
           vector<vector<double> > (odtP->nStat,
           vector<double>(max_ngrd ,0.0) ) ) );
    // cavg(m,k,l,i):     i = coordinate
    //                    m = u, u^2, du
    //                    k = u, v, w, T
    //                    l = average periode
    cavg = vector<vector<vector<vector<double> > > > (3,
           vector<vector<vector<double> > > (4,
           vector<vector<double> > (odtP->nStat,
           vector<double>(max_ngrd ,0.0) ) ) );
    // phavg(n,l,i):    i = coordinate
    //                  n = phase (corresponding to phases(n))
    //                  l = average periode
    phavg = vector<vector<vector<double> > > (phases.size(),
            vector<vector<double> > (odtP->nStat,
            vector<double>(max_ngrd ,0.0) ) );
    
    std::vector< std::vector<double> > tke;
    std::vector< std::vector<double> > tv;
    std::vector< std::vector<double> > d;
    std::vector< std::vector<double> > d2;
    std::vector< std::vector<double> > ta;
    std::vector< std::vector<double> > p;
    std::vector< std::vector<double> > bal;
    std::vector< std::vector<double> > delt;
    std::vector< std::vector<double> > dadp;
    std::vector<double> temp;
    
    // tke(i,l):        i = coordinate
    //                  l = average periode
    // tv, d, ta, p have the same as tke
    tke  = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd ,0.0) );
    tv   = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd ,0.0) );
    d    = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd ,0.0) );
    d2   = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd ,0.0) );
    ta   = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd ,0.0) );
    p    = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd ,0.0) );
    bal  = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd ,0.0) );
    delt = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd ,0.0) );
    dadp = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd ,0.0) );
    temp = vector<double> (max_ngrd, 0.0);
    
    //*proc.ostrm << endl << "BSnap :: calculate phavg";
    //*proc.ostrm << endl << "BSnap :: odtP->nStat  = " << odtP->nStat;
    //*proc.ostrm << endl << "BSnap :: ngrd_av.size = " << ngrd_av.size();
    //*proc.ostrm << endl << "BSnap :: phases.size  = " << phases.size();
    *proc.ostrm << endl;
    for(int k = 0; k < odtP->nStat; k++){
        for(int i = 0; i < (int)phases.size(); i++){
            for(int j = 0; j < ngrd_av[k]; j++){
                phavg[i][k][j] = phstat[i][k][j] / ctime[k];
                //*proc.ostrm << "k,j,i = " << k << " "<< j << " " << i << " " << ngrd_av[k] << endl;
            }
        }
    }
    
    *proc.ostrm << "BSnap :: calculate cavg" << endl;
    for(int i = 0; i < odtP->nStat; i++){
        //*proc.ostrm << "----- periode = " << k << " ----- ctime = " << ctime[k] << endl;
        for(int k = 0; k < 3; k++){
            for(int l = 0; l < 4; l++){
                for(int j = 0; j < ngrd_av[i]; j++){
                    cavg[k][l][i][j] = cstat[k][l][i][j] / ctime[i];
                }
            }
        }
    }
    
    *proc.ostrm << "BSnap :: calculate eavg" << endl;
    for(int k = 0; k < odtP->nStat; k++){
        //*proc.ostrm << "----- periode = " << k << " ----- ctime = " << ctime[k] << endl;
        for(int i = 0; i < 8; i++){
            for(int l = 0; l < 4; l++){
                for(int j = 0; j < ngrd_av[k]; j++){
                    eavg[i][l][k][j] = edstat[i][l][k][j] / ctime[k];
                }
            }
        }
    }
    
    // normalization of tauWstats
    mTauWlo /= odtP->tEnd; mTauWlo2 /= odtP->tEnd;
    mTauWup /= odtP->tEnd; mTauWup2 /= odtP->tEnd;
    for(int ii=0; ii<(int)tau_wstat[0].size(); ii++)
    {
        tau_wstat[0][ii] /= odtP->tEnd; // stats at lower boundary
        tau_wstat[1][ii] /= odtP->tEnd; // stats at upper boundary
    }
    
    //--------------------------------------------------------------------------
    // Calculation of the butget terms of the turbulent kinetic energy
    //
    // --> assumptions: statistically stationary problem like channel flow
    //
    // ToDo: The derivations at the boundaries are currently calculated as if 
    //       the grid is a finite differences grid -> 1st point is on the wall.
    //--------------------------------------------------------------------------
    
    // begin calculation of the viscous transport
    for (int k = 0; k < odtP->nStat; k++){
        N        = ngrd_av[k];
        dz       = Ldomain_av[k] / (N*1.0);
        fluxfac2 = odtP->visc_0/odtP->rho_0/(dz*dz);
        
        for (int j = 0; j < N; j++){
            tke[k][j]  = cavg[1][0][k][j] -pow(cavg[0][0][k][j], 2); // <u'^2>
            tke[k][j] += cavg[1][1][k][j] -pow(cavg[0][1][k][j], 2); // <v'^2>
            tke[k][j] += cavg[1][2][k][j] -pow(cavg[0][2][k][j], 2); // <w'^2>
        }
        for (int j = 1; j < N-1; j++){
            tv[k][j] = 0.5 * fluxfac2 * (tke[k][j+1] +tke[k][j-1] -2.0 * tke[k][j]);
        }
        tv[k][0]   = 0.5 * fluxfac2 * (tke[k][0]   -2.0 * tke[k][1]   +tke[k][2]);
        tv[k][N-1] = 0.5 * fluxfac2 * (tke[k][N-1] -2.0 * tke[k][N-2] +tke[k][N-3]);
    }
    *proc.ostrm << "BSnap :: calculated tke and tv" << endl;
    
    // end calculation of the viscous transport
    //----------------------------------
    // start calculation of the dissipation
    
    for (int k = 0; k < odtP->nStat; k++){
        for (int j = 0; j < ngrd_av[k]; j++){
            d[k][j]  = -0.5 * eavg[3][0][k][j] +cavg[0][0][k][j] * eavg[2][0][k][j] +tv[k][j];
            d[k][j] += -0.5 * eavg[3][1][k][j] +cavg[0][1][k][j] * eavg[2][1][k][j];
            d[k][j] += -0.5 * eavg[3][2][k][j] +cavg[0][2][k][j] * eavg[2][2][k][j];
        } 
    }
    *proc.ostrm << "BSnap :: calculated d" << endl;
    
    // end calculation of the dissipation
    //---------------------------------
    // start calculation of the shear production
    
    for (int k = 0; k < odtP->nStat; k++){
        N = ngrd_av[k];
        int i = 0;
        temp[N-1] = 0.5 * eavg[0][i][k][N-1];
        for (int j = N-2; j >= 0; j--){
            temp[j] = temp[j+1] +0.5*(eavg[0][i][k][j]+eavg[0][i][k][j+1]);
        }
        for (int j = 1; j < N-1; j++){
            p[k][j] = -0.5*(cavg[0][i][k][j+1]-cavg[0][i][k][j-1])*temp[j];
        }
        p[k][0]   = 0.0;
        p[k][N-1] = 0.0;
        
        for (int l = 1; l <= 2; l++){
            temp[N-1] = 0.5 * eavg[0][l][k][N-1];
            for (int j = N-2; j >= 0; j--){
                temp[j] = temp[j+1] +0.5*(eavg[0][l][k][j] + eavg[0][l][k][j+1]);
            }
            for (int j = 1; j < N-1; j++){
                p[k][j] += -0.5 * (cavg[0][l][k][j+1] - cavg[0][l][k][j-1]) * temp[j];
            }
        }
    }
    *proc.ostrm << "BSnap :: calculated p" << endl;
    
    // end calculation of the shear production
    //---------------------------------
    // start calculation of the advective transport
    
    for (int k = 0; k < odtP->nStat; k++){
        for (int j = 0; j < ngrd_av[k]; j++){
            ta[k][j]  = 0.5 * eavg[1][0][k][j] -cavg[0][0][k][j] * eavg[0][0][k][j] -p[k][j];
            ta[k][j] += 0.5 * eavg[1][1][k][j] -cavg[0][1][k][j] * eavg[0][1][k][j];
            ta[k][j] += 0.5 * eavg[1][2][k][j] -cavg[0][2][k][j] * eavg[0][2][k][j];
        }
    }
    *proc.ostrm << "BSnap :: calculated ta" << endl;
    
    // end calculation of the advective transport
    //---------------------------------
    // start calculating d2
    
    for (int k = 0; k < odtP->nStat; k++){
        N        = ngrd_av[k];
        dz       = Ldomain_av[k] / (N*1.0);
        fluxfac2 = odtP->visc_0/odtP->rho_0;
        
        d2[k][0]  = cavg[2][0][k][0]-pow(((3*cavg[0][0][k][0]-4*cavg[0][0][k][1]+cavg[0][0][k][2])/(2.0*dz)),2);
        d2[k][0] += cavg[2][1][k][0]-pow(((3*cavg[0][1][k][0]-4*cavg[0][1][k][1]+cavg[0][1][k][2])/(2.0*dz)),2);
        d2[k][0] += cavg[2][2][k][0]-pow(((3*cavg[0][2][k][0]-4*cavg[0][2][k][1]+cavg[0][2][k][2])/(2.0*dz)),2);
        d2[k][0] *= fluxfac2;
        
        for(int j=1; j<N-1; j++){
            d2[k][j]  = cavg[2][0][k][j]-pow(((cavg[0][0][k][j+1]-cavg[0][0][k][j-1])/(dz*2.0)),2);
            d2[k][j] += cavg[2][1][k][j]-pow(((cavg[0][1][k][j+1]-cavg[0][1][k][j-1])/(dz*2.0)),2);
            d2[k][j] += cavg[2][2][k][j]-pow(((cavg[0][2][k][j+1]-cavg[0][2][k][j-1])/(dz*2.0)),2);
            d2[k][j] *= fluxfac2;
        }
        d2[k][N-1]  = cavg[2][0][k][N-1]-pow(((3*cavg[0][0][k][N-1]-4*cavg[0][0][k][N-2]+cavg[0][0][k][N-3])/(dz*2.0)),2);
        d2[k][N-1] += cavg[2][1][k][N-1]-pow(((3*cavg[0][1][k][N-1]-4*cavg[0][1][k][N-2]+cavg[0][1][k][N-3])/(dz*2.0)),2);
        d2[k][N-1] += cavg[2][2][k][N-1]-pow(((3*cavg[0][2][k][N-1]-4*cavg[0][2][k][N-2]+cavg[0][2][k][N-3])/(dz*2.0)),2);
        d2[k][N-1] *= fluxfac2;
    }
    *proc.ostrm << "BSnap :: calculated d2" << endl;
    
    // end calculating d2
    //---------------------------------
    // start calculating bal
    
    for (int k = 0; k < odtP->nStat; k++){
        for (int j = 0; j < ngrd_av[k]; j++){
            bal[k][j] = tv[k][j] -d[k][j] +p[k][j] +ta[k][j];
        }
    }
    *proc.ostrm << "BSnap :: calculated bal" << endl;
    
    // end calculating bal
    //---------------------------------
    // start calculating delt

    // delt is the left hand side term d(tke)/dt
    for (int k = 0; k < odtP->nStat; k++){
        for (int j = 0; j < ngrd_av[k]; j++){
            delt[k][j]  = 0.5 * eavg[5][0][k][j] +cavg[0][0][k][j] * eavg[4][0][k][j];
            delt[k][j] += 0.5 * eavg[5][1][k][j] +cavg[0][1][k][j] * eavg[4][1][k][j];
            delt[k][j] += 0.5 * eavg[5][2][k][j] +cavg[0][2][k][j] * eavg[4][2][k][j];
        }
    }
    *proc.ostrm << "BSnap :: calculated delt" << endl;
    
    // end calculating delt
    //---------------------------------
    // start calculating dadp
    
    // dadp is the left hand side term d(tke)/dt
    for (int k = 0; k < odtP->nStat; k++){
        for (int j = 0; j < ngrd_av[k]; j++){
            dadp[k][j]  = 0.5 * eavg[7][0][k][j] +cavg[0][0][k][j] * eavg[6][0][k][j];
            dadp[k][j] += 0.5 * eavg[7][1][k][j] +cavg[0][1][k][j] * eavg[6][1][k][j];
            dadp[k][j] += 0.5 * eavg[7][2][k][j] +cavg[0][2][k][j] * eavg[6][2][k][j];
        }
    }
    *proc.ostrm << "BSnap :: calculated delt" << endl;
    
    // end calculating dadp
    //---------------------------------
    // start calculating budgets for T
    
    std::vector< std::vector<double> > tkeT;
    std::vector< std::vector<double> > dT;
    std::vector< std::vector<double> > d2T;
    std::vector< std::vector<double> > tvT;
    std::vector< std::vector<double> > taT;
    std::vector< std::vector<double> > pT;
    std::vector< std::vector<double> > balT;
    std::vector< std::vector<double> > deltT; 
    std::vector< std::vector<double> > dadpT;
    tkeT   = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd, 0.0) );
    dT     = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd, 0.0) );
    d2T    = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd, 0.0) );
    tvT    = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd, 0.0) );
    taT    = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd, 0.0) );
    pT     = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd, 0.0) );
    balT   = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd, 0.0) );
    deltT  = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd, 0.0) );
    dadpT  = vector<vector<double> > (odtP->nStat, vector<double>(max_ngrd, 0.0) );
    if (odtP->LheatedChannel){
        for (int k = 0; k < odtP->nStat; k++){
            N  = ngrd_av[k];
            dz = Ldomain_av[k] / (N*1.0);
            double fluxfac2T = odtP->lambda_0/(dz*dz);
            
            for (int j = 0; j < N; j++)
                tkeT[k][j] =  cavg[1][3][k][j] -pow(cavg[0][3][k][j], 2); // <T'^2>
            for (int j = 1; j < N-1; j++)
                tvT[k][j] = 0.5 * fluxfac2T * (tkeT[k][j+1] +tkeT[k][j-1] -2.0 * tkeT[k][j]);
            tvT[k][0]   = 0.5 * fluxfac2T * (tkeT[k][0]   -2.0 * tkeT[k][1]   +tkeT[k][2]);
            tvT[k][N-1] = 0.5 * fluxfac2T * (tkeT[k][N-1] -2.0 * tkeT[k][N-2] +tkeT[k][N-3]);
            
            for (int j = 0; j< ngrd_av[k]; j++)
                dT[k][j] = -0.5 * eavg[3][3][k][j] +cavg[0][3][k][j] * eavg[2][3][k][j] +tvT[k][j];
            
            temp[N-1] = 0.5 * eavg[0][3][k][N-1];
            for (int j = N-2; j >= 0; j--)
                temp[j] = temp[j+1] +0.5*(eavg[0][3][k][j]+eavg[0][3][k][j+1]);
            for (int j = 1; j < N-1; j++)
                pT[k][j] = -0.5*(cavg[0][3][k][j+1]-cavg[0][3][k][j-1])*temp[j];
            pT[k][0]   = 0.0;
            pT[k][N-1] = 0.0;
            
            for (int j = 0; j < ngrd_av[k]; j++)
                taT[k][j]  = 0.5 * eavg[1][3][k][j] -cavg[0][3][k][j] * eavg[0][3][k][j] -pT[k][j];
            
            d2T[k][0] = cavg[2][3][k][0]-pow(((3*cavg[0][3][k][0]-4*cavg[0][3][k][1]+cavg[0][3][k][2])/(2.0*dz)),2);
            for (int j = 1; j < N-1; j++)
                d2T[k][j]  = cavg[2][3][k][j]-pow(((cavg[0][3][k][j+1]-cavg[0][3][k][j-1])/(dz*2.0)),2);
            d2T[k][N-1]  = cavg[2][3][k][N-1]-pow(((3*cavg[0][3][k][N-1]-4*cavg[0][3][k][N-2]+cavg[0][3][k][N-3])/(dz*2.0)),2);
            for (int j = 0; j < N-1; j++)
                d2T[k][j] *= fluxfac2T;
            
            for (int j = 0; j < N-1; j++)
                balT[k][j] = tvT[k][j] -dT[k][j] +pT[k][j] +taT[k][j];
            
            for (int j = 0; j < ngrd_av[k]; j++)
                deltT[k][j]  = 0.5 * eavg[5][3][k][j] +cavg[0][3][k][j] * eavg[4][3][k][j];
            
            for (int j = 0; j < ngrd_av[k]; j++)
                dadpT[k][j]  = 0.5 * eavg[7][3][k][j] +cavg[0][3][k][j] * eavg[6][3][k][j];
        }
    }
    
    // end calculating budgets for T
    //---------------------------------
    // write data file for tauW
    
    ofstream file_tauW;
    stringstream out; string outs;
    out.clear(); outs.clear();
    out << proc.dataDir << "tauWstats.dat";
    out >> outs;
    file_tauW.open(outs.c_str());
    double temp2 = 0.0;
    double temp3 = 0.0;
    
    file_tauW << scientific << setprecision(16);
    file_tauW << "# meanTauW  = " << setw(25) << meanTauW << endl;
    file_tauW << "# mTauWlo   = " << setw(25) << mTauWlo  << endl;
    file_tauW << "# mTauWup   = " << setw(25) << mTauWup  << endl;
    file_tauW << "# tauWrmyLo = " << setw(25) << sqrt( mTauWlo2 - pow(mTauWlo,2.0) ) << endl;
    file_tauW << "# tauWrmyUp = " << setw(25) << sqrt( mTauWup2 - pow(mTauWup,2.0) ) << endl;
    file_tauW << "# tau_w                     pdf_lo                   pdf_up" << endl;
    temp2 = tau_wstat[3][0]-tau_wstat[2][0];
    for(int ii=0; ii<(int)tau_wstat[0].size(); ii++)
    {
        temp3 += tau_wstat[0][ii]/temp2;
        file_tauW << setw(25) << (tau_wstat[2][ii]+tau_wstat[3][ii])/2.0
                  << setw(25) << tau_wstat[0][ii]/temp2
                  << setw(25) << tau_wstat[1][ii]/temp2 << endl;
    }
    //file_tauW << "# summ of pdf = " << temp3;
    file_tauW.close();
    
    // end write data file for tauW
    //---------------------------------
    // start data output
    
    *proc.ostrm << "BSnap :: start writing output" << endl;
    // -- writing ASCII output
    ofstream file_avg;
    ofstream file_stat;
    temp2 = 0.0;
    temp3 = 0.0;
    
    for (int k = 0; k < odtP->nStat; k++){
        N = ngrd_av[k];
        stringstream out; string outs;
        out.clear(); outs.clear();
        out << proc.dataDir << "Data-" << myid << "-average-output_" << k << ".dat";
        out >> outs;
        file_avg.open(outs.c_str());
        file_avg << "# grid points = " << N;
        file_avg << "\n# Domain Size = " << odtP->domainLength;
        // could be used for single values to be in the output file
        // the values to be written in the order of the previous line
        file_avg << "\n#          l_min          l_ave"
                 <<    "          l_max            rho"
                 <<    "          molec           visc"
                 <<    "         lambda       lamdaRho"
                 <<    "              t";
        file_avg << setprecision(8);
        file_avg << "\n# " << setw(14) << odtP->Lmin * odtP->nsgrd
                 <<    " " << setw(14) << odtP->Lp   * odtP->nsgrd
                 <<    " " << setw(14) << odtP->Lmax * odtP->nsgrd
                 <<    " " << setw(14) << odtP->rho_0
                 <<    " " << setw(14) << odtP->visc_0
                 <<    " " << setw(14) << odtP->visc_0/odtP->rho_0
                 <<    " " << setw(14) << odtP->lambda_0
                 <<    " " << setw(14) << odtP->lambda_0/odtP->rho_0
                 <<    " " << setw(14) << odtP->tEnd * (k+1) / odtP->nStat;
        file_avg <<  "\n#  y                      "
                 <<     "  u                      "
                 <<     "  v                      "
                 <<     "  w                      "
                 <<     "  u2                     "
                 <<     "  v2                     "
                 <<     "  w2                     "
                 <<     "  eavg                   " // mean(u'v')
                 <<     "  tke                    "
                 <<     "  d                      "
                 <<     "  tv                     "
                 <<     "  p                      "
                 <<     "  ta                     "
                 <<     "  bal                    "
                 <<     "  d2                     "
                 <<     "  delt                   "
                 <<     "  dadp                   ";
        if (odtP->LheatedChannel) {
            file_avg << "  T                      "
                     << "  T2                     "
                     << "  eavgT                  " // mean(T'v')
                     << "  tkeT                   "
                     << "  dT                     "
                     << "  tvT                    "
                     << "  pT                     "
                     << "  taT                    "
                     << "  balT                   "
                     << "  deltT                  "
                     << "  dadpT                  ";
        }
        file_avg <<     "  eddy1                  "
                 <<     "  eddy2                  "
                 <<     "  eddyInfluenceRegion    ";
        for(int i=0; i<(int)phases.size(); i++){
            file_avg <<"  phase_" << (int)phases[i] << "               ";
        }
        file_avg << scientific;
        file_avg << setprecision(16);
        temp2 = 0.0;
        temp3 = 0.0;
        for (int j = 0; j < N; j++){
            temp2 += dz * eavg[0][0][k][j];
            file_avg << endl;
            file_avg << setw(25) << pos[j]
                     << setw(25) << cavg[0][0][k][j]
                     << setw(25) << cavg[0][1][k][j]
                     << setw(25) << cavg[0][2][k][j]
                     << setw(25) << ( cavg[1][0][k][j] -pow(cavg[0][0][k][j],2) )
                     << setw(25) << ( cavg[1][1][k][j] -pow(cavg[0][1][k][j],2) )
                     << setw(25) << ( cavg[1][2][k][j] -pow(cavg[0][2][k][j],2) )
                     << setw(25) << temp2
                     << setw(25) << tke[k][j]
                     << setw(25) << d[k][j]
                     << setw(25) << tv[k][j]
                     << setw(25) << p[k][j]
                     << setw(25) << ta[k][j]
                     << setw(25) << bal[k][j]
                     << setw(25) << d2[k][j]
                     << setw(25) << delt[k][j]
                     << setw(25) << dadp[k][j];
            if (odtP->LheatedChannel) {
                temp3 += dz * eavg[0][3][k][j];
                file_avg << setw(25) << cavg[0][3][k][j]
                         << setw(25) << (cavg[1][3][k][j] -pow(cavg[0][3][k][j],2))
                         << setw(25) << temp3
                         << setw(25) << tkeT[k][j]
                         << setw(25) << dT[k][j]
                         << setw(25) << tvT[k][j]
                         << setw(25) << pT[k][j]
                         << setw(25) << taT[k][j]
                         << setw(25) << balT[k][j]
                         << setw(25) << deltT[k][j]
                         << setw(25) << dadpT[k][j];
            }
            file_avg << setw(25) << cstat[0][4][k][j]
                     << setw(25) << cstat[1][4][k][j]
                     << setw(25) << cstat[2][4][k][j];
            for(int i=0; i<(int)phases.size(); i++){
                file_avg << setw(25) << phavg[i][k][j];
            }
        }
        file_avg.close();
        
        out.clear(); outs.clear();
        out << proc.dataDir << "Data-" << myid << "-stat-output_" << k << ".dat";
        out >> outs;
        file_stat.open(outs.c_str());
        file_stat << "# grid points = " << N;
        file_stat << "\n# Domain Size = " << odtP->domainLength;
        file_stat << "\n#    y                      "
                  <<   "  ctime                  "
                  <<   "  cstat[u]               "
                  <<   "  cstat[v]               "
                  <<   "  cstat[w]               "
                  <<   "  cstat[T]               "
                  <<   "  cstat[eddy1]           "
                  <<   "  cstat[u2]              "
                  <<   "  cstat[v2]              "
                  <<   "  cstat[w2]              "
                  <<   "  cstat[T2]              "
                  <<   "  cstat[eddy2]           "
                  <<   "  cstat[du]              "
                  <<   "  cstat[dv]              "
                  <<   "  cstat[dw]              "
                  <<   "  cstat[dT]              "
                  <<   "  cstat[eddyInfluRegio]  "
                  <<   "  edstat[del_eddy,u]     "
                  <<   "  edstat[del_eddy^2,u]   "
                  <<   "  edstat[del_flow,u]     "
                  <<   "  edstat[del_flow^2,u]   "
                  <<   "  edstat[del_all,u]      "
                  <<   "  edstat[del_all^2,u]    "
                  <<   "  edstat[del_adpt,u]     "
                  <<   "  edstat[del_adpt^2,u]   "
                  <<   "  edstat[del_eddy,v]     "
                  <<   "  edstat[del_eddy^2,v]   "
                  <<   "  edstat[del_flow,v]     "
                  <<   "  edstat[del_flow^2,v]   "
                  <<   "  edstat[del_all,v]      "
                  <<   "  edstat[del_all^2,v]    "
                  <<   "  edstat[del_adpt,v]     "
                  <<   "  edstat[del_adpt^2,v]   "
                  <<   "  edstat[del_eddy,w]     "
                  <<   "  edstat[del_eddy^2,w]   "
                  <<   "  edstat[del_flow,w]     "
                  <<   "  edstat[del_flow^2,w]   "
                  <<   "  edstat[del_all,w]      "
                  <<   "  edstat[del_all^2,w]    "
                  <<   "  edstat[del_adpt,w]     "
                  <<   "  edstat[del_adpt^2,w]   "
                  <<   "  edstat[del_eddy,T]     "
                  <<   "  edstat[del_eddy^2,T]   "
                  <<   "  edstat[del_flow,T]     "
                  <<   "  edstat[del_flow^2,T]   "
                  <<   "  edstat[del_all,T]      "
                  <<   "  edstat[del_all^2,T]    "
                  <<   "  edstat[del_adpt,T]     "
                  <<   "  edstat[del_adpt^2,T]   ";
        file_stat << scientific;
        file_stat << setprecision(16);
        temp2 = 0.0;
        for (int j = 0; j < N; j++){
            temp2 += dz * eavg[0][0][k][j];
            file_stat << endl;
            file_stat << setw(25) << (j*dz)
                      << setw(25) << ctime[k];
            for (int m = 0; m < 3; m++){ // u, u^2, du
                for (int i = 0; i < 5; i++){ // u, v, w, T, eddy
                    file_stat << setw(25) << cstat[m][i][k][j];
                }
            }
            for (int i = 0; i < 4; i++) { // u, v, w, T
                for (int l = 0; l < 8; l++) { // eddy, eddy^2, flow, flow^2, all, all^2, adpt, adpt^2
                    file_stat << setw(25) << edstat[l][i][k][j];
                }
            }
        }
        file_stat.close();
        
    }
    *proc.ostrm << "BSnap :: output in ASCII files finished" << endl;
    
// Will be included later
//#ifdef NETCDF
//    N = max_ngrd;
//    
//    *proc.ostrm << "BSnap :: start writing NetCDF file" << endl;
//    // -- creating new NetCDF-File replacing the one that my be existing
//    stringstream out; string outs;
//    out.clear(); outs.clear();
//    out << proc.dataDir << "Data-" << myid << "-average-output" << ".nc";
//    out >> outs;
//    NcFile file_avgNC(outs.c_str(), NcFile::Replace);
//    *proc.ostrm << "      :: after opening file" << endl;
//    if (!file_avgNC.is_valid()){
//        *proc.ostrm << "ERROR:\nNetCDF-File is not created!" << endl;
//        return;
//    }
//    *proc.ostrm << "      :: after checking file" << endl;
//    // Adding information to file
//    
//    // Adding dimensions to file
//    NcDim* xDim   = file_avgNC.add_dim("y", N); // the location
//    *proc.ostrm << "      :: after dim y" << endl;
//    NcDim* perDim = file_avgNC.add_dim("period", odtP->nStat); // the period 
//    *proc.ostrm << "      :: after dim period" << endl;
//    
//    // Adding variables to file
//    // Declaration
//    NcVar *xNC    = file_avgNC.add_var("y",    ncDouble, xDim);
//    *proc.ostrm << "      :: after var y" << endl;
//    NcVar *timeNC = file_avgNC.add_var("time", ncDouble, perDim);
//    *proc.ostrm << "      :: after var time" << endl;
//    NcVar *uNC    = file_avgNC.add_var("u",    ncDouble, xDim, perDim);
//    *proc.ostrm << "      :: after var u" << endl;
//    NcVar *vNC    = file_avgNC.add_var("v",    ncDouble, xDim, perDim);
//    *proc.ostrm << "      :: after var v" << endl;
//    NcVar *wNC    = file_avgNC.add_var("w",    ncDouble, xDim, perDim);
//    *proc.ostrm << "      :: after var w" << endl;
//    NcVar *u2NC   = file_avgNC.add_var("u2",   ncDouble, xDim, perDim);
//    *proc.ostrm << "      :: after var u2" << endl;
//    NcVar *v2NC   = file_avgNC.add_var("v2",   ncDouble, xDim, perDim);
//    *proc.ostrm << "      :: after var v2" << endl;
//    NcVar *w2NC   = file_avgNC.add_var("w2",   ncDouble, xDim, perDim);
//    *proc.ostrm << "      :: after var w2" << endl;
//    NcVar *eavgNC = file_avgNC.add_var("eavg", ncDouble, xDim, perDim);
//    *proc.ostrm << "      :: after var eavg" << endl;
//    NcVar *pNC    = file_avgNC.add_var("p",    ncDouble, xDim, perDim);
//    *proc.ostrm << "      :: after var p" << endl;
//    NcVar *taNC   = file_avgNC.add_var("ta",   ncDouble, xDim, perDim);
//    *proc.ostrm << "      :: after var ta" << endl;
//    NcVar *tvNC   = file_avgNC.add_var("tv",   ncDouble, xDim, perDim);
//    *proc.ostrm << "      :: after var tv" << endl;
//    NcVar *dNC    = file_avgNC.add_var("d",    ncDouble, xDim, perDim);
//    *proc.ostrm << "      :: after var d" << endl;
//    NcVar *balNC  = file_avgNC.add_var("bal",  ncDouble, xDim, perDim);
//    *proc.ostrm << "      :: after var bal" << endl;
//    NcVar *tkeNC  = file_avgNC.add_var("tke",  ncDouble, xDim, perDim);
//    *proc.ostrm << "      :: after var tke" << endl;
//    NcVar *d2NC  = file_avgNC.add_var("d2",  ncDouble, xDim, perDim);
//    *proc.ostrm << "      :: after var d2" << endl;
//    NcVar *deltNC  = file_avgNC.add_var("delt",  ncDouble, xDim, perDim);
//    *proc.ostrm << "      :: after var delt" << endl;
//    
//    // Filling
//    pNC->put(   &p[0][0],   N, odtP->nStat);
//    *proc.ostrm << "      :: after put p" << endl;
//    taNC->put(  &ta[0][0],  N, odtP->nStat);
//    *proc.ostrm << "      :: after put ta" << endl;
//    tvNC->put(  &tv[0][0],  N, odtP->nStat);
//    *proc.ostrm << "      :: after put tv" << endl;
//    dNC->put(   &d[0][0],   N, odtP->nStat);
//    *proc.ostrm << "      :: after put d" << endl;
//    tkeNC->put( &tke[0][0], N, odtP->nStat);
//    *proc.ostrm << "      :: after put tke" << endl;
//    d2NC->put( &d2[0][0], N, odtP->nStat);
//    *proc.ostrm << "      :: after put d2" << endl;
//    deltNC->put( &delt[0][0], N, odtP->nStat);
//    *proc.ostrm << "      :: after put delt" << endl;
//    balNC->put(&bal[0][0], N, odtP->nStat);
//    *proc.ostrm << "      :: after put bal" << endl;
//    
//    // using p as dummy variable to put the data to the file
//    // writing u2
//    for (int k = 0; k < odtP->nStat; k++){
//        for (int j = 0; j < N; j++){
//            p[j][k] = cavg[j][4][k];
//        }
//    }
//    u2NC->put(&p[0][0], N, odtP->nStat);
//    *proc.ostrm << "      :: after put u2" << endl;
//    
//    // writing v2
//    for (int k = 0; k < odtP->nStat; k++){
//        for (int j = 0; j < N; j++){
//            p[j][k] = cavg[j][5][k];
//        }
//    }
//    v2NC->put(&p[0][0], N, odtP->nStat);
//    *proc.ostrm << "      :: after put v2" << endl;
//    
//    // writing w2
//    for (int k = 0; k < odtP->nStat; k++){
//        for (int j = 0; j < N; j++){
//            p[j][k] = cavg[j][6][k];
//        }
//    }
//    w2NC->put(&p[0][0], N, odtP->nStat);
//    *proc.ostrm << "      :: after put w2" << endl;
//    
//    // writing v
//    for (int k = 0; k < odtP->nStat; k++){
//        for (int j = 0; j < N; j++){
//            p[j][k] = cavg[j][2][k];
//        }
//    }
//    vNC->put(&p[0][0], N, odtP->nStat);
//    *proc.ostrm << "      :: after put v" << endl;
//    
//    // writing w
//    for (int k = 0; k < odtP->nStat; k++){
//        for (int j = 0; j < N; j++){
//            p[j][k] = cavg[j][3][k];
//        }
//    }
//    wNC->put(&p[0][0], N, odtP->nStat);
//    *proc.ostrm << "      :: after put w" << endl;
//    
//    // writing u
//    for (int k = 0; k < odtP->nStat; k++){
//        for (int j = 0; j < N; j++){
//            p[j][k] = cavg[j][1][k];
//        }
//    }
//    uNC->put(&p[0][0], N, odtP->nStat);
//    *proc.ostrm << "      :: after put u" << endl;
//    
//    // writing eavg
//    for (int k = 0; k < odtP->nStat; k++){
//        p[N-1][k] = dz * eavg[N-1][0][0][k];
//        for (int j = N-2; j >= 0; j--){
//            p[j][k] = p[j+1][k] +dz * eavg[j][0][0][k];
//        }
//    }
//    eavgNC->put(&p[0][0], N, odtP->nStat);
//    *proc.ostrm << "      :: after put eavg" << endl;
//    *proc.ostrm << "BSnap :: output in NetCDF files finished" << endl;
//#endif
    
    return;
}


///////////////////////////////////////////////////////////////////////////////
/**
 *  BLoadStats is only used during a mpi run. After calculating BSnap for each 
 *  processor, this function loads the written stats files and add them to the
 *  stats of processor zero. Afterwards, BSnap is called again for the over
 *  all CPUs gathered data.
 * 
 *  @Param myid \input the id of the processor whos data should be read
 */

// ToDo: - including phstat
void stats::BLoadStats(const int myid){
    
    int          N = odtP->nsgrd; 
    ifstream     ifile;
    string       input_s;
    stringstream input_ss;
    double temp2 = 0.0;
    
    //std::vector<double> temp;
    //temp = vector<double>(36 ,0.0);
    
    for (int k = 0; k < odtP->nStat; k++){
        
        input_ss.clear(); input_s.clear();
        input_ss << proc.dataDir << "../data_" << myid << "/Data-" << myid << "-stat-output_" << k << ".dat";
        input_ss >> input_s;
        ifile.open(input_s.c_str(), ifstream::in);
        if (!ifile.is_open()){
            *proc.ostrm << endl << "ERROR: File not opened!" << endl;
            return;
        }
        
        getline(ifile, input_s); // read line "# grid point = 100"
        getline(ifile, input_s); // read line "# Domain Size = 2"
        getline(ifile, input_s); // read line "#  1_pos                    2_u ..."
        
        for (int j = 0; j <= N-1; j++){
            //getline(ifile, input_s);
            ifile >> scientific;
            ifile >> setprecision(16);
            //*proc.ostrm << j << "  ";
            ifile >> temp2; // location y
            //*proc.ostrm << temp2 << "  ";
            ifile >> temp2; // ctime
            //*proc.ostrm << temp2 << "  ";
            if (j == 0){    // only once
                ctime[k] += temp2;
            }
            for (int m = 0; m < 3; m++) {
                for (int i = 0; i < 5; i++) {
                    ifile >> temp2;
                    //*proc.ostrm << temp2 << "  ";
                    cstat[m][i][k][j] += temp2;
                }
            }
            for (int i = 0; i < 4; i++) {
                for (int l = 0; l < 8; l++) {
                    ifile >> temp2;
                    //*proc.ostrm << temp2 << "  ";
                    edstat[l][i][k][j] += temp2;
                }
            }
            //*proc.ostrm << endl;
            
        }
        ifile.close();
    }
    return;
}

#endif






