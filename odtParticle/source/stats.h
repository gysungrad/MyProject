/**
 * @file stats.h
 * Header file for classes stats 
 */

#ifndef STATS_H
#define STATS_H

#include "anyline.h"
#include "odtline.h"
#include "odtParam.h"
#include <string>

///////////////////////////////////////////////////////////////////////////////

/** Class implementing statistics gathering on the odt line such as means and variances.
 *  
 *  @author David O. Lignell
 */

class stats {

    ////////////////////// DATA MEMBERS /////////////////////

    public:

        int                 ngrd;          ///< Number of grid points. Constant grid assumed here
        int                 ngrdf;         ///< Number of grid faces. Constant grid assumed here
        double              Ldomain;       ///< stat domain size (same as odt)

        odtParam            *odtP;         ///< Pointer to #odtParam object of current realization

        std::vector<double> pos;           ///< vector of cell center positions
        std::vector<double> posf;          ///< vector of cell face  positions

        std::vector<double> uMean;         ///< mean u velocities
        std::vector<double> vMean;         ///< mean v velocities
        std::vector<double> wMean;         ///< mean w velocities

        std::vector<double> uMsqr;         ///< mean squares u
        std::vector<double> vMsqr;         ///< mean squares v
        std::vector<double> wMsqr;         ///< mean squares w

        std::vector<double> vTrans;        ///< a vector to tranfer odt grid to stat grid
    
#ifdef NEWSTATS
        /** The index of the uniform mesh korresponds to the locations 
         *  x(i) = dom * (i-1) / (N-1)
         * edstat(j,k,l,i)  ::  i = 0:N-1     ::  the index of the uniform mesh
         *                      j = 0:7       ::  the difference caused by eddy,
         *                                        by diffusion, by both, or by
         *                                        eddy adaption
         *                                        (del_eddy, del_eddy^2,
         *                                        del_fiff, del_diff^2,
         *                                        del_all, del_all^2,
         *                                        del_adpt, del_adpt^2)
         *                      k = 0:3       ::  the variable (u,v,w,T)
         *                      l = 0:Nstat-1 ::  the index of the data
         *                                        gathering periode
         * cstat(m,k,l,i)   ::  i = 0:N-1     ::  same as edstat
         *                      m = 0:2       ::  gathered value (u,u^2,du)
         *                      k = 0:4       ::  gathered value (u,v,w,T,eddy)
         *                      l = 0:Nstat-1 ::  same as edstat
         * ctime(l)         ::  l = 0:Nstat-1 ::  same as edstat
         * phstat(i,n,l)    ::  i = 0:N-1     ::  same as edstat
         *                      n = 1:Nphases ::  the index of the current phase
         *                      l = 0:Nstat-1 ::  same as edstat
         * phases(n)        ::  n = 1:Nphases ::  the index of the current phase
         * oldVars(i,k)     ::  i = 0:N-1     ::  same as edstat
         *                      k = 0:2       ::  same as edstat
         */
        
        std::vector< std::vector< std::vector< std::vector<double> > > > edstat;
        // data accumulation for calculation of budget terms     
        
        std::vector< std::vector< std::vector< std::vector<double> > > > cstat;
        // second variable for data accumulation                 
        
        std::vector<double>                                              ctime;
        // accumulation of time; used for calculating the average
        
        std::vector< std::vector<double> >                               oldVars;
        // variable for saving the old state to calculate the 
        // differenc caused by an eddy or a diffusion step
        
        std::vector< std::vector< std::vector<double> > >                phstat;
        // variable for saving changes of the phases in each stat grid
        // cell. Each phase gets its own index.
        // phstat = 1 : 100% of this phase
        // phstat = 0 :   0% of this phase
        
        std::vector<double>                                              phases;
        // variable to save the numbers of the phases used in odtline
        // temporary stats for statistics during diffusion
        
        /** 
         * cstats(m,k,ii)   :: ii = 1:ngrd    :: the index of the adaptive grid
         *                      m = 0:2       :: same as cstat
         *                      k = 0:3       :: variable (u, v, w, T)
         * phstats(n,ii)    :: ii = 1:ngrd    :: same as cstats
         *                      n = 1:Nphases :: same as phstat
         */
        std::vector< std::vector< std::vector<double> > >                cstats;
        // temporary variable for statisic written to cstat
        std::vector< std::vector<double> >                               phstats;
        // temporary variable for statisic written to phstat
        std::vector< std::vector<double> >                               tau_wstat;
        double  meanTauW, mTauWlo, mTauWlo2, mTauWup, mTauWup2;
        // statistic of wall stress. collects the statistic of tau_w.
        
        // variables used for change of domain size during the simulation
        std::vector<int>     ngrd_av;
        // saves ngrd for each periode
        int                  max_ngrd;
        // saves the maximum of ngrd_av
        std::vector<double>  Ldomain_av;
        // saves the domain sizes for each period
        double               max_Ldomain;
        // saves the maximum of the Domain sizes
#endif
        
    ////////////////////// MEMBER FUNCTIONS  /////////////////////
        
        void updateMeans(odtline &odtl, const double &t, const double &dt);
        void outputProperties(std::string fname, const int iStat);
        void initStats(odtline *odtl);
        double getUBulk(int i, int istat);
        
#ifdef NEWSTATS
        void initBStats(odtline &odtl);
        void BStats(odtline &odtl, const double tStep, const int iStat);
        void BStatsTime(odtline &odtl, const double tStep, const int iStat);
        void BStatsSpace(odtline &odtl, const double tStep, const int iStat);
        void BStats2statGrd(odtline &odtl, const int iStat);
        void BSetOld(odtline &odtl);
        void BChange(odtline &odtl, const int &jj, const int &iStat);
        void BSnap(const int myid);
        void BLoadStats(const int myid);
        
        void ChangeStatGrid(const int iStat);
#endif
        void odtGrd2statGrd(std::vector<double> odtposf,
                      std::vector<double> odtVec, double pJump=0.0);
    private:
        void wrapAround(std::vector<double> &odtposf,                             
                      std::vector<double> &odtvec, double &pJump);
        
        void odtGrd2statGrd_c(std::vector<double> odtposf, std::vector<double> odtpos,
                      std::vector<double> odtVec, double pJump=0.0,
                      int BCl=1, double al=0.0, double bl=0.0, double cl=0.0,
                      int BCu=1, double au=0.0, double bu=0.0, double cu=0.0);
        
        void addToVec(std::vector<double> &vecBase,
                      std::vector<double> &vecToAdd,
                      const double &t, const double &dt);
        void addToVec2(std::vector<double> &vecBase,
                      std::vector<double> &vecToAdd,
                      const double &t, const double &dt);
        
        
    ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:

        stats(odtParam &odtpp, double Ld=1.0, int npts=100);               // constructor
#ifdef NEWSTATS
        stats(const stats *oldStats, const int iStat, double newLdomain);  // 2nd constructor for change file usage
        void operator=(const stats &newStats);  // =-Operator used for 2nd constructor
#endif
        ~stats(){ odtP=0; }



};

#endif
