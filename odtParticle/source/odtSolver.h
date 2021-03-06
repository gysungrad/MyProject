/**
 * @file odtSolver.h
 * Header file for class odtSolver
 */

#ifndef ODTSOLVER_H
#define ODTSOLVER_H

#include <iostream>
#include <vector>
#ifdef BOOST
#include <boost/ptr_container/ptr_vector.hpp>
#endif
#include <fstream>
#include <string>
#include "anyline.h"
#include "odtParam.h"
#include "eddy.h"
#include "diffuser.h"
#include "stats.h"
#include "randomGenerator.h"
#include "streams.h"
#include "particles.h"
#include "ETA.h"
#include "ETA_Aq.h"
#include "MOM.h"

#include "point.h"
#include "subdomain.h"
#include <algorithm>
#include "tripletmap.h"
#include "eddy_functions.h"
#include "stopevent.h"
#include "sdgroup.h"
#include <queue>
#include "walltime.h"


#ifdef DOCANTERA
#include "IdealGasMix.h"
#include "transport.h"
#endif

#include "cantera_shell_functions.h"

using Cantera::Transport;
#ifdef CANTERA18
using Cantera_CXX::IdealGasMix;
#else
using Cantera::IdealGasMix;
#endif

using namespace std;

///////////////////////////////////////////////////////////////////////////////

/** A class implementing the ODT code solution.  This object is created in main
 *  and solution is performed with calculateSoln.  This class is the main driver
 *  program, containing data members of most of the other classes.  The main function
 *  sets up and calls this object.
 * 
 *  For spatial cases, a "time" increment is really a "space" increment.
 *  
 *  @author David O. Lignell
 */

class odtSolver {

  ////////////////////// DATA MEMBERS /////////////////////

    public:

        /// Don't change the order of these variables here (for the constructor init list)!

        odtParam             odtP;            ///< parameters object
        odtline              *mainLine;       ///< points to odtline for temporal, odtline for spatial
#ifdef COMPSCI
        odtline              *eddyLine;
#endif
        randomGenerator      rndmgn;          ///< random generator object
        diffuser             *lineDiffuser;   ///< diffuser/reactor object
        eddy                 ed;              ///< eddy object
#ifdef COMPSCI
        eddy                 ed3;             ///< testThirds test eddy
#endif
        stats                odtStats;        ///< statistics object 
        particles            *part;           ///< particles object
        table                *luTable;        ///< pointer to lookup table
        ETA*                 etaTools;        ///< Pointer to eta object
        MOM*                 momTools;        ///< Pointer to mom object
        
        double               dtSmean;         ///< mean eddy sample time
        double               PaSum;           ///< sum of Pa of eddies
        int                  nPaSum;          ///< number going into PaSum
        int                  neddies;         ///< number of eddies accepted
        
        double               PaSumC;          ///< sum of Pa of eddies
        int                  nPaSumC;         ///< number going into PaSum
        
        vector<double>       lastDA;          ///< constant (unif) mesh to list time of last adapt
        
        bool                 Ltune;           ///< tune mingradfrac for mesh adptn for varing turb. intensity
        int                  fewestCells;     ///< for tuning 
        int                  iLES;            ///< number of Eddy suppressions 
        bool                 breakReali;      ///< flag for breaking the current reali and start with a new one
        vector<vector<double> > breakUpStats; ///< 
#ifdef NEWSTATSF
        vector<double>       eddyRegion;
#endif
	
	bool		    LhasVel;  
//variables for subdomain decomposition
        double              maxsimplemapsize;
        double              limit;
        double              solregsize, guardsize, subdomainsize;
        int                 Nsubdomains;
        int                 stopnumber;
        double              eddycount;
        double              timeinsd, timeincontrol;
        double              t1,t2,clockzero,stopeventtime;
        double              lambda;
        int                 distributed_maps;
        //        
//         odtParam            middlesdparam, rboundarysdparam, lboundarysdparam;
         std::vector<std::priority_queue<tripletmap,std::vector<tripletmap>, tripletmap::greater > > tripletmaps;                                              //(for every subdomain, all of the tripletmaps which will be applied to the domain)
         std::vector<std::priority_queue<point,std::vector<point>, point::greater > > points;          //(for every subdomain, two heaps (left,right) in
        //which the points of already inserted triplet maps are stored)
         std::vector<stopevent>       stopevents;
#ifdef BOOST        
	 boost::ptr_vector<subdomain>       subdomains_comb;
#endif	 
 	std::list<double>  acquired_stopnumbers; 

	walltime               wt;
	
        ////////////////////// MEMBER FUNCTIONS  /////////////////////
        
        void   calculateSoln();                
	
	//subdomain decomposition
	
         void apply_event(int stopnumber);
        
    private:
        
        bool   sampleEddy(int &iStart, int &iEnd, double t0, double time, const int &iStat);
        double computeDtSmean();                                                           
        double sampleDt();                                                                 
        void   raiseDtSmean(int &iEtrials);                                                
        void   lowerDtSmean(double time);                                                             
        double computeDtCUmax();                                                           

        void   updateDA(const double &time, double &tLastDA,                               
                        int &cLastDA, int iStart, int iEnd);   
        bool   adaptAfterSufficientDiffTime (const double &time,                           
                                             double       &tLastDA,
                                             int          &cLastDA,                
                                             double       &dtCUmax ); 
        void   adaptEddyRegionOfMesh(int &iStart, int &iEnd, double &dtCUmax,              
                                     const double &time, double &tLastDA, int &cLastDA);

        void   outputHeader();                                                             
        void   outputProgress(double &time, double &t0, int &neddies, int &iEtrials);      
        void   dumpLine(string fname,
                        ofstream &gnufile,
                        ofstream &gnufilePart);
        void   dumpStats(string fname, ofstream &gnufile, const int iStat);
        void   dumpStatsEndIstat(const int &iStat);                                        

        // void   setRestartState(string odtlRstFile, string partRstFile);
        void   setRestartState(string odtlRstFile);

        void   tuneMinGradFrac(int iStart, int iEnd);

        bool   testThirds();          //large eddy suppression
        bool   testLES_elapsedTime(double time, double tauEddy); // large eddy suppression
        bool   testLES_fracDomain(double eSize);
        bool   testLES_integralLength(double time, double eSize);

        // needed for periRestart
        void   dumpStatsToRestart(void);
        void   setPeriRestartState(std::string odtlRstFile, std::string partRstFile);
        
        // needed for reading change file
        void   ReadChangeFile(const int iStat);
        void   IncreaseDomain(ifstream *ifile, char LUC);
        void   ChangeParameter(ifstream *ifile);
        void   SetUMean(ifstream *ifile);

        // initializing odtSolver in the constructor and re-set for new realizations
        void init();
        // some private pointers needed for init()
        streams              *combStrms;
        ETA_Aq               *ETA_Aq_obj;
        IdealGasMix          *gas;
        Transport            *tran;        

	
	//subdomain decomposition
	void    create_subdomains_from_mainline(); //subdomain decomposition

        bool insert_map(double leftEdge, double rightEdge, double time);
        bool check_guard_region(std::priority_queue<point,std::vector<point>, point::greater > &points,double start, double end);
        bool is_dependent_map(double lStart, double lEnd, vector<int> &domains_affected);
        void return_subdomains(double lStart, double lEnd, std::vector<int> &domains_affected, bool neighbors=true);
        void insert_stop_event(tripletmap &map, vector<int> &sdindex, bool timeevent=false,int type=0);                                                  //is necessary
        void gather_sol(odtline &cmbl, double lstart,double lend, vector<int> sdindex);
        void gather_sol(odtline &cmbl,vector<int> sdindex);

        void distribute_sol(sdgroup &group, bool wrap=false);
	
        void distribute_maps(subdomain &sd);
        void update_mainLine(int sd);
        double modL(const double &inL);
        bool is_in_guardregion(const double &inp, const subdomain *sd, bool left);
        void vector_out(std::vector<double> &v);
        bool check_time_contraint(int sdindex, double time, bool &leftok, bool &rightok);
        void updt_guardregion(int sdindex, double time,bool guardleftok, bool guardrightok);
	
	void shiftFlame(/*double x_0, double t_0, double t, std::ofstream &sL_file, int iEtrials*/);
	
	
        ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:

        odtSolver(string odtParamFile,
                  string bcFile,
                  string odtlRstFile,
//                  string partRstFile,
                  string gasSetupFile,
                  string partSetupFile,
                  std::string sootTableFile,
                  streams     &combStrms,
                  IdealGasMix &gas,
                  Transport *tran);                                                        

        ~odtSolver(){ 
            if(mainLine) delete mainLine; 
#ifdef COMPSCI
            if(eddyLine) delete eddyLine;
#endif
            if(lineDiffuser) delete lineDiffuser;
            if(part) delete part;
            if(luTable) delete luTable;
            if(etaTools) delete etaTools;         
        };                       // destructor                                 


};

#endif
