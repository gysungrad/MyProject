/**
 * @file odtSolver.cc
 * Header file for class odtSolver
 */

#include "odtSolver.h"
#include "inputFile.h"
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <string>
#include <iomanip>
#include "processor.h"
#include "ETA_tableLookup.h"
#include "ETA_default.h"
#include "MOM_soot.h"             
#include "MOM_default.h"
#include "ETA_Aq.h"
#include "mom_Aq.h"

using namespace std;
#ifdef BOOST
using namespace boost;
#endif

using Cantera::Transport;
#ifdef CANTERA18
using Cantera_CXX::IdealGasMix;
#else
using Cantera::IdealGasMix;
#endif

extern processor proc;

///////////////////////////////////////////////////////////////////////////////
/** Constructor. 
 *
 *  @param odtParamFile  \input main parameter input file name.
 *  @param bcFile        \input boundary condition input file name.
 *  @param odtlRstFile   \input odtline restart file name.
 *  @param partRstFile   \input particle restart file name.
 *  @param caseSetupFile \input profile file (mixture fraction field for init).
 *  @param partSetupFile \input particle setup file
 *  @param combStrms     \input streams object for pointer for odtline.
 *  @param gas           \input Cantera gas object for odtline.
 *  @param tran          \input Cantera transport object for odtline object.
 */
 
odtSolver::odtSolver(string        odtParamFile,
                     string        bcFile,
                     string        odtlRstFile,
//                      string        partRstFile,
                     string        caseSetupFile,
                     string        partSetupFile,
                     string        sootTableFile,
                     streams       &_combStrms,
                     IdealGasMix   &_gas,
                     Transport     *_tran) :
                     
    //---------- constructor initialization list

    odtP         (odtParamFile, bcFile),
#ifdef COMPSCI
	eddyLine     (0),
#endif
    rndmgn       (odtP.seed +proc.myid*1000), // needed to set a differend seed for each CPU
    odtStats     (odtP, odtP.domainLength, odtP.nsgrd),
    combStrms	 (&_combStrms),
    gas          (&_gas),
    tran	 (_tran),
    tripletmaps  (odtP.numSubdomains),//subdomain decomposition
    points       (2*odtP.numSubdomains)//subdomain decomposition
{
    //-------------------------------------------Initialization
    mainLine     = 0;
    part         = 0;
    etaTools     = 0;
    momTools     = 0;
    luTable      = 0;
    ETA_Aq_obj   = 0;
    momTools     = 0;
    lineDiffuser = 0;
    
    init();

#ifdef NEWSTATSF
    eddyRegion = vector<double>(mainLine->ngrd, 0.0);
#endif
    
    iLES       = 0;

    if(odtP.LperiRestart){
        dumpStatsToRestart(); // writing initial periRestart file
    }

    breakReali = false;
    // breakUpStats(i,j):   i = number of realisation;
    //                      j = (time of break up, size of eddy, size of droplet)
    breakUpStats = vector<vector<double> >
        (odtP.nOdtReals, vector<double>(3, 0.0));
}

///////////////////////////////////////////////////////////////////////////////
/** init() initializes the odtSolver. It is basically a copy of the old
 * constructor for odtSolver and is called from the constructor. In addition,
 * init() serves as a re-init of the odtSolver for a new realization. 
 * It does that by completely (hopefully) deleting the old odtSolver 
 * structures (e.g. mainline, part, etc.) and re-allocating them
 */
void odtSolver::init()
{ 
    string odtlRstFile   = "../input/"+ inputFileDir + odtP.restartOdtL;
    string partRstFile   = "../input/"+ inputFileDir + odtP.restartPart;
    string caseSetupFile = "../input/"+ inputFileDir + odtP.caseInp;
    string partSetupFile = "../input/"+ inputFileDir + odtP.partInp;
    string sootTableFile = "../input/"+ inputFileDir + odtP.tableInp;

    //-----------------------------------------------------------------------------
    //-------- Specify line definitions
    //-------- Spatial is just one line with velocity and with or without reaction
    //-------- Temporal is one line or two, with reaction if two and with or without rxn if one
    ed.rnd = &rndmgn; //doldb

    LhasVel=!odtP.Llem || odtP.Iparticles;
    if (mainLine) delete mainLine;
    mainLine = new odtline(odtP.ngrd_0, odtP.domainLength, gas, tran, combStrms, 
                           &odtP, caseSetupFile, LhasVel, odtP.Lrxn);

    //-------------------------------------------Initialization

    part = 0;
    etaTools = 0;
    momTools=0 ;
    luTable = 0;
    
    if(odtP.Lrestart) //GSDB
       setRestartState(odtlRstFile);
//        setRestartState(odtlRstFile, partRstFile);
    
    if(odtP.Iparticles)
        part = new particles(odtP.Iparticles, 
                &odtP, mainLine, partSetupFile);
    
    // ---------------------------Create ETA object Set the required pointers to line, diffuser and eta   

    if (odtP.Ieta && odtP.IetaType == odtP.IETATYPE_DEFAULT) {
        if (etaTools) delete etaTools;
        etaTools = new ETA_default(mainLine);
    } else if (odtP.Ieta && odtP.IetaType == odtP.IETATYPE_TABLELOOKUP) {
        if (luTable) delete luTable;
        luTable = new table(mainLine);
        if (etaTools) delete etaTools;
        etaTools = new ETA_tableLookup(mainLine, luTable);
        etaTools->updateOdtLineVecs(true);
    }
    else if (odtP.Ieta && odtP.IetaType == odtP.IETATYPE_AQUEOUS) {
        if (etaTools) delete etaTools;
        etaTools = new ETA_Aq(mainLine);
    }
    
    // ---------------------------Set the lines etaTools to be this etaTools
    mainLine->etaTools = etaTools;
    
    if(odtP.Imom && odtP.ImomType==odtP.IMOMTYPE_DEFAULT) {
        if (momTools) delete momTools;
        momTools = new MOM_default(mainLine);
    }
    else if(odtP.Imom && odtP.ImomType==odtP.IMOMTYPE_AQUEOUS) {
        if (ETA_Aq_obj) delete ETA_Aq_obj;
        ETA_Aq_obj =   new  ETA_Aq(mainLine)   ;
        if (momTools) delete momTools;
        momTools = new   mom_Aq(ETA_Aq_obj, mainLine);
    }
    else if(odtP.Imom && odtP.ImomType==odtP.IMOMTYPE_SOOT) {
        if (momTools) delete momTools;
        momTools = new MOM_soot(mainLine, sootTableFile);
    }
    
    // ---------------------------Set the line's momTools to be this momTools
    mainLine->momTools = momTools;
    
    // ---------------------------Construct particles, particles may need diffuser for construction, 
    //                            so pass it as a parameter
    
    if (odtP.Iparticles)
    {
        if (part) delete part;
        part = new particles(odtP.Iparticles,
            &odtP, mainLine, partSetupFile);
    }
    // ---------------------------Construct diffuser
    if (lineDiffuser) delete lineDiffuser;
	  lineDiffuser = new diffuser(mainLine, odtP,part);

    // ---------------------------Set the particles's diffuser
    if(part) {
        part->diff = lineDiffuser;
    }
    
    // ---------------------------Read in initial particle profile if flagged to.
    
    if (part && part->LInitProfile) {
        part->readInitialProfile();
    }
    
    // ---------------------------Set the diffuser to etaTools
    if(odtP.Ieta)
        etaTools->diff = lineDiffuser;
    if(odtP.Imom) 
        momTools->diff = lineDiffuser;
        
    // --------------------------initialize subdomains
    
#ifdef BOOST
    if(odtP.Lsubdomain) {
      Nsubdomains=odtP.numSubdomains;
        for (int i=0; i<Nsubdomains; i++)
        {
            subdomains_comb.push_back(new subdomain(0,1,0,1,0.0,i, odtP, this, gas, tran, combStrms));
            //cout << "odtS C           : &sd("<<i<<")   :  "<< &(subdomains_comb.at(i)) << endl; cout.flush();
	    
        }
        stopnumber = 0;
    }
#endif    
    
    //-----------------------------------------------------------------------------

    PaSum      = 0.0;
    nPaSum     = 0;
    dtSmean    = computeDtSmean();
    neddies    = 0;
    PaSumC     = 0.0;
    nPaSumC    = odtP.Llem ? 1 : 0; 
    Ltune      = false;

    // Reading restart file, setting this state and 
    // writing initial periRestart file, if needed.
//     if(odtP.Lrestart){
//         setRestartState(odtlRstFile, partRstFile);
//     }
    lastDA     = vector<double>(odtP.sLastDA,  0.0);

    if(odtP.Lperiodic && proc.myid==0) 
        if(odtP.Lrxn && proc.myid==0) {
            cout << endl << "***** ERROR, periodic domains won't work with combustion yet" << endl;
            exit(0);
        }
}


///////////////////////////////////////////////////////////////////////////////
/** Solve the odt line.  This is the main driver program.
 * 4 nested loops: Loop over realizations (1), over stat gathering periods (2),
 * over sub intervals (3), and eddy advancement within a sub interval (4). 
 * tMark goes with (3).  dtSample goes with (4).
 *\vc{
 *      2              2              2              2              2
 *      |    3    3    |    3    3    |    3    3    |    3    3    |
 *  1   |::::|::::|::::|::::|::::|::::|::::|::::|::::|::::|::::|::::|
 *      |4444 4444 4444|4444 4444 4444|4444 4444 4444|4444 4444 4444|
 *
 *      |              |              |              |              |
 *  1   |::::|::::|::::|::::|::::|::::|::::|::::|::::|::::|::::|::::|
 *      |              |              |              |              |
 * }
 * Advance in time only sampling eddies till get one (EE).  Then diffuse the system
 * until you catch up to the eddy time:
 * <code><pre>
 *     last EE    this EE       diffuse to catch up
 * .....|..........|           ---------------------->  ................|
 *      t0         time                                                t0,time
 *                          (last EE)   this EE         next EE
 * Then advance some more ......;.........|................|         etc.
 *                                        t0              time
 * </pre></code>
 * Eddy timesteps are smaller than required diffusive steps, so if we were to
 * lock-step the two processes we would do too much diffusive work (i.e., take
 * more diffusive steps than needed)
 */
 
void odtSolver::calculateSoln() {

    if(!odtP.Llem) odtStats.initStats(mainLine);

    string fileName;
    fileName = proc.dataDir+"plotmeans.gnu";    ofstream gnufileMeans(fileName.c_str());
    fileName = proc.dataDir+"plot_odt.gnu";     ofstream gnufile(fileName.c_str());
    fileName = proc.dataDir+"plot_part.gnu";    ofstream gnufilePart;
    if(odtP.Iparticles)
        gnufilePart.open(fileName.c_str());
    if(odtP.Lrxn) 
        gnufile << "plot '" << "odt_init.dat' us 1:8; pause -1;" << endl;
    else
        gnufile << "plot '" << "odt_init.dat' us 1:6; pause -1;" << endl;

    //---------- set tuning of smallGradFrac for mesh adaption //tune

    if(odtP.smallGradFrac < 0.)  { //negative flags tuning of smallGradFrac
        Ltune = true;
        odtP.smallGradFrac *= -1.;
    }

    double posLower, posUpper;

    
    //------------------------------ loop over realizations

    for(int iReals=1; iReals<=odtP.nOdtReals; iReals++) { 
        
        breakReali      = false;
        if (iReals > 1)
        {
          cout << endl
            << "#-----------------------------------------------------------------"  
            << endl
            << "# resetting odtSolver for the realization no.: " << iReals << endl
            << "#-----------------------------------------------------------------" 
            << endl;

          init();
        }

        dtSmean         = computeDtSmean();   // initial mean eddy sample time
        double dtMark   = odtP.tEnd/(odtP.nTseg*odtP.nStat);
        double tMark    = 0.0;                // statistics: sub interval time
        int    iTime    = 0;                  // statistics: sub interval counter
        double time     = odtP.Llaminar? dtMark : sampleDt();         // current eddy sample time
        double tLastDA  = 0.0;                // time of last diffusive mesh adaptation
        double t0       = 0.0;                // current diffused time
        double dtDiffuse= 0.0;                // diffusion time (catch up to eddy)
        bool   LEA      = 0;                  // flag eddy acceptance
        int    iEtrials = 0;                  // number of eddy trials
        double dtCUmax  = 0.0;                // max time before catch up diff/eddy
        int    iStart   = 0;                  // eddy start index
        int    iEnd     = 0;                  // eddy end index

        neddies         = 0;
        PaSum           = 0.0;                // sum of Pa of eddies
        nPaSum          = 0;                  // number going into PaSum
        dtCUmax         = computeDtCUmax();

        stringstream ss1;  string s1;         // eddy number (for filename, # -> string)
        stringstream ss2;  string s2;         // realization number (for filename, # -> string)
        stringstream ss3;  string s3; // temp

        double timeOffset;                    // time - timeOffset = t since iStat interval


	// Warning: mainLine should be reinitialized for each realization here!
	
        //---------- variables for adaptation after sufficient time

        int cLastDA     = 0;
        lastDA          = vector<double>(odtP.sLastDA,  0.0);

        fewestCells = mainLine->ngrd; //tune: initialize fewestCells for this iteration


        //---------- do an initial mesh adaption
        
        *proc.ostrm << endl << "#*********** Adapting the initial condition before beginning" << endl;
        if(mainLine->LhasRxn) mainLine->setTempVec(); 
        
        else if(odtP.ItableLookup) mainLine->etaTools->updateOdtLineVecs();

        mainLine->meshAdapter.adaptGrid(0,mainLine->ngrd-1);
        dumpLine("init_adpt.dat", gnufile, gnufilePart);

        if(odtP.ItableLookup){
            mainLine->etaTools->updateOdtLineVecs();
        }
        if (odtP.Iparticles) {
            mainLine->setVoidFrac(this->part);
        }
        //--------------------------------------------------
        // reading periRestart file
        cout << endl << "Check restart: " << (odtP.LperiRestart && !(odtP.Llem)) << " " << odtP.Lrestart;
        if (odtP.LperiRestart && !(odtP.Llem)){
            //ss1.clear(); ss1 << iReals; ss1 >> s1;
            //mainLine->outputProperties("01_Before_SetPeriRestarti_" + s1 + ".dat");
            ss1.clear(); ss1 << proc.myid; ss1 >> s1;
            setPeriRestartState("../input/restart_odtl_" + s1 + ".dat",
                                "../input/restart_particles.dat");
            //ss1.clear(); ss1 << iReals; ss1 >> s1;
            //mainLine->outputProperties("02_After_SetPeriRestart_" + s1 + ".dat");
            if (odtP.domainLength != mainLine->Ldomain)
                odtP.domainLength = mainLine->Ldomain;
        }
        
#ifdef NEWSTATSF
        // creating file for eddy sequence
        ss1.clear(); ss1 << iReals; ss1 >> s1;
        fileName = proc.dataDir+"EddySequence_"+s1+".dat";
        ofstream eddyfile(fileName.c_str());
        eddyfile << "#               N             time"
                 << "         leftEdge        rightEdge";
#endif
        

//---------- subdomain decomposition
#ifdef BOOST	
	if(odtP.Lsubdomain) create_subdomains_from_mainline(); 
	
	//---------- variables for subdomain decomposition
	bool privdone;
	int countapply=0;
	double simplemapsaverage=0;

	time_t calendar_time;
	struct tm todays_date;

	ptr_vector<subdomain>::iterator it;
	bool alldone=false;
	bool overtheedge=false;
	int outcount =0;
	vector<int> domains_comb(Nsubdomains);
	//ptr_vector<subdomain> domainspointer_comb(Nsubdomains); // -> use domains_comb instead
	int i;
	//ptr_vector<subdomain> stopdomains_comb(Nsubdomains); // -> use stopdomains_index instead
	vector<int> stopdomains_index(Nsubdomains); 
	tripletmap newmap;
	for(i=0;i<Nsubdomains;i++){
	    domains_comb.at(i)=i;
	}
#endif     
        
        //------------------------------ loop over statistics periods

        for(int iStat=1; iStat<=odtP.nStat; iStat++) {
            
            *proc.ostrm << endl << "#***************** BEGINNING iStat region " << iStat << endl;
            iStart = 0; iEnd = 0;
            timeOffset = odtP.tEnd/odtP.nStat*(iStat-1);  // for stats

#ifdef NEWSTATS
            if(odtStats.ngrd_av.size() != 1 && !(odtP.Llem)){
                *proc.ostrm << endl << "#**************** changing stats";
                odtStats.ChangeStatGrid(iStat);
            }
            // reading change file if existing
            if(!(odtP.Llem)){
                ReadChangeFile(iStat);
            }
            //ss1.clear(); ss1 << iReals; ss1 >> s1;
            //ss3.clear(); ss3 << iStat; ss3 >> s3;
            //mainLine->outputProperties("03_After_ReadChangeFile_" + s1 +"_"+ s3 + ".dat");
#endif
            //------------------------------ loop over statistics sub-intervals

            for(int iTseg=1; iTseg<=odtP.nTseg; iTseg++) {

              *proc.ostrm << endl << "#----------------- Stat interval " << iTseg << endl;

                tMark += dtMark;
                if(odtP.Llaminar) time = tMark; // to bypass while-loop
                iTime++;

                //------------------------------ proceed to tMark

                *proc.ostrm << endl << "# dtSmean = " << dtSmean << endl;

                while(time <= tMark) {

		  if(odtP.Lsubdomain) {
#ifdef BOOST		   
		      double dtime= time;
		      int mapstotal =0;
		      distributed_maps=0;
		      double dtStep= odtP.dtGatherSubdomains;
		      int dtstart = t0 /dtStep;
		      for(int dt = dtstart; dt< tMark/dtStep; dt++){
			  for(i=0;i<Nsubdomains;i++){
			      subdomains_comb.at(i).curtime=t0;
			      subdomains_comb.at(i).t0left=t0;
			      subdomains_comb.at(i).t0right=t0;

			  }
			  //subdomains_comb.at(5).printcounter=printcounter_old;


			  dtMark=(dt+1) * dtStep;
			  t1=wt.walltime_return(&clockzero);

			  while(dtime < dtMark && !overtheedge) {
			      double dtSample = sampleDt();
			      ed.sampleEddySize(odtP,rndmgn);
			      cout <<  endl <<"eddy size " << ed.eddySize << endl; cout.flush();
			      ed.sampleEddyPosition(odtP, *mainLine, rndmgn);   


			      overtheedge=insert_map(ed.leftEdge,ed.rightEdge,dtime);
			      mapstotal++;


			      if(dtime>=dtMark&&!overtheedge){ // this is never true ! // D.Nolte
				  newmap=tripletmap(0,0,dtime);
				  newmap.applymap=false;              
				  for(i=0;i<Nsubdomains;i++){
				      //                        stopdomains_comb.at(i)=(subdomains_comb.at(i));
				      stopdomains_index.at(i) = i; // = subdomains_comb.at(i).sdnumber                        
				  }               
				  insert_stop_event(newmap,stopdomains_index,true,4);
			      }    

			      dtime+=dtSample;
			  }
			  overtheedge=false;


			  //list<tripletmap>::iterator lit;
			  // int j=0;
			  for (i=0; i<subdomains_comb.size(); i++) {
				
			      distribute_maps(subdomains_comb.at(i));
			      
			      
			  }
			
			  //distributed_maps+=it->maps.size();
			  //                                 if(it->sdnumber== int(Nsubdomains/2)) {for(lit = it->maps.begin(); lit!=it->maps.end(); lit++){
			  //                                     cout<<endl << "sub " <<it->sdnumber <<" map "<< j << " mapstart;mapend "<<lit->lStart<< " " << lit->lEnd <<" time " <<lit->time;
			  //                                     j++;
			  //                                 }}

#ifdef NEWSTATSF
                        eddyfile << scientific << setprecision(8) << endl
                            << setw(17) << neddies
                            << setw(17) << time
                            << setw(17) << ed.leftEdge
                            << setw(17) << ed.rightEdge;
#endif
			  *proc.ostrm << endl << "maps distributed";
			  alldone=false;
			  while(!alldone){
			    
			    
			      
			      alldone=true;
			      //apply_maps returns true when the subdomain maplist is empty
			      //thus, when all subdmains are out of maps notalldone will be
			      //true at the end of the loop.
			      int i;
			      timeincontrol+=wt.walltime_return(&t1);
			      cout << endl << "alldone " << alldone; cout.flush();
	      #pragma omp parallel for private (i,privdone) //schedule(dynamic,1)
			      
			      
			      for(i=0;i<Nsubdomains;i++){
				
				  //cout << "-----------------------------------------------------" << endl;
				  *proc.ostrm << endl << "diffuse and apply maps; time; mapstotal; distributed_maps: " << dtime << " " << mapstotal<< " " <<distributed_maps; 

				  *proc.ostrm << endl << "subdomain " << i<< " Stopnumber "<< subdomains_comb.at(i).stopnumber<< endl;(*proc.ostrm).flush();
				  if( subdomains_comb.at(i).stopnumber<0){
					
			      
				      privdone=subdomains_comb.at(i).apply_maps();
				      
				      				    
				      cout << endl << "privdone " << privdone; cout.flush();
				      cout << endl << "alldone " << alldone; cout.flush();
				      simplemapsaverage+=subdomains_comb.at(i).nmapsindependent;
				      cout << endl<< "maps applied"<<endl;
				      
// 				     

				  }
				  else {

				      privdone=false;
				  }
				  //    #pragma omp critical (crit2)
				  {
				      alldone*=privdone;
				  }



				  countapply++;
				  *proc.ostrm << endl << "alldone "<< alldone<< " countapply: " << countapply << endl; (*proc.ostrm).flush();
			      } // end for i<Nsubdomains
			      

			  } // end while(!alldone)

			  cout <<endl<< "catch up diffusion ";
			  for(i=0;i<Nsubdomains;i++){
			      subdomains_comb.at(i).diffuse(dtime);
			      cout <<endl << "sub, dtime, curtime " << i << " " <<  dtime << " " << subdomains_comb.at(i).curtime;
			      subdomains_comb.at(i).curtime=dtime;
			  		      
			  }


			  countapply=0;
			  //                        
			  //            domainspointer_comb = get_locations(domains_comb); // deprecated !! 
			  //            at this point, domains_comb == 0...Nsubdomains-1   (cf. l.258)
			  gather_sol(*mainLine,mainLine->posf[0],mainLine->posf[mainLine->ngrd],domains_comb);
			  //DEBUG
			
			  
			  cout <<  "gather sol" << endl; cout.flush(); 
			  mainLine->setTempVec();
			  cout <<  "set TempVec" << endl; cout.flush();
			  cout<<" mainLine->ngrd "<<mainLine->ngrd<<endl; cout.flush(); 
			  ss1.clear();
			  ss1 << dtime; ss1 >> s1; 


			  
			  if(odtP.shiftMethod>0) {
			      *proc.ostrm << endl << "shift flame";(*proc.ostrm).flush();
			      shiftFlame(/*x0,t0,dtime, sLfile, iEtrials*/);
// 			      if(odtP.binary && iEtrials%odtP.modDisp==0) {		//NetCDF output; TO DO
// 				  mainLine->writeNcTimeData_flameSpeed(nc_T_flameSpeed, dtime);
// 				  ncflameSpeed.sync();
// 			      }
			  }


			  mainLine->meshAdapter.adaptGrid(0, mainLine->ngrd-1);

			  cout<<" mainLine->ngrd "<<mainLine->ngrd<<endl; cout.flush(); 
			  if(iEtrials%odtP.modDump==0) dumpLine(s1+"_adapt.dat", gnufile, gnufilePart);
			  //subdomains_comb.at(5).print_subdomain();
			  //int printcounter_old=subdomains_comb.at(5).printcounter;

			  t0 = dtime;
			  
			  //---------------NetCDF output; TO DO
// 			  if(odtP.binary && iEtrials%odtP.modDisp==0) {//&& dmpr.n_time > 0 && dmpr.checkDumpTime(dtime)) {
// 			      *proc.ostrm << endl << "write data; time: "<< dtime;
// 			      mainLine->writeNcTimeData_props(nc_XT_propsComb, dtime,odtP.dxmin/odtP.domainLength, &odtP);
// 			      ncpropsComb.sync();
// 
// 			      //                             subdomains_comb.at(1).odtl.writeNcTimeData_props(nc_XT_propsComb_sub, dtime,odtP.dxmin/odtP.domainLength/Nsubdomains);
// 			      //                             ncpropsComb_sub.sync();
// 			      mainLine->writeNcTimeData_spec(nc_XT_specComb, dtime,odtP.dxmin/odtP.domainLength);
// 			      ncspecComb.sync();
// 
// 
// 			  }
			  
			   
			  if(odtP.shiftMethod>0) create_subdomains_from_mainline();

			  iEtrials++; 
			  } // end for (dt = dtstart; dt< tMark/dtStep; dt++)


// 			  calendar_time = time(NULL);
// 			  todays_date = *localtime(&calendar_time);
// 			  cout << " WRAPPED " <<  dtime << " at time "<< todays_date.tm_hour << ":" << todays_date.tm_min << ":" << todays_date.tm_sec << endl;


			  time=dtime;
#endif
		      }  // end if (odtP.Lsubdomain)
		  else{
			if(time-t0 >= dtCUmax) {           // diffusion catchup block;
							  // do if no eddies in a long time

    #ifdef NEWSTATS
			    if(!(odtP.Llem)) odtStats.BSetOld(*mainLine);
    #endif

			    dtDiffuse = time-t0;
			    *proc.ostrm << endl << "# Catching up diffusion to eddies: time = " << time; 

			    posLower = mainLine->posf[iStart];   // to update iStart, iEnd below
			    posUpper = mainLine->posf[iEnd+1];
			    lineDiffuser->diffuseOdtLine(dtDiffuse, t0-timeOffset, odtStats, iStat);
			    iStart = mainLine->linePositionToIndex(posLower, true); 
			    iEnd   = mainLine->linePositionToIndex(posUpper, false);
			    
 #ifdef NEWSTATSF
			    if(!(odtP.Llem)){
				odtStats.BChange(*mainLine, 2, iStat);
				odtStats.BChange(*mainLine, 4, iStat);
				odtStats.BSetOld(*mainLine);
			    }
#endif
			    
			    
			    adaptAfterSufficientDiffTime(time, tLastDA, cLastDA, dtCUmax);
			    
			    if(odtP.ItableLookup){
				mainLine->etaTools->updateOdtLineVecs();
			    }

			    t0 = time;
    
#ifdef NEWSTATSF
                        if(!(odtP.Llem)){
                            odtStats.BChange(*mainLine, 6, iStat);
                            //odtStats.BChange(*mainLine, 4, iStat);
                        }
#else
  #ifdef NEWSTATS
                        if(!(odtP.Llem)){
                            odtStats.BChange(*mainLine, 2, iStat);
                            odtStats.BChange(*mainLine, 4, iStat);
                        }
  #endif
#endif
                    }

                    //---------- sample next eddy

                    // Since eddies are implemented at t0, we will pass t0 to testLES_elapsedTime()
                    // This also helps with the particle history timing in tripletMapParticles()
                    // JCH changed on 7/30/12
                    // DOL changed by adding parameter t0, so John's comment applies inside function with t0 or time
//                    bool overlappedEddies = false;
                    if(!odtP.Llaminar)
                        LEA = sampleEddy(iStart, iEnd, t0, time, iStat);  // may reduce dtSmean; sets Pa; does eddy
                    
                    else LEA = false;
                    iEtrials++;

                    //---------- diffuse from t0 to time of implemented eddy (catch up)


                    if(LEA) {
                        // if the eddy was accepted ...
#ifdef NEWSTATS
                        if(breakReali){
                            breakUpStats[iReals-1][0] = time;
                            breakUpStats[iReals-1][1] = ed.eddySize;
                            // the first phase of the stats::phases array will be taken as the phase of the jet
                            iStart = mainLine->linePositionToIndex(ed.leftEdge,  true); // number of left face of the eddy region
                            iEnd   = mainLine->linePositionToIndex(ed.rightEdge, false); // number of right face of the eddy region
                            if (iStart == 0 || iEnd == mainLine->ngrd){
                                cout << "\n\nERROR:\n";
                                cout << "iStart = " << iStart << "   iEnd = " << iEnd << endl << endl;
                                exit(0);
                            }
                            //
                            // TODO:
                            // 
                            // if cases are currently wrong !!!
                            //
                            if (mainLine->phase[iStart-1] == odtStats.phases[0]
                               &&   mainLine->phase[iEnd] != odtStats.phases[0]){
                                // 2-phase eddy on the right edge of the core jet
                                while(mainLine->phase[iStart] == odtStats.phases[0])
                                    iStart++; // find first cell not belonging to the core jet
                                while(mainLine->phase[iStart] != odtStats.phases[0])
                                    iStart++; // find first cell of droplet
                                while(mainLine->phase[iEnd-1] != odtStats.phases[0])
                                    iEnd--;   // find last cell of droplet
                            }
                            else if (mainLine->phase[iStart-1] != odtStats.phases[0]
                                    &&   mainLine->phase[iEnd] == odtStats.phases[0]){
                                // 2-phase eddy on the left edge of the core jet
                                while(mainLine->phase[iEnd-1] == odtStats.phases[0])
                                    iEnd--;   // find last cell not belonging to the jet
                                while(mainLine->phase[iEnd-1] != odtStats.phases[0])
                                    iEnd--;   // find last cell of the droplet
                                while(mainLine->phase[iStart] != odtStats.phases[0])
                                    iStart++; // find first cell of the droplet
                            }
                            else{
                                // 2-phase eddy containing the whole core jet; final break up
                                while(mainLine->phase[iStart] != odtStats.phases[0])
                                    iStart++; // find first cell belonging to the first droplets
                                iEnd = iStart +1; // start surch of last cell of the first droplet
                                while(mainLine->phase[iEnd] == odtStats.phases[0])
                                    iEnd++;
                            }
                            breakUpStats[iReals-1][2] = mainLine->posf[iEnd] - mainLine->posf[iStart];
                            cout << "\n Falko breakUpStats = " << breakUpStats[iReals-1][0] << " " << breakUpStats[iReals-1][1] << " " << breakUpStats[iReals-1][2] << endl;
                            break; // while: advance to tMark
                        }
#endif
                        if(neddies % (odtP.modDisp*50) == 0) outputHeader();

                        neddies++;                             // increment total eddies implemented
                        if(neddies%odtP.modDisp==0) outputProgress(time, t0, neddies, iEtrials);

                        ss1.clear(); ss2.clear();

                        ss1 << neddies; ss1 >> s1; 
                        ss2 << iReals;  ss2 >> s2;
#ifdef NEWSTATSF
                        eddyfile << scientific << setprecision(8) << endl
                            << setw(17) << neddies
                            << setw(17) << time
                            << setw(17) << ed.leftEdge
                            << setw(17) << ed.rightEdge;
#endif
                        // Dump step 1: dump eddy information

                        // For a single realization and a positive modDump...
//                        if (odtP.nOdtReals == 1 && odtP.modDump > 0) {
//                            if (neddies % odtP.modDump == 0) {
//                                dumpLine (s1 + "_eddy.dat", gnufile, gnufilePart);     // comment
//                            }
//                        }

                        // For multiple realizations and a positive modDump...
//                        if (odtP.nOdtReals > 1 && odtP.modDump > 0) {
//                            if (neddies % odtP.modDump == 0) {
//                                dumpLine (s2 + "_" + s1 + "_eddy.dat", gnufile, gnufilePart); //comment
//                            }
//                        }
                        
                        // Dump step 2: dump adapted eddy information
#ifdef NEWSTATSF
                        if(!(odtP.Llem)) odtStats.BSetOld(*mainLine);
#endif
                        adaptEddyRegionOfMesh(iStart, iEnd, dtCUmax, time, tLastDA, cLastDA);   
#ifdef NEWSTATSF
                        if(!(odtP.Llem)){
                            odtStats.BChange(*mainLine, 6, iStat);
                            //odtStats.BChange(*mainLine, 4, iStat);
                        }
#endif
                        if(odtP.ItableLookup){
                            mainLine->etaTools->updateOdtLineVecs();
                        }
                        
                        // For a single realization and a positive modDump...
                        if (odtP.nOdtReals == 1 && odtP.modDump > 0) {
                            if (neddies % odtP.modDump == 0) {
                                dumpLine (s1 + "_adptEd.dat", gnufile, gnufilePart);
                            }
                        }

                        // For multiple realizations and a positive modDump..
                        if (odtP.nOdtReals > 1 && odtP.modDump > 0) {
                            if (neddies % odtP.modDump == 0) {
                                dumpLine (s2 + "_" + s1 + "_adptEd.dat", gnufile, gnufilePart);
                            }
                        }

                        // Dump step 3: dump diffused information

                        dtDiffuse = time - t0;

                        if(dtDiffuse > 0.0)  {           // if stmt only really needed for eddies
                                                         // right after the catchup block
#ifdef NEWSTATS
                            if(!(odtP.Llem)) odtStats.BSetOld(*mainLine);
#endif
                            posLower = mainLine->posf[iStart];   // to update iStart, iEnd below
                            posUpper = mainLine->posf[iEnd+1];
                            lineDiffuser->diffuseOdtLine(dtDiffuse, t0-timeOffset, odtStats, iStat);
                            iStart = mainLine->linePositionToIndex(posLower, true); 
                            iEnd   = mainLine->linePositionToIndex(posUpper, false); 

#ifdef NEWSTATSF
                            if(!(odtP.Llem)){
                                odtStats.BChange(*mainLine, 2, iStat);
                                odtStats.BChange(*mainLine, 4, iStat);
                                odtStats.BSetOld(*mainLine);
                            }
#endif

                            // For a single realization and a positive modDump...
//                            if (odtP.nOdtReals == 1 && odtP.modDump > 0) {
//                                if (neddies % odtP.modDump == 0) {
//                                    dumpLine (s1 + "_diffuse.dat", gnufile, gnufilePart); //comment
//                                }
//                            }

                            // For multiple realizations and a positive modDump...
//                            if (odtP.nOdtReals > 1 && odtP.modDump > 0) {
//                                if (neddies % odtP.modDump == 0) {
//                                    dumpLine(s2 + "_" + s1 + "_diffuse.dat", gnufile, gnufilePart); //comment
//                                }
//                            }

                            //iStart = mainLine.linePositionToIndex(mainLine.posf[iStart], true); 
                            //iEnd   = mainLine.linePositionToIndex(mainLine.posf[iEnd+1], false);

                            // Dump step 4: dump adapted diffused information
                            
                            adaptEddyRegionOfMesh(iStart, iEnd, dtCUmax, time, tLastDA, cLastDA);

                            if(odtP.ItableLookup){
                                mainLine->etaTools->updateOdtLineVecs();
                            }
                            
                            // For a single realization and a positive modDump...
//                            if (odtP.nOdtReals == 1 && odtP.modDump > 0) {
//                                if (neddies % odtP.modDump == 0) {
//                                    dumpLine (s1 + "_adptDiff.dat", gnufile, gnufilePart); //comment
//                                }
//                            }
                            
                            // For multiple realizations and a positive modDump...
//                            if (odtP.nOdtReals > 1 && odtP.modDump > 0) {
//                                if (neddies % odtP.modDump == 0) {
//                                    dumpLine (s2 + "_" + s1 + "_adptDiff.dat", gnufile, gnufilePart); //comment
//                                }
//                            }
                            
                            adaptAfterSufficientDiffTime(time, tLastDA, cLastDA, dtCUmax);

                            if(odtP.ItableLookup){
                               mainLine->etaTools->updateOdtLineVecs();
                            }
#ifdef NEWSTATSF
                            if(!(odtP.Llem)){
                                odtStats.BChange(*mainLine, 6, iStat);
                                //odtStats.BChange(*mainLine, 4, iStat);
                            }
#else
  #ifdef NEWSTATS
                            if(!(odtP.Llem)){
                                odtStats.BChange(*mainLine, 2, iStat);
                                odtStats.BChange(*mainLine, 4, iStat);
                            }
  #endif
#endif
                            t0 = time;
                        }

                        // Dump statistics for a single realization
//                        if(neddies%odtP.modDump==0) {
//                            dumpStats("means_" + s1 + ".dat", gnufileMeans, iStat);
//                        }

                        // Dump statistics for multiple realizations
                        if(neddies%odtP.modDump==0) {
//                            dumpStats("means_" + s2 + "_" + s1 + ".dat", gnufileMeans, iStat);
                        }
                        
                    }  

                    time += sampleDt();             // advance the time

                    //---------- raise dtSmean if needed

                    raiseDtSmean(iEtrials);         // may reset PaSum, nPaSum, dtSmean
                    
		  }
                }   //****************************** end while: advance to tMark
                
                if (breakReali)
                    break; // loop over sub intervals
                dtDiffuse = tMark-t0;
                if(dtDiffuse > 0.0)  {
                    *proc.ostrm << endl << "# diffusing to time tMark " << tMark; 
                    (*proc.ostrm).flush();

#ifdef NEWSTATS
                    if(!(odtP.Llem)) odtStats.BSetOld(*mainLine);
#endif
                    posLower = mainLine->posf[iStart];   // to update iStart, iEnd below
                    posUpper = mainLine->posf[iEnd+1];
                    lineDiffuser->diffuseOdtLine(dtDiffuse, t0-timeOffset, odtStats, iStat);
                    iStart = mainLine->linePositionToIndex(posLower, true); 
                    iEnd   = mainLine->linePositionToIndex(posUpper, false);
#ifdef NEWSTATSF
                    if(!(odtP.Llem)){
                        odtStats.BChange(*mainLine, 2, iStat);
                        odtStats.BChange(*mainLine, 4, iStat);
                        odtStats.BSetOld(*mainLine);
                    }
#endif
                    adaptAfterSufficientDiffTime(time, tLastDA, cLastDA, dtCUmax); //wait

                    if(odtP.ItableLookup){
                      mainLine->etaTools->updateOdtLineVecs();
                    }
                    
                    t0 = tMark;

#ifdef NEWSTATSF
                    if(!(odtP.Llem)){
                        odtStats.BChange(*mainLine, 6, iStat);
                        //odtStats.BChange(*mainLine, 4, iStat);
                    }
#else
  #ifdef NEWSTATS
                    if(!(odtP.Llem)){
                        odtStats.BChange(*mainLine, 2, iStat);
                        odtStats.BChange(*mainLine, 4, iStat);
                    }
  #endif
#endif

                }
                if(odtP.Llaminar){
                    ss1.clear(); 
                    ss1 << time; ss1 >> s1; 
                    dumpLine (s1 + ".dat", gnufile, gnufilePart);
                }
                

            }   //****************************** end loop over sub intervals
            
            if (breakReali)
                break; // loop over stat periode
            dumpStatsEndIstat(iStat);

            //ss1.clear(); ss1 << iStat; ss1 >> s1;
            //ss3.clear(); ss3 << iReals; ss3 >> s3;
            //mainLine->outputProperties("EndePeri_"+s3+"_"+s1+".dat");
            if(!(odtP.Llem))
                if (odtP.LperiRestart == iStat)
                    dumpStatsToRestart(); // writing new periRestart file

        }   //***************************** end loop over stat periods
        
        if( odtP.Iparticles && part->Lhistories ) {
            // part->outputHistories(proc.dataDir+"hist_part_" );
            part->outputHistories(proc.dataDir+"hist_part.dat" );
        }
        if (!breakReali){
            breakUpStats[iReals-1][0] = -1.0;
            breakUpStats[iReals-1][1] = -1.0;
            breakUpStats[iReals-1][2] = -1.0;
        }
#ifdef NEWSTATF
        eddyfile.close();
#endif
    }   //****************************** end loop over realizations

    gnufile.close();
    if(gnufilePart) gnufilePart.close();
//    gnufileMeans.close();
        
    *proc.ostrm << endl << "# number of Large eddy suppressions = " << iLES << endl;
#ifdef NEWSTATS
    *proc.ostrm << endl << "\n# break up statistik";
    *proc.ostrm << endl << "i           time            eddySize            dropletSize";
    *proc.ostrm << scientific; *proc.ostrm << setprecision(16);
    double meanSize = 0.0, meanTime = 0.0, meanTime2 = 0.0;
    double maxTime  = 0.0, minTime  = odtP.tEnd;
    int    iBreakUps = 0,  iNonBreakUps = 0;
    for(int iReals=0; iReals<odtP.nOdtReals; iReals++) {
        if (breakUpStats[iReals][0] != -1.0) {
            if (maxTime < breakUpStats[iReals][0]-0.02)
                maxTime = breakUpStats[iReals][0]-0.02;
            if (minTime > breakUpStats[iReals][0]-0.02)
                minTime = breakUpStats[iReals][0]-0.02;
            iBreakUps++;
            meanTime  += breakUpStats[iReals][0]-0.02;
            meanTime2 += (breakUpStats[iReals][0]-0.02)*(breakUpStats[iReals][0]-0.02);
            meanSize  += breakUpStats[iReals][2];
            *proc.ostrm << endl << iReals+1 << " " << " " << setw(25) << breakUpStats[iReals][0]-0.02
                << " " << setw(25) << breakUpStats[iReals][1] << " " << setw(25) << breakUpStats[iReals][2];
        }
        else
            iNonBreakUps++;
    }
    *proc.ostrm << endl << "mean break up time = " << setw(25) << meanTime/(1.0*iBreakUps);
    *proc.ostrm << endl << "rms break up time =  " << setw(25) << sqrt( (meanTime2/(1.0*iBreakUps)) - (meanTime/(1.0*iBreakUps)) * (meanTime/(1.0*iBreakUps)) );
    *proc.ostrm << endl << "min break up time =  " << setw(25) << minTime;
    *proc.ostrm << endl << "max break up time =  " << setw(25) << maxTime;
    *proc.ostrm << endl << "mean droplet size  = " << setw(25) << meanSize/(1.0*iBreakUps);
    *proc.ostrm << endl << "percentage of breaking realisations:     " << 100*iBreakUps/odtP.nOdtReals << "%";
    *proc.ostrm << endl << "percentage of non-breaking realisations: " << 100*iNonBreakUps/odtP.nOdtReals << "%";
#endif
}

///////////////////////////////////////////////////////////////////////////////
/**Computes dtCUmax (as the name suggests).  This variable is the time 
 * increment of eddy trial time advancement before we diffuse to catch
 * up to that time in the event of no eddy before that time.
 *
 * @return advancement catch up time increment.
 */

double odtSolver::computeDtCUmax() {

    double dxmin = 1.0E10;
    double d1;
    for(int i=1; i<mainLine->ngrdf; i++) {
        d1 = mainLine->posf[i]-mainLine->posf[i-1];
        if(d1 < dxmin)
            dxmin = d1;
    }
    if(!odtP.Lspatial){
        //doldb return 50.0*dxmin*dxmin/odtP.visc_0*odtP.rho_0;
        return odtP.tdfac * dxmin * dxmin / odtP.visc_0 * odtP.rho_0;
        //return 1.0 * dxmin * dxmin / odtP.visc_0 * odtP.rho_0;
        //return (val<0.05)?val:0.05;
    }
    else
        return 50.0 * odtP.tdfac * dxmin;

}

///////////////////////////////////////////////////////////////////////////////

/**dtSmean is computed, which is the mean eddy sample time.  The Poisson 
 * process draws dt's with this mean.
 *
 * @return dtSmean, the mean sample time increment.
 */

double odtSolver::computeDtSmean() {

    if(odtP.Llem)
        return 1.0/(odtP.lemRateParam * odtP.domainLength);
    
    else if(!odtP.Lspatial)
        return 0.1*odtP.Pav * odtP.domainLength * odtP.domainLength * odtP.rho_0 
                    / odtP.visc_0 / mainLine->ngrd / mainLine->ngrd / mainLine->ngrd;
    else 
        return 10.*odtP.Pav * odtP.domainLength / mainLine->ngrd / mainLine->ngrd;
   

}

///////////////////////////////////////////////////////////////////////////////

/**Sample an eddy size and position.  Create the eddy as an odtline.
 * Then triplet map the eddy, compute the eddy timescale, then the
 * acceptance probability.  Roll the dice and if you win (rand nu < prob)
 * then accept the eddy.  This means, apply velocity kernels, then
 * insert the eddy into the old line: delete the old eddy region and
 * insert the new.  If the domain is periodic and the eddy splits the
 * domain, the eddy is tacked entirely on the end of the domain, and 
 * the starting position is no longer zero (the whole line moves to the right).
 * A flag is returned if the eddy was accepted/implemented
 *
 * @param iStart \inout index range of the eddy.
 * @param iEnd   \inout index range of the eddy.
 * @param t0     \input eddy implementation time.  Added 7/31/12 DOL, based on JCH's mod
 * @param time   \input current sampling time.
 * @param iStat  \input number of iStat period. Needed for stats-object update
 * @return true if the sampled eddy was implemented.
 */

bool odtSolver::sampleEddy(int &iStart, int &iEnd, double t0, double time, const int &iStat) {

    bool L_eddyImplemented = false;
    
    ed.LmultiPhaseEddy = false;
    ed.sampleEddySize(odtP);
    if(!testLES_fracDomain(ed.eddySize) && !odtP.Llem) {
        iLES++;
        ed.Pa=0.0;
        return false;
    }
    // ed.sampleEddyPosition(odtP, *mainLine, rndmgn);
    ed.sampleEddyPosition(odtP, *mainLine);

    //---------- extract the eddy segment from mainLine into eddyLine

    if(!ed.LperiodicEddy) {
        iStart = mainLine->linePositionToIndex(ed.leftEdge,  true);
        iEnd   = mainLine->linePositionToIndex(ed.rightEdge, false);
        if( iEnd - iStart < odtP.eddyMinCells && !odtP.Llem) { //only for odt because numerical shear can be >> physical shear
            ed.Pa = 0.0;
            return false;
        }
        if(!testThirds() && !odtP.Llem) {                     // large eddy supression test
            iLES++;
            ed.Pa = 0.0;
            return false;
        }
    }
    else {
        iStart = mainLine->linePositionToIndex(ed.leftEdge,  true);
        iEnd   = mainLine->linePositionToIndex(ed.rightEdge, false);
        if( (mainLine->ngrd-iStart)+(iEnd)+1 < odtP.eddyMinCells  && !odtP.Llem) { 
            ed.Pa = 0.0;
            return false;
        }
        // dolcheck large eddy suppression test here
    }

#ifndef COMPSCI
    odtline *eddyLine;
    eddyLine = new odtline(mainLine, iStart, iEnd, ed.LperiodicEddy, false);
#else
    if (!eddyLine) {
		eddyLine = new odtline(mainLine, iStart, iEnd, ed.LperiodicEddy, false);
	} else {
		eddyLine->initFromCopy(mainLine, iStart, iEnd, ed.LperiodicEddy, false);
	}
#endif

    //---------- invoke the eddy

    ed.tripMap(*eddyLine);                       // apply trip map to eddyLine
    
    if(!odtP.Llem){
        ed.fillKernel(*eddyLine);                    

        if(odtP.LconstProp)  {
            if(!ed.eddyTauCP(*eddyLine, odtP, -1)) { 
                ed.Pa = 0.0;
#ifndef COMPSCI
                delete eddyLine;
#endif
                return false; 
            }              
        }
        else {
            if(!ed.eddyTau(*eddyLine, odtP, -1)) {       
                ed.Pa = 0.0;
#ifndef COMPSCI
                delete eddyLine;
#endif
                return false;
            }               
        }

        double timeSpaceProgress = (odtP.Lspatial ? ed.eddySize : 1.0/ed.invTauEddy);
        //DOL if(!testLES_elapsedTime(t0, timeSpaceProgress)) {
        if(!testLES_elapsedTime(time, timeSpaceProgress)) {
            iLES++;
            ed.Pa = 0.0;
#ifndef COMPSCI
            delete eddyLine;
#endif
            return false;
        }

        if(!testLES_integralLength(time, ed.eddySize)) { //GSDB
            iLES++;
            ed.Pa = 0.0;
#ifndef COMPSCI
            delete eddyLine;
#endif
            return false;
        }

        ed.computeEddyAcceptanceProb(odtP, dtSmean);

        //---------- lower dtSample if Pa too high; changes Pa, dtSmean

        lowerDtSmean(time);    
    }
    //---------- apply the eddy if accepted

    double rnd01 = rndmgn.getRand();

    if(rnd01 <= ed.Pa || odtP.Llem) {
        if (ed.LmultiPhaseEddy)
            breakReali = true;
        
        if(ed.LperiodicEddy)
            *proc.ostrm << endl << "#   periodic eddy ";

        if(Ltune) tuneMinGradFrac(iStart, iEnd); //tune

        //--------- do particles
        if (odtP.Iparticles)
            ed.tripletMapParticles(odtP, *mainLine, *eddyLine, *part, t0, time);
        
        ed.applyVelocityKernels(*eddyLine, odtP);
                
        //--------- 

#ifdef NEWSTATS
        if(!(odtP.Llem)) odtStats.BSetOld(*mainLine);
#endif

        mainLine->insertEddy(*eddyLine, iStart, iEnd, 
                            ed.leftEdge, ed.rightEdge, ed.LperiodicEddy);

        if(odtP.Iparticles) {
            if(part->PeddyType == 2 || part->PeddyType == 3) {  // ----- TypeC or TypeIC
                randomGenerator     rr(-1);
                double r = rr.getRand();
                int sign = 1; 
                if (r >= 0 && r < 0.5) sign = -1;
                double Pvel = sign * 2 * ed.eddySize * 0.19245 / part->ParamEddylife[0] * ed.invTauEddy; //GYSun

                lineDiffuser->updateActiveEddy(ed.leftEdge, ed.rightEdge, ed.eddySize, ed.eddyLife[0], 
                        ed.uEddyAvg, ed.vEddyAvg, ed.wEddyAvg, time+ed.eddyLife[0], Pvel, ed.iPartInEddy,
                        ed.iPartInEddyIC, ed.xPartPos, ed.zPartPos//GYSun
                        );
            }
        }
        
#ifdef NEWSTATS
        if(!(odtP.Llem)){
            odtStats.BChange(*mainLine, 0, iStat);
            odtStats.BChange(*mainLine, 4, iStat);
        }
        
        // adding eddy to eddy size statistic
        double dist = min(ed.leftEdge,mainLine->Ldomain-ed.rightEdge);
        if (dist >= ed.eddySize){
            odtStats.cstat.at(0).at(4).at(iStat-1).at(int((ed.eddySize/mainLine->Ldomain)*(odtStats.ngrd/3.0))+1)++;
        }
        else{ // near wall eddy
            odtStats.cstat.at(1).at(4).at(iStat-1).at(int((ed.eddySize/mainLine->Ldomain)*(odtStats.ngrd/3.0))+1)++;
        }
  #ifdef NEWSTATSF
        eddyRegion.resize(mainLine->ngrd);
        std::fill(eddyRegion.begin(), eddyRegion.end(), 0.0);
        std::fill(eddyRegion.begin()+iStart, eddyRegion.begin()+iEnd, 1.0);
        odtStats.odtGrd2statGrd(mainLine->posf, eddyRegion, 0.0);
        for (int i = 0; i < (int)odtStats.vTrans.size(); i++)
            odtStats.cstat[2][4][iStat-1][i] += odtStats.vTrans[i];
  #endif
#endif

        L_eddyImplemented = true;
    }
    
#ifndef COMPSCI
    delete eddyLine;
#endif
    return L_eddyImplemented;

}

///////////////////////////////////////////////////////////////////////////////

/**Reduce dtSmean if it is resulting in too high an acceptance probability. 
 * @param time \input current run time
 */

void odtSolver::lowerDtSmean(double time) {

    if(odtP.Llem) return;
    
    if(ed.Pa == 0.0) return;

    if(ed.Pa > odtP.Pmax) {
        dtSmean *= odtP.Pmax/ed.Pa;
        ed.Pa = odtP.Pmax;
        *proc.ostrm << endl << "#   reducing dtSmean to " << dtSmean;
        *proc.ostrm << endl << "#   ed.Pa, odtP.Pmax = " << ed.Pa << " " << odtP.Pmax;
        *proc.ostrm << endl << "#   time " << time;
    }

    //---------- update properties used to raise dtSample

    PaSum += ed.Pa;
    nPaSum++;

    PaSumC += ed.Pa;
    nPaSumC++;

}

///////////////////////////////////////////////////////////////////////////////

/**Routine by ARK
 * Called if Ltune is true
 * Periodically reset odtP.smallGradFrac based on number of cells in the smallest
 * eddies.  smallGradFrac scales as \fun{ Re^{-1/4} }, so it will change with different \c Re
 * flows, and within a given temporally evolving problem.
 * If the smallest accepted eddies are adequately resolved, but no over-resolved, then
 * resolution should be near-optimal.
 * Set target so its not too close or too far from eddyMinCells.
 * Currently tuneWait and target are hardcoded here, but could be read in from
 * odtParam.inp
 *
 * @param iStart \input beginning index of region of interest.
 * @param iEnd   \input beginning index of region of interest.
 */

void odtSolver::tuneMinGradFrac(int iStart, int iEnd) {


    // The following are hard-wired but should instead be read from input

    const int tuneWait = 1000;   // number of accepted eddies before doing adjustment
    const double target = 10.;   // target value of fewestCells

    int cells = iEnd - iStart;
    if(ed.LperiodicEddy) cells += mainLine->ngrd + 1;
    fewestCells = min(fewestCells, cells);

    if((neddies+1) % tuneWait != 0)         // neddies has not yet been incremented
        return;                             // only do the adjustment once in a while

    if(target < odtP.eddyMinCells) {
        *proc.ostrm << "\nError in tunMinGradFrac. Ltune is on with target < eddyMinCells" << endl;
        exit(0);
    }

    // Coarsen mesh if too many cells, refine mesh if too few cells (relative to target)
    odtP.smallGradFrac *= static_cast<double>(fewestCells) / target;

    *proc.ostrm << endl << "neddies, fewestCells, smallGradFrac: " << neddies << " " <<
        fewestCells << " " << odtP.smallGradFrac << endl;

    fewestCells = mainLine->ngrd;   // Initialize fewestCells for the next adjustment cycle

}

///////////////////////////////////////////////////////////////////////////////

/**Every once in a while (nDtSmeanWait) test the mean acceptance probability.
 * Increase dtSmean if its too small.
 *
 * @param iEtrials \input number of eddy trials so far.
 */

void odtSolver::raiseDtSmean(int &iEtrials) {  

    if(odtP.Llem) 
	return;
    
    if(iEtrials % odtP.nDtSmeanWait != 0)
        return;                             // only do this once in a while

    if(nPaSum > 0)
        PaSum /= nPaSum;                    // PaSum is now average Pa of eddies

    if(PaSum < odtP.Pav) {
       if(PaSum < odtP.Pav/odtP.dtfac)
           dtSmean *= odtP.dtfac;          // increase by at most dtfac
       else
           dtSmean *= odtP.Pav/PaSum;      // increase dtSmean to target value

       *proc.ostrm << endl << "#  increasing dtSmean to " << dtSmean 
           << " (neddies = " << neddies << "  eddy trails = " << iEtrials << ")";
    } 

    PaSum  = 0.0;                          // reset values
    nPaSum = 0;
    if (iEtrials > 10000*odtP.nDtSmeanWait) {
        *proc.ostrm << endl << "#  reset iEtrials, PaSumC, nPaSumC after "
            << 10000*odtP.nDtSmeanWait << " eddy trials.";
        iEtrials = 0;
        PaSumC   = 0.0;
        nPaSumC  = 0;
    }

}

///////////////////////////////////////////////////////////////////////////////

/**Output title of properties displayed to screen. */

void odtSolver::outputHeader() {

    *proc.ostrm << endl << "#--------------------------------------------------"
                 << "--------------------------------------------------------------------";
    *proc.ostrm << endl;
    *proc.ostrm << setw(5) << "# EE,"
        << setw(12) << "time,"
        << setw(12) << "t-t0,"
        << setw(10) << "nEtry,"
        << setw(6)  << "ngrd,"
        << setw(12) << "edSize,"
        << setw(12) << "edPos,"
        << setw(12) << "edPa,"
        << setw(12) << "nEposs,"
        << setw(12) << "PaAvg,"
        << setw(12) << "invTauEddy"
        ;

    *proc.ostrm << endl << "#--------------------------------------------------"
                 << "--------------------------------------------------------------------";
}

///////////////////////////////////////////////////////////////////////////////

/**Outputs the data corresponding to outputHeader.
 * After a given number of accepted eddies, output this info.
 *
 * @param time \input current eddy sample time.
 * @param t0 \input current diffused time.
 * @param neddies \input number of accepted eddies.
 * @param iEtrials \input number of attempted eddies.
 */

void odtSolver::outputProgress(double &time, double &t0, int &neddies, int &iEtrials) {


    double dmb = 0.5*(ed.leftEdge+ed.rightEdge);
    dmb = (dmb > mainLine->posf[mainLine->ngrd]) ? dmb-mainLine->Ldomain : dmb;

    *proc.ostrm << scientific << setprecision(3) << endl;
    *proc.ostrm << setw(5) << neddies   //  1: EE
        << setw(12) << time             //  2: time
        << setw(12) << time-t0          //  3: t-t0
        << setw(10) << iEtrials         //  4: nEtry
        << setw(6)  << mainLine->ngrd   //  5: ngrd
        << setw(12) << ed.eddySize      //  6: edSize
        << setw(12) << dmb              //  7: edPos
        << setw(12) << ed.Pa            //  8: edPa
        << setw(12) << nPaSumC          //  9: nEposs
        << setw(12) << PaSumC/nPaSumC   // 10: PaAvg
        << setw(12) << ed.invTauEddy       // 11: invTauEddy
        ;
    (*proc.ostrm).flush();
}

///////////////////////////////////////////////////////////////////////////////

/**Used for selective mesh adaptation after a sufficient time.
 * Write the current value of time into the cells of lastDA that overlap or are 
 * contained in the mesh adaption interval.  If cLastDA is one of those cells, scan the lastDA 
 * array to find the smallest time value in the array.  Set tLastDA equal to that 
 * smallest time value and set cLastDA equal to the index of the cell containing 
 * that smallest value.
 *
 * @param time \input current time.
 * @param tLastDA \inout time of last diffusion advancement.
 * @param cLastDA \inout cell index.
 * @param iStart \input starting cell index to check.
 * @param iEnd \input ending cell index to check.
 */

void odtSolver::updateDA(const double &time, double &tLastDA, int &cLastDA, 
			             int iStart, int iEnd) {

    //----------- find the range of affected lastDA cells

    int leftCell  = static_cast<int>(odtP.sLastDA * mainLine->posf[iStart] / odtP.domainLength);
    if(leftCell==odtP.sLastDA) leftCell--;
    int rightCell = static_cast<int>(odtP.sLastDA * mainLine->posf[iEnd+1] / odtP.domainLength);
    if(rightCell == odtP.sLastDA) rightCell--;

    if(leftCell > rightCell && !odtP.Lperiodic) 
        return;            // no affected cells of LastDA array

    //---------- set the lastDA values in this range equal to the current time
    
    if(leftCell > rightCell) {                    // region like | * | * |   |   | * | * |
        for(int i=leftCell; i<odtP.sLastDA; i++)
            lastDA[i] = time;
        for(int i=0; i<=rightCell; i++)
            lastDA[i] = time;
        if( !(cLastDA >= leftCell) && !(cLastDA <= rightCell) )
            return;                               // oldest cell is unaffected
    }
    else {                                        // region like |   | * | * | * | * |   |
        for(int i=leftCell; i<=rightCell; i++) 
            lastDA[i]=time;
        if(!( (cLastDA >= leftCell) && (cLastDA <= rightCell) )) 
            return;                               // oldest cell is unaffected
    }

    //---------- find oldest cell: cLastDA, and earliest last adapt time

    cLastDA = 0;                                  // initialize cLastDA
    for(int i=1; i<odtP.sLastDA; i++) 
        if(lastDA[i] < lastDA[cLastDA]) 
            cLastDA = i;                          // find oldest cell
    tLastDA = lastDA[cLastDA];                    // earliest last-adapt time

}

///////////////////////////////////////////////////////////////////////////////

/**Adapt regions of the mesh depending on age since last adaption.
 * Start at cell index cLastDA in the "lastDA" domain (grid) and search right and left
 * for old cells.  Update these cells (tag them with current time), then adapt them.
 * Then cLastDA (from Update) is the next oldest cell not just adapted.  Search left
 * and right for a new region, and so on until nothing on the grid is too old.
 *
 * @param time \input current time.
 * @param tLastDA \inout time of last diffusion advancement.
 * @param cLastDA \inout cell index.
 * @param dtCUmax \inout the time increment of eddy trial time advancement before we diffuse
 */

bool odtSolver::adaptAfterSufficientDiffTime (const double &time, 
                                              double       &tLastDA,
                                              int          &cLastDA,            
                                              double       &dtCUmax) {

    if(time - tLastDA < odtP.DAtimeFac * dtCUmax) 
        return false;    

    int    iStart, iEnd;

    double youngsters = odtP.DAtimeFac*dtCUmax;

    while(time - tLastDA >= youngsters) {  

        //---------- initialize searches

        int leftSide  = cLastDA;
        int rightSide = (cLastDA == odtP.sLastDA-1) ? cLastDA : cLastDA+1;

        //---------- search for leftSide

        while(time - lastDA[leftSide] >= 0.5 * youngsters){          
            if(leftSide == 0) 
                break;
            leftSide--;
        }

        //---------- correct for possible overshoot

        if(time - lastDA[leftSide] < 0.5 * youngsters)
            leftSide++; 

        //---------- search for rightSide

        while(time - lastDA[rightSide] >= 0.5 * youngsters){
            if(rightSide == odtP.sLastDA-1) 
                break;
            rightSide++;
        }

        //---------- correct for possible overshoot

        if(time - lastDA[rightSide] < 0.5 * youngsters)
            rightSide--;

        //---------- find the corresponding range of indices of the posf array

        iStart = mainLine->linePositionToIndex(leftSide *odtP.domainLength/
                                              odtP.sLastDA, true);
        iEnd   = mainLine->linePositionToIndex(rightSide*odtP.domainLength/
                                              odtP.sLastDA, false);

        updateDA(time, tLastDA, cLastDA, iStart, iEnd);

        //---------- adapt the mainLine

        if(mainLine->LhasRxn) mainLine->setTempVec(); 
        
        else if(odtP.ItableLookup) mainLine->etaTools->updateOdtLineVecs();
        
        mainLine->meshAdapter.adaptGrid(iStart, iEnd);
        
        if (odtP.ItableLookup) {
            mainLine->etaTools->updateOdtLineVecs();
        }
        if (odtP.Iparticles) {
            mainLine->setVoidFrac(this->part);
        }

    }

    //---------- update dtCUmax 

    dtCUmax = computeDtCUmax();

    return true;                               
}

///////////////////////////////////////////////////////////////////////////////

/** @param iStart \inout starting index of region to adapt.
 *  @param iEnd \inout ending index of region to adapt.
 *  @param dtCUmax \inout the time increment of eddy trial time advancement before we diffuse
 *  @param time \input current time.
 *  @param tLastDA \inout time of last diffusion advancement.
 *  @param cLastDA \inout cell index.
 */

void odtSolver::adaptEddyRegionOfMesh(int &iStart, int &iEnd, double &dtCUmax,
                                      const double &time, 
                                      double &tLastDA, int &cLastDA) {  

    if(iStart > 0)             iStart--;
    if(iEnd < mainLine->ngrd-1) iEnd++;

    double posLower = mainLine->posf[iStart];   // to update iStart, iEnd below
    double posUpper = mainLine->posf[iEnd+1];

    updateDA(time, tLastDA, cLastDA, iStart, iEnd);  

    if(mainLine->LhasRxn) mainLine->setTempVec(); 
    
    else if(odtP.ItableLookup) mainLine->etaTools->updateOdtLineVecs();
    //--------------------------------------------------

    
    mainLine->meshAdapter.adaptGrid(iStart, iEnd);
    
    if (odtP.ItableLookup) {
        mainLine->etaTools->updateOdtLineVecs();
    }
    if (odtP.Iparticles) {
        mainLine->setVoidFrac(this->part);
    }

    //---------- update dtCUmax 

    dtCUmax = computeDtCUmax();

    //---------- update iStart, iEnd

    iStart = mainLine->linePositionToIndex(posLower, true); 
    iEnd   = mainLine->linePositionToIndex(posUpper, false); 

}

///////////////////////////////////////////////////////////////////////////////

/**Sample the eddy trial time step with mean dtSmean. */

double odtSolver::sampleDt() {

    return -dtSmean*log( max(1.0e-14, rndmgn.getRand()) );
    
}

///////////////////////////////////////////////////////////////////////////////

/**Output a data file for the mainLine object.
 * Note, odtline::outputProperties is more general, and is read with restarts,
 * while outputProperties2 is nice for viewing (sometimes) in that
 * it plots the tophat profiles in cells that are assumed in mesh adaption 
 * and merge/split operations.
 *
 * @param fname       \input name of file to dump.
 * @param gnufile     \input name of gnuplot file with list of dumped odtline files.
 * @param gnufilePart \input name of gnuplot file where the particles are dumped.
 */

void odtSolver::dumpLine(string fname,
                         ofstream &gnufile,
                         ofstream &gnufilePart) {
    mainLine->outputProperties(proc.dataDir+"odt_"+fname);                  
    if(odtP.Lrxn)
        gnufile << "plot '" << "odt_"+fname << "' us 1:8; pause -1;" << endl; 
    else
        gnufile << "plot '" << "odt_"+fname << "' us 1:6; pause -1;" << endl; 

    //--------------- do particles

    if(odtP.Iparticles){
        part->outputProperties(proc.dataDir+"part_"+fname);
        gnufilePart << "plot '" << "part_"+fname << "' us 2:6 with impulses;pause-1;" << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////

/**Output a data file for the statistics object.  Similar to that for an #anyline.
 *
 * @param fname   \input name of statistics file to dump.
 * @param gnufile \input name of gnuplot file with list of dumpted stats files.
 */

void odtSolver::dumpStats(string fname, ofstream &gnufile, const int iStat) {

    odtStats.outputProperties(proc.dataDir+fname, iStat);                  
    gnufile << "plot '" << fname << "' us 1:3; pause -1;" << endl; 
}

///////////////////////////////////////////////////////////////////////////////

/**Output a data file for the statistics object at the end of a stats gathering period.
 *
 * @param iStat \input which statistic interval is it: goes into the file name.
 */

void odtSolver::dumpStatsEndIstat(const int &iStat) {

    if(odtP.Llem)
      return;
    stringstream ss1;
    string s1;
    *proc.ostrm << endl << "# Writing Stats for iStat interval " << iStat << endl;
    ss1.clear(); ss1 << iStat; ss1 >> s1;  
    odtStats.outputProperties(proc.dataDir+"stats_interval_" + s1 + ".stat", iStat);
    odtStats.initStats(mainLine);
}

///////////////////////////////////////////////////////////////////////////////
/** Overwrite the restart file for this processor. Though, you get a
 *  circulating condition.
 */
void odtSolver::dumpStatsToRestart(void){
    stringstream  ss1;
    string        s1;
    *proc.ostrm << endl << "#***************** dump stats to restart";
    *proc.ostrm << endl << "# Writing new peri restart file." << endl;
    ss1.clear(); ss1 << proc.myid; ss1 >> s1;
    mainLine->outputProperties("../input/restart_odtl_" + s1 + ".dat");
}


///////////////////////////////////////////////////////////////////////////////

/**Read a restart file to set the mainLine properties. Also update the input 
 * file properties for consistency.  That is, When restart, check/fix 
 * inconsistencies between the restart file and input file
 *
 * @param odtlRstFile \input name of odtline restart file.
 * @param partRstFile \input name of odtline restart file for particles.
 */

// void odtSolver::setRestartState(string odtlRstFile, string partRstFile) {
void odtSolver::setRestartState(string odtlRstFile) {

    if(odtP.Lspatial) {
        mainLine->readProperties(odtlRstFile);
        if(odtP.Lrxn) {
            mainLine->setTempVec();
            mainLine->setRhoVec();
            mainLine->setViscosity();
            mainLine->setMixfVec();
            mainLine->setChiVec();
        }
    }
    else {
        mainLine->readProperties(odtlRstFile);
        mainLine->setTempVec();
        mainLine->setRhoVec();
        mainLine->setViscosity();
        mainLine->setMixfVec();
        mainLine->setChiVec();
    }
 
    if(odtP.LconstProp) {     // for restarting with constant properties

        bool LredoRho = false;
        bool LredoMu  = false;

        for(int i=0; i<mainLine->ngrd; i++)
            if(mainLine->molec[i] != odtP.visc_0) {
                *proc.ostrm << "\n\n************* WARNING *************";
                *proc.ostrm <<   "\n************* LconstProp but visc_0 != molec from restart";
                *proc.ostrm <<   "\n************* Changing molec array to equal visc_0" << endl;
                LredoMu = true;
                break;
            }

        for(int i=0; i<mainLine->ngrd; i++)
            if(mainLine->rho[i] != odtP.rho_0) {
                *proc.ostrm << "\n\n************* WARNING *************";
                *proc.ostrm <<   "\n************* LconstProp but rho_0 != rho from restart";
                *proc.ostrm <<   "\n************* Changing rho array to equal rho_0" << endl;
                LredoRho = true;
                break;
            }

        if(LredoRho) 
            for(int i=0; i<mainLine->ngrd; i++)
                mainLine->rho[i] = odtP.rho_0;
        if(LredoMu) 
            for(int i=0; i<mainLine->ngrd; i++)
                mainLine->molec[i] = odtP.visc_0;

    }

    else {                    // for restarting with variable properties

        *proc.ostrm << "\n\n************ WARNING *************";
        *proc.ostrm <<   "\n************ Restarting variable prop case; using visc_0, rho_0";
        *proc.ostrm <<   "\n************    in input file for sampling operations" << endl;

    }

//    if(odtP.Iparticles)
//        part->readProperties(partRstFile);

//    if(!odtP.Llem){
//        if(odtP.LperiRestart > 0){
//            *proc.ostrm << "\n\n#***********************************************************";
//            *proc.ostrm << "\n#***************** odtP.LperiRestart = " << odtP.LperiRestart;
//            dumpStatsToRestart(); // writing initial periRestart file
//            *proc.ostrm << "#***********************************************************\n\n";
//        }
//    }

}

///////////////////////////////////////////////////////////////////////////////
/** Apply the large-eddy suppression test. 
 */
 
bool odtSolver::testThirds() {

	if(odtP.LES_type != 1)  
		return true;

#ifndef COMPSCI
	 eddy ed3;
#endif
	 ed3.eddySize = ed.eddySize/3.0;
	 ed3.LperiodicEddy=false;
 
     double leftThird = ed.leftEdge;
     double rightThird = leftThird + ed.eddySize/3.0;
 
     for(int third=0; third<3; third++) {

		 ed3.leftEdge  = leftThird;
		 ed3.rightEdge = rightThird;
 
         int iStart = mainLine->linePositionToIndex(leftThird,  true);
         int iEnd   = mainLine->linePositionToIndex(rightThird, false);
 
         odtline eddyLine(mainLine, iStart, iEnd);

         //---------- invoke the eddy
 
         ed3.tripMap(eddyLine);                       // apply trip map to eddyLine
 
         ed3.fillKernel(eddyLine);
 
         if(odtP.LconstProp)
             if(!ed3.eddyTauCP(eddyLine, odtP, odtP.Z_third)) 
                 return false;
         else
             if(!ed3.eddyTau(eddyLine, odtP, odtP.Z_third)) 
                 return false; 
 
         leftThird = rightThird;
         rightThird = leftThird + ed.eddySize/3.0;
 
     }
 
     return true; // Eddy is allowed (not suppressed).
 
 }

///////////////////////////////////////////////////////////////////////////////
/** Apply a large eddy suppression test based on elapsed time (Echekki 2001) 
 *  @param time \input current time.
 *  @param tauEddy \input eddy timescale, or in spatial cases, the eddy size
 */

bool odtSolver::testLES_elapsedTime(double time, double tauEddy) {

//      if(time > 0.003)
//          return false;
//
//      return true;

    if(odtP.LES_type != 2)   
        return true;               

    if( time < odtP.Z_third * tauEddy )
        return false;
    
    return true;

}

///////////////////////////////////////////////////////////////////////////////
/** Large eddy suppression test on fraction of domain size
 *  @param eSize \input eddy size
 */
bool odtSolver::testLES_fracDomain(double eSize) {
    if(odtP.LES_type != 3)   
        return true;               
    if(eSize/odtP.domainLength > odtP.Z_third)
        return false;
    return true;

}

///////////////////////////////////////////////////////////////////////////////
/** Apply a large eddy suppression test based on integral length scale 
 *  integral length scale L = L0 * (t/t0)^(1-n/2)
 *  t is elapsed time;
 *  n is between 1.15 and 1.45, usually 1.3
 *  @param time \input current time.
 *  @param eSize \input eddy size
 *  Guangyuan Sun 06/2013
 */

bool odtSolver::testLES_integralLength(double time, double eSize) {
    
    if(odtP.LES_type !=4)
        return true;
//    double n = 1.3;
//    double n = odtP.Z_third;
    double n = 1.1;
    double t0 = 0.15899;
//     double uprime = 0.47;
//     double integralLength = t0 * sqrt(0.5)*uprime/n*pow(time/t0, 1-0.5*n);
    
    double L0 = 0.028323;
    double integralLength = odtP.Z_third * L0 * pow(time/t0, 1-0.5*n);

    if(eSize > integralLength)
        return false;
    return true;

}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


/** This function does the same as setRestartState except skiping the writing
 *  of the processor bounded restart file.
 *
 * @param odtlRstFile input: name of odtline restart file.
 * @param partRstFile input: name of odtline restart file.
 */
void odtSolver::setPeriRestartState(std::string odtlRstFile, std::string partRstFile){
    int save = odtP.LperiRestart;
    odtP.LperiRestart = 0;
    // setRestartState(odtlRstFile, partRstFile);
    setRestartState(odtlRstFile);
    odtP.LperiRestart = save;
}

///////////////////////////////////////////////////////////////////////////////
/** This function reads a given odtChange_*.inp file located in the input folder
 *  before the next period starts. The * is a place holder for the periode which
 *  would be changed. This file could change parameters and can increase the 
 *  domain size for the remainder of the simulation.
 *  NEEDS FLAG NEWSTATS
 *  
 *  @Param iStat    The periode for which the change file will be read
 */
#ifdef NEWSTATS
void odtSolver::ReadChangeFile(const int iStat){
    
    ifstream     ifile;
    string       input_s;
    stringstream input_ss;
    double       Ldomain = odtStats.Ldomain;
    
    input_ss.clear();
    input_ss << "../input/odtChange_" << iStat << ".inp";
    input_ss >> input_s;
    
    *proc.ostrm << endl << "# input_s = " << input_s;
    
    ifile.open(input_s.c_str(), ifstream::in);
    if (!ifile.is_open()){
        *proc.ostrm << endl << "# WARNING: File " << input_ss.str() << " not opened!" << endl;
        *proc.ostrm << endl << "***********************************************************";
        return;
    }
    
    while (!ifile.eof()){
        
        getline(ifile, input_s);
        // options for comment lines
        if (input_s.size() == 0){
            //*proc.ostrm << endl << "#--- Line is empty !!!!";
            continue;
        }else if (input_s.compare(0,2,"--") == 0){
            //*proc.ostrm << endl << "# comment line with '--'";
            continue;
        }else if (input_s.compare(0,2,"//") == 0){
            //*proc.ostrm << endl << "# comment line with //";
            continue;
        }else if (input_s.compare(0,3,"/**") == 0){
            //*proc.ostrm << endl << "# comment block with '/**'";
            size_t found;
            found = input_s.find("*/");
            while (!ifile.eof() && found == string::npos){
                getline(ifile, input_s);
                //*proc.ostrm << endl << "# next line: " << input_s;
                found = input_s.find("*/");
            }
            *proc.ostrm << endl << "# found end of comment block";
            continue;
        } // options for increasing domain
        else if (string::npos != input_s.find("IncreaseDomainLower")){
            *proc.ostrm << endl << "# found order to increase domain at lower bc";
            IncreaseDomain(&ifile, 'L');
        }else if (string::npos != input_s.find("IncreaseDomainUpper")){
            *proc.ostrm << endl << "# found order to increase domain at upper bc";
            IncreaseDomain(&ifile, 'U');
        }else if (string::npos != input_s.find("IncreaseDomainCenter")){
            *proc.ostrm << endl << "# found order to increase domain half at both bcs";
            IncreaseDomain(&ifile, 'C');
        } // options to change parameters
        else if (string::npos != input_s.find("ChangeParameter")){
            *proc.ostrm << endl << "# found order to change several parameter";
            ChangeParameter(&ifile);
        }
        else if (string::npos != input_s.find("SetUMean")){
            *proc.ostrm << endl << "# found order to set the mean of the u velocity";
            SetUMean(&ifile);
        } // Line could not be read
        else{
            *proc.ostrm << endl << "# ERROR: Line in ../input/odtChange_" 
            << iStat << ".inp could not be read";
            *proc.ostrm << endl << "Line: >" << input_s << "<";
            exit(0);
        }
        
    }
    *proc.ostrm << endl << "# Reading change file has finished";
    if(Ldomain != mainLine->Ldomain && mainLine->Ldomain != odtStats.max_Ldomain){
        *proc.ostrm << endl << "# Creating a new odtStats variable keeping the old information";
        *proc.ostrm << endl << "# LDomain    = " << Ldomain;
        *proc.ostrm << endl << "# LDomain_ML = " << mainLine->Ldomain;
        stats newStats(&odtStats, iStat, mainLine->Ldomain);
        odtStats = newStats;
        odtStats.initStats(mainLine);
    }
    else{
        *proc.ostrm << endl << "# Reinitialize odtStats";
        odtStats.initStats(mainLine);
    }
    mainLine->meshAdapter.adaptGrid(0,mainLine->ngrd-1);
    *proc.ostrm << endl << "# Reading change file has finished";
    *proc.ostrm << endl << "#***********************************************************";
}

///////////////////////////////////////////////////////////////////////////////
/** This function increases the domain of the ODT-line at the lower and/or the
 *  upper boundary.
 *  
 *  @Param  *ifile  the opened input file
 *  @Param  LUC     L expansion at lower boundary, U expansion at upper boundary
 *                  C expansion with old domain centered
 */
void odtSolver::IncreaseDomain(ifstream *ifile, char LUC){
    
    *proc.ostrm << endl << "#--- start increase domain";
    
    if (mainLine->LhasRxn || mainLine->LhasMom || mainLine->LhasEta){
        *proc.ostrm << endl << "ERROR: the current code of "
        "odtSolver::IncreaseDomain is not programmed for reacting flows, "
        "moments, and ETA. Please update this function to work fine." << endl;
        exit(0);
    }
    
    stringstream input_ss;
    string  input_s;
    string  input;
    int     phase  = 0;
    double  rho    = 0.0;
    double  visc   = 0.0;
    double  lambda = 0.0;
    double  u      = 0.0;
    double  v      = 0.0;
    double  w      = 0.0;
    input_ss.clear();
    
    // extracting delz
    getline(*ifile, input_s); // this has to be the line with "delz"
    
    *proc.ostrm << endl << "#--- extracting delz for increase";
    *proc.ostrm << endl << "#--- got line: " << input_s;
    
    if (string::npos == input_s.find("delz")){
        *proc.ostrm << endl << "ERROR: First line after IncreaseDomain is not the 'delz' line!" << endl;
        exit(0);
    }
    input = input_s.substr(0,15);
    double delz;
    input_ss.str(input);
    input_ss >> delz;
    input_ss.clear();
    
    *proc.ostrm << endl << "#--- delz = " << delz << " extracted";
    *proc.ostrm << endl << "#--- call mainLine->expand";
    
    // increasing odtline with values of next neighbor
    //mainLine->outputProperties("03.1_vorExpandMainline.dat");
    mainLine->expand(delz, LUC);
    //mainLine->outputProperties("03.2_nachExpandMainline.dat");
    
    *proc.ostrm << endl << "#--- reading further data";
    while(!(ifile->eof())){
        getline(*ifile, input_s);
        //*proc.ostrm << endl << "#--- getline: >>" << input_s << "<<";
        if (input_s.size() == 0){
            continue;
        }else if (string::npos != input_s.find("endIncreaseDomain")){
            break;
        }else if (string::npos != input_s.find("phase")){
            //*proc.ostrm << endl << "#--- value for phase";
            input_ss.str(input_s.substr(0,15));
            input_ss >> phase;
            input_ss.clear();
            if(odtStats.phases.size() != 1 || odtStats.phases.at(0) != phase){
                odtP.LmultiPhase = true;
                if(odtP.LconstProp == true){
                    *proc.ostrm << endl << endl << "ERROR:" << endl;
                    *proc.ostrm << "LconsProp = 1 and multiple phases" << endl;
                    *proc.ostrm << "A multi phase simulation needs LconstProp = 0!\n";
                    exit(0);
                }
            }
        }else if (string::npos != input_s.find("rho")){
            //*proc.ostrm << endl << "#--- value for rho";
            input_ss.str(input_s.substr(0,15));
            input_ss >> rho;
            input_ss.clear();
        }else if (string::npos != input_s.find("visc")){
            //*proc.ostrm << endl << "#--- value for visc";
            input_ss.str(input_s.substr(0,15));
            input_ss >> visc;
            input_ss.clear();
        }else if (string::npos != input_s.find("lambda")){
            //*proc.ostrm << endl << "#--- value for lambda";
            input_ss.str(input_s.substr(0,15));
            input_ss >> lambda;
            input_ss.clear();
        }else if (string::npos != input_s.find("uvel")){
            //*proc.ostrm << endl << "#--- value for uvel";
            input_ss.str(input_s.substr(0,15));
            input_ss >> u;
            input_ss.clear();
        }else if (string::npos != input_s.find("vvel")){
            //*proc.ostrm << endl << "#--- value for vvel";
            input_ss.str(input_s.substr(0,15));
            input_ss >> v;
            input_ss.clear();
        }else if (string::npos != input_s.find("wvel")){
            //*proc.ostrm << endl << "#--- value for wvel";
            input_ss.str(input_s.substr(0,15));
            input_ss >> w;
            input_ss.clear();
        }else{
            *proc.ostrm << endl << "# ERROR: Line in ../input/odtChange_"
            << "*.inp could not be interpreted "
            << "for increasing domain.";
            *proc.ostrm << endl << "Line: >" << input_s << "<";
            exit(0);
        }
        
    }
    
    // setting values to mainLine
    if (LUC == 'L' || LUC == 'C'){
        *proc.ostrm << endl << "#--- change data of lower expanded area";
        mainLine->uvel.at(0)   = u;
        mainLine->vvel.at(0)   = v;
        mainLine->wvel.at(0)   = w;
        mainLine->rho.at(0)    = rho;
        mainLine->molec.at(0)  = visc;
        mainLine->lambda.at(0) = lambda;
        mainLine->phase.at(0)  = phase;
    }
    //mainLine->outputProperties("03.3_nachSetVarsLowerBC.dat");
    
    if (LUC == 'U' || LUC == 'C'){
        *proc.ostrm << endl << "#--- change data of upper expanded area";
        mainLine->uvel.at(mainLine->ngrd-1)   = u;
        mainLine->vvel.at(mainLine->ngrd-1)   = v;
        mainLine->wvel.at(mainLine->ngrd-1)   = w;
        mainLine->rho.at(mainLine->ngrd-1)    = rho;
        mainLine->molec.at(mainLine->ngrd-1)  = visc;
        mainLine->lambda.at(mainLine->ngrd-1) = lambda;
        mainLine->phase.at(mainLine->ngrd-1)  = phase;
    }
    //mainLine->outputProperties("03.4_nachSetVarsUpperBC.dat");
}


///////////////////////////////////////////////////////////////////////////////
/** This function reads a part of the change file. It changes some parameters
 *  if given.
 */
void odtSolver::ChangeParameter(ifstream *ifile){
    
    stringstream input_ss;
    string       input_s;
    bool         changePP = false;
    double       factor = 1.0;
    input_ss.clear();
    
    *proc.ostrm << endl << "#--- reading further data to change parameter";
    while(!(ifile->eof())){
        getline(*ifile, input_s);
        if (input_s.size() == 0){
            continue;
        }else if (string::npos != input_s.find("endChangeParameter")){
            *proc.ostrm << endl << "#--- finished changing parameter";
            break;
            //return;
        }else if (string::npos != input_s.find("dPdx")){
            input_ss.str(input_s.substr(0,15));
            cout << endl << "#---------- dPdx before = " << odtP.dPdx;
            //input_ss >> input_s;
            cout << endl << "#---------- the read line >" << input_s;
            input_ss >> odtP.dPdx;
            cout << endl << "#---------- dPdx after  = " << odtP.dPdx;
            input_ss.clear();
        }else if (string::npos != input_s.find("Pmax")){
            input_ss.str(input_s.substr(0,15));
            input_ss >> odtP.Pmax;
            factor = 0.1;
            input_ss.clear();
        }else if (string::npos != input_s.find("Pav")){
            input_ss.str(input_s.substr(0,15));
            input_ss >> odtP.Pav;
            factor = 0.1;
            input_ss.clear();
        }else if (string::npos != input_s.find("Z_param")){
            input_ss.str(input_s.substr(0,15));
            input_ss >> odtP.Z_param;
            factor = 0.1;
            input_ss.clear();
        }else if (string::npos != input_s.find("C_param")){
            input_ss.str(input_s.substr(0,15));
            input_ss >> odtP.C_param;
            factor = 0.1;
            input_ss.clear();
        }else if (string::npos != input_s.find("A_param")){
            input_ss.str(input_s.substr(0,15));
            input_ss >> odtP.A_param;
            factor = 0.1;
            input_ss.clear();
        }else if (string::npos != input_s.find("dtFactor")){
            input_ss.str(input_s.substr(0,15));
            input_ss >> factor;
            input_ss.clear();
        }else if (string::npos != input_s.find("Lp")){
            input_ss.str(input_s.substr(0,15));
            input_ss >> odtP.Lp;
            input_ss.clear();
            changePP = true;
        }else if (string::npos != input_s.find("Lmax")){
            input_ss.str(input_s.substr(0,15));
            input_ss >> odtP.Lmax;
            input_ss.clear();
            changePP = true;
        }else if (string::npos != input_s.find("Lmin")){
            input_ss.str(input_s.substr(0,15));
            input_ss >> odtP.Lmin;
            input_ss.clear();
            changePP = true;
        }else if (string::npos != input_s.find("bcType")){
            input_ss.str(input_s.substr(0,15));
            input_ss >> odtP.bcType;
            mainLine->setBCprops();
            input_ss.clear();
        }else{
            *proc.ostrm << endl << "# ERROR: Line in ../input/odtChange_"
            << "*.inp could not be interpreted "
            << "for changing parameter";
            *proc.ostrm << endl << "Line: >" << input_s << "<";
            exit(0);
        }
    }
    if (changePP){
        odtP.computeEddySizeDistParams();
    }
    cout << "\n#---------- Change dtSmean from " << dtSmean << " to ";
    dtSmean *= factor;
    cout << dtSmean;
    return;
}


///////////////////////////////////////////////////////////////////////////////
/** This function sets the bulk u velocity to the given value by multiplying
 *  with a factor
 *  @param input *ifile the given odtChange_*.inp file.
 */
void odtSolver::SetUMean(ifstream *ifile){
    stringstream input_ss;
    string       input_s;
    double       integ = 0.0;
    input_ss.clear();
    
    for (int i=0; i<mainLine->ngrd; i++){
        integ += mainLine->uvel.at(i)*(mainLine->posf.at(i+1)-mainLine->posf.at(i));
    }
    integ /= mainLine->Ldomain;
    //cout << "\nSetUMean_v: umean = " << integ;
    double fac = 0.0;
    
    while(!(ifile->eof())){
        getline(*ifile, input_s);
        if (input_s.size() == 0){
            continue;
        }else if (string::npos != input_s.find("endSetUMean")){
            *proc.ostrm << endl << "#--- finished setting Umean";
            for (int i=0; i<mainLine->ngrd; i++){
                mainLine->uvel.at(i) *= fac;
            }
            //integ = 0.0;
            //for (int i=0; i<mainLine->ngrd; i++){
            //    integ += mainLine->uvel.at(i)*(mainLine->posf.at(i+1)-mainLine->posf.at(i));
            //}
            //integ /= mainLine->Ldomain;
            //cout << "\nSetUMean_n: umean = " << integ;
            return;
        }else if (string::npos != input_s.find("Umean")){
            double input = 0.0;
            input_ss.str(input_s.substr(0,15));
            input_ss >> input;
            input_ss.clear();
            fac = input/integ;
            //cout << "\nSetUMean: fac   = " << fac;
        }else{
            *proc.ostrm << endl << "# ERROR: Line in ../input/odtChange_"
            << "*.inp could not be interpreted "
            << "for changing parameter";
            *proc.ostrm << endl << "Line: >" << input_s << "<";
            exit(0);
        }
    }
}
#endif


#ifdef BOOST
void odtSolver::create_subdomains_from_mainline(){

    /*********************************************************************************************************************************
     *create new sub-domains from the contents of the mainline.                                                                       *
     *--------------------------------------------------------------------------------------------------------------------------------*
     *---------------------------------------------Modifies---------------------------------------------------------------------------*
     *creates new sub-domains, thus modifies the vector sub-domains, as well as the double guardregionsize, subdomainsize, solregsize *
     *--------------------------------------------------------------------------------------------------------------------------------*
     *the odtlines will have the same values as the mainline.                                                                         *
     *                                                                                                                                *
     *********************************************************************************************************************************/


    double left, right,templ;
    int iStart,iEnd;
    int periodic;
    int tempngrd=2;
    double guardleft,guardright;
    double cellsize=mainLine->posf[1]-mainLine->posf[0];
    //to shift the left/right parameters a little bit inwards of the sub-domain
    //when using linePositionToIndex (in case they fall on a cell face boundary)

    double tol=0.00001*cellsize;
    ptr_vector<subdomain>::iterator it, it2;
//     odtline templine_comb, templine_comb2;
    solregsize=odtP.domainLength/Nsubdomains;
    guardsize=0.5*solregsize;
    limit=0.5*guardsize;
    subdomainsize=solregsize+2*guardsize;
    periodic=(odtP.bcType==1)?1:0;
    maxsimplemapsize=solregsize*0.5;

    //initialize the paramters for the subdomains,
    //i.e. the objects of class odtParam: odtP,odtP,odtP

    //Initialize the first, leftmost sub-domain
    it2 = subdomains_comb.begin();



    left=mainLine->posf[0]-periodic*guardsize;


    right=left+subdomainsize - (1-periodic)*guardsize;
    guardleft=left+periodic*guardsize;
    guardright=right-guardsize;

    iStart=mainLine->linePositionToIndex(left+tol,true);

    iEnd=mainLine->linePositionToIndex(right-tol,false);
    odtline templine_comb = odtline(mainLine,iStart,iEnd,iStart>iEnd,false);
    if(periodic)   templine_comb.addL(-odtP.domainLength);
    templ=right-left;
    odtline templine_comb2 = odtline(tempngrd,templ, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn);
    templine_comb2.addL(left);
    iStart=templine_comb2.linePositionToIndex(left,true);
    iEnd=templine_comb2.linePositionToIndex(right,false);
    templine_comb2.insertEddy(templine_comb,iStart,iEnd,left,right,false);

    //    subdomains_comb[0]=subdomain(left,right,guardleft,guardright,odtP.visc_0,0.0,0,odtP, this, mainLine->gas, mainLine->tran, mainLine->strm);


    subdomains_comb[0].set_subdomain(left,right,guardleft,guardright,0.0,0,&odtP, this, mainLine->gas, mainLine->tran, mainLine->strm);

    subdomains_comb[0].set_odtline(templine_comb2);
    subdomains_comb[0].sddiffuser.posNextSub_start= right - solregsize;
    subdomains_comb[0].sddiffuser.dx_exp= 0;

    

    //shift the newline down again (the odtline constructor has pos[iStart]=left+odtP.domainLength
    //and the guardregion was taken from the right side of the domain if jp bc are used

    //treat the middle subdomains.
    int i;
    for(i=1;i<=Nsubdomains-2;i++){
        left=i*solregsize-guardsize;

        right=left+subdomainsize;
        guardleft=left+guardsize;
        guardright=right-guardsize;

        iStart=mainLine->linePositionToIndex(left+tol,true);
        iEnd=mainLine->linePositionToIndex(right-tol,false);
        templine_comb= odtline(mainLine,iStart,iEnd,iStart>iEnd,false);
        templ=right-left;
        templine_comb2 =odtline(tempngrd,templ, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn);//, mainLine->gas);
        templine_comb2.addL(left);
        iStart=templine_comb2.linePositionToIndex(left,true);
        iEnd=templine_comb2.linePositionToIndex(right,false);
        templine_comb2.insertEddy(templine_comb,iStart,iEnd,left,right,false);
        //        subdomains_comb[i]=subdomain(left,right,guardleft,guardright,odtP.visc_0,0.0,i,&odtP,this, mainLine->gas, mainLine->tran, mainLine->strm);
        subdomains_comb[i].set_subdomain(left,right,guardleft,guardright,0.0,i,&odtP,this, mainLine->gas, mainLine->tran, mainLine->strm);
        subdomains_comb[i].set_odtline(templine_comb2);
        subdomains_comb[i].sddiffuser.posNextSub_start= right - solregsize;
        subdomains_comb[i].sddiffuser.dx_exp= 0;

	
    }



    left+=solregsize;
    right=left+subdomainsize- (1-periodic)*guardsize;

    guardleft=left+guardsize;
    guardright=right-periodic*guardsize;


    iStart=mainLine->linePositionToIndex(left+tol,true);
    iEnd=mainLine->linePositionToIndex(right-tol,false);


    templine_comb= odtline(mainLine,iStart,iEnd,iStart>iEnd,false);
    templ=right-left;
    templine_comb2 = odtline(tempngrd,templ, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn);//, mainLine->gas);
    templine_comb2.addL(left);
    iStart=templine_comb2.linePositionToIndex(left,true);
    iEnd=templine_comb2.linePositionToIndex(right,false);
    templine_comb2.insertEddy(templine_comb,iStart,iEnd,left,right,false);
    //    subdomains_comb[Nsubdomains-1]=subdomain(left,right,guardleft,guardright,odtP.visc_0,0.0,Nsubdomains-1,&odtP,this, mainLine->gas, mainLine->tran, mainLine->strm);
    subdomains_comb[Nsubdomains-1].set_subdomain(left,right,guardleft,guardright,0.0,Nsubdomains-1,&odtP,this, mainLine->gas, mainLine->tran, mainLine->strm);
    subdomains_comb[Nsubdomains-1].set_odtline(templine_comb2);
    subdomains_comb[Nsubdomains-1].sddiffuser.posNextSub_start= right - solregsize;
    subdomains_comb[Nsubdomains-1].sddiffuser.dx_exp= 0;

    

}





///////////////////////////////////////////////////////////////////////////////
/**
 * inserts a triplet map from lStart to lEnd into the appropriate sub-domains. 
 *                                                                                                        
 * -----------------------------------------Input---------------------------------------------------------
 * the start coordinate of the map, the end coordinate of the domain,                                     
 * the time of the map                                                                                    
 * ----------------------------------------Output---------------------------------------------------------
 * A boolean, true if the map extends over the right edge of the domain                                   
 * ---------------------------------------Modifies--------------------------------------------------------
 * inserts maps into the vector tripletmaps                                                               
 * -------------------------------------------------------------------------------------------------------
 *                                                                                                        
 * The map can be either simple or dependent. It is dependent when the number of guard regions            
 * affected by the domain is greater or equal two.                                                        
 * This implies that the number of affected sub-domains is greater or equal to two.                       
 * If the map is dependent, create a numbered stopevent. Save a map with this stopnumber                  
 * in the queue of each sub-domain which is affected by the map.                                          
 *                                                                                                        
 */
bool odtSolver::insert_map(double lStart, double lEnd, double time){

    bool insert_ok=true;
    vector<int> domains_affected;
    vector<int>::iterator it;

    subdomain *tsub;
    bool ok,leftok,rightok;
    int i;

    tripletmap newmap(lStart,lEnd,time); 

    //get all the sub-domains affected by the map


    return_subdomains(lStart,lEnd,domains_affected);


    //check the time contraints for all sub-domains
    //if necessary, insert the necessary stop-events in sub-domains
    //to ensure the sub-domain does not violate the time-contraint for the guard region
    //at the time of the map


    for(it=domains_affected.begin();it!=domains_affected.end();it++){

        ok=check_time_contraint(*it,time,leftok,rightok);

        if(!ok){
            updt_guardregion(*it,time,leftok,rightok);

        }
    } 
    if(lEnd>(mainLine->posf[0]+odtP.domainLength)){
        //if the map extends across the edge, all subdomain must return
        domains_affected.clear();

        for (int i=0; i<subdomains_comb.size(); i++) 
            domains_affected.push_back(subdomains_comb.at(i).sdnumber); // = i

        insert_stop_event(newmap,domains_affected,false,3);

        tsub=0;
        return true;
    }


    //check if the map is dependent
    if(is_dependent_map(lStart,lEnd,domains_affected)){
        //    else{
        insert_stop_event(newmap,domains_affected,false,1);

        tsub=0;
        return false;
        //}

    }
    /*********************************************************************************************************************************
     * In the case of a simple map, check the following cases
     * 1. Do the triplet maps on the left guardregion cover more than the limit size of the guard region (e.g. 0.5*guardregionsize)?
     * 2. Same for the right guardregion.
     *
     * //The size of domains is maxmimum 2 here, and for the two guard-regions which could be affected by the map,
     * the sub-domain needed for the update is always in the list
     **********************************************************************************************************************************/

    for(it=domains_affected.begin();it!=domains_affected.end();it++){
        tsub=&subdomains_comb.at(*it);


        //is the map in the left guard region of the current sub-domain?
        //to cover the case of the left-most sub-domain, the coordinates of
        //the left guardregion are taken modulo the domainlength.
        //since no left guard regions ever extend over the domain boundary,
        //modL will never change any other coordinates than the left guard-region
        //of the left most sub-domain
        if(lStart > modL(tsub->leftedge) && lEnd < modL(tsub->guardleft)){
            insert_ok=check_guard_region(points.at(2*(tsub->sdnumber)),lStart,lEnd);
            //return here
        }
        //since lStart E [0,Ldomain], there is no need to take
        //periodicity into account
        else if(lStart > tsub->guardright && lEnd < tsub->rightedge){
            insert_ok=check_guard_region(points.at(2*(tsub->sdnumber+1)),lStart,lEnd);
            //return here
        }
    }
    if(!insert_ok){
        insert_stop_event(newmap,domains_affected,false,5);
    }
    else{

        for(it=domains_affected.begin();it!=domains_affected.end();it++){
            tsub=&subdomains_comb.at(*it);
            //if map is in the guardregion, don't insert

            if(!(is_in_guardregion(newmap.lStart, tsub,true) || is_in_guardregion(newmap.lEnd, tsub,false)|| newmap.lStart<tsub->leftedge || newmap.lEnd>tsub->rightedge)){ // || newmap.lStart<tsub->leftedge || newmap.lEnd>tsub->rightedge) ???
                tripletmaps.at(tsub->sdnumber).push(newmap);
            }
        }

    }

    tsub=0;

    return false;

}


/**
 * Inserts a new stop-event in the the list of stopevents in the controller, insert stop events for all the sub-domains
 * in "domains," updates the subdomains (reset the values when the guard-regions are updated the latest
 * -----------------------------------------------------------------------------------------------------------------------------
 * ----------------------------------------------Input--------------------------------------------------------------------------
 * the tripletmap which triggered the stop-event, all sub-domains affected by it
 * timevent is false by default, if true the stopevent has a marker which tells the controller not to gather
 * the solution to save time 
 *
 * -------------------------------------Modifies---------------------------------------------------------------------------------
 * Insert a new stopevent into the vector stopevent, inserts stopevent maps into the heaps of the affeced sub-domains
 * clears the heaps of the inserted points (vector points) of the sub-domains in domains, resets the timelimits in the
 * subdomains.
 * ------------------------------------------------------------------------------------------------------------------------------
 *
 * the stopevents are numbered and the number corresponds to the number in the vector for thestop events for easy access.
 * The number will be returned during execution time by the subdomain.
 * which gets to the corresponding time the latest.
 * After the  new event has been created, insert_stop_event inserts a stop-map into the list of tripletmaps
 * of  the subdomain, the map has the corresponding stop-number
 */
void odtSolver::insert_stop_event(tripletmap &map, vector<int> &sdindex,bool timeevent, int type) {
    // sdindex -> domains_affected
    vector<int>::iterator it;
    //    vector<int> sdindex;
    int tsd;
    subdomain *sd;
    //  for(it=sdindex.begin();it!=sdindex.end();it++){
    //      sd=&(subdomains_comb.at(*it));
    //      sdindex.push_back(sd->sdnumber); 
    //  }

    stopevent se(map.lStart,map.lEnd,map.time,sdindex);
    se.applymap=map.applymap;
    if(timeevent){
        se.donotgather=true;
    }
    stopevents.push_back(se);
    tripletmap stopmap = map;
    stopmap.subdomains=sdindex;
    stopmap.stopnumber=stopnumber;
    stopmap.type=type;
    stopnumber++;


    //the first subdomain will only be affected on the right side, the last on the left
    //clear the points of the triplet-maps already inserted in the guard-regions,
    //reset the time of the last update of the guard-regions, insert a new stopmap
    //for all the sub-domains


    for(it=sdindex.begin(); it!=sdindex.end();it++){
        sd=&(subdomains_comb.at(*it));
        if(it==sdindex.begin()){
            sd->t0right=map.time;
            points.at((*it)*2+1)=priority_queue<point,std::vector<point>, point::greater>();
        }
        if(it ==sdindex.end()-1){
            sd->t0left=map.time;
            points.at((*it)*2)=priority_queue<point,std::vector<point>, point::greater>();

        }
        if(it!=sdindex.begin()&& it!=sdindex.end()-1){
            sd->t0left=map.time;
            sd->t0right=map.time;
            points.at((*it)*2)=priority_queue<point,std::vector<point>, point::greater>();
            points.at((*it)*2+1)=priority_queue<point,std::vector<point>, point::greater>();
        }

        tripletmaps.at(*it).push(stopmap);


    }
    sd=0;
}

bool odtSolver::is_dependent_map(double lStart, double lEnd, vector<int> &domains_affected){

    /*********************************************************************************************************
     * check if the map between lStart,lEnd is a depdentent map                                               *
     *                                                                                                        *
     * --------------------------Input------------------------------------------------------------------------*
     * the start and end coordinate of the map (lStart and lEnd), all the domains affected by the map         *
     * --------------------------Output-----------------------------------------------------------------------*
     * a boolean, true if the map is dependent, false if not -------------------------------------------------*
     * -------------------------------------------------------------------------------------------------------*
     * Count the number of total guard regions and solution regions the map affects                           *
     * if the map is bigger than some treshhold maxsimplemapsize, return true immediatly.                     *
     * maxsimplemaxsize is assigned by create_subdomain and depends on how big the guard regions              *
     * relative to the solution regions are.                                                                  *
     * This routine counts the number of guardregions affected since this number is greater or equal          *
     * to the number of solution regions                                                                      *
     *********************************************************************************************************/


    //ptr_vector<subdomain>::iterator it;
    vector<int>::iterator it;
    subdomain *sd;

    if((lEnd-lStart)>maxsimplemapsize){

        return true;
    }

    //since lStart < ldomain, if lEnd is bigger, then the map will be overlapping at least two guard regions anyway

    if(lEnd>odtP.domainLength+mainLine->pos[0]){

        return true;
    }

    int numbofgr=0;

    for(it=domains_affected.begin();it!=domains_affected.end();it++){
        sd=&subdomains_comb.at(*it);
        if( is_in_guardregion(lStart,sd,true) || is_in_guardregion(lEnd,sd,false))
        {
            numbofgr++;
            if(numbofgr>1){

                return true;
            }
        }
    }


    sd=0;
    return false;
}


void odtSolver::return_subdomains(double lStart, double lEnd,  vector<int> &domains_affected,bool neighbors){

    /****************************************************************************************************************
     * returns all the sub-domains (guard and sub-domain regions) which are affected by the map between lStart,lEnd  *
     * --------------------------------------------------------------------------------------------------------------*
     * -------------------------------------Input--------------------------------------------------------------------*
     * the start/end coordinates of the map, a reference to a vector which stores the pointers to the sub-domains    *
     * if neighbors is set to true, domains will contain the neighborung sub-domains which are not affected in       *
     * the solution region but only in the guardregions                                                              *
     * -------------------------------------Output-------------------------------------------------------------------*
     * No output, the vector domains is filled with the subdomains                                                   *
     * --------------------------------------------------------------------------------------------------------------*
     * The coordinate of a a point on the domain can be expressed as                                                 *
     *  p=posf[0]+j/N*L+x. j is the index of the solutionregion of the jth domain. domains are numbered as           *
     * 0,1,2,3..... x is the rest from the last domain to the one where the point is.                                *
     * if the point is lStart, then j=floor(lStart-posf[0]/domainlength*Nsubdomains).                                *
     * if x>solsize-guardsize means that the next subdomain is also affected.                                        *
     * if x< guardsize, that means the left subdomain is also affected.                                              *
     *                                                                                                               *
     ****************************************************************************************************************/


    //    ptr_vector<subdomain>::iterator it;
    //periodic case: all the solution regions have the same size, thus the approach described above works
    if(odtP.bcType==1){
        int j1,j2;
        double x1,x2;
        subdomain *sd1,*sd2;

        ptr_vector<subdomain> updategroup;
        domains_affected.clear();

        j1=floor((lStart-mainLine->posf[0])/odtP.domainLength*Nsubdomains);

        j2=ceil((lEnd-mainLine->posf[0])/odtP.domainLength*Nsubdomains)-1;

        x1=fmod(lStart,solregsize);
        x2=fmod(lEnd,solregsize);
        /**************************************************************************************************************************************
         *    check if the range affects also the neighbors of the subdomains (if flag is set)                                                 *
         *    x is 0 at the subdomains guardregionleft corrdinate, i.e. it is how                                                              *
         *    far the overshoot into the solution region is. if this value x is                                                                *
         *    smaller than guardregionsize, the left sub-domain is affected to (in it's right guard-region                                     *
         *    if this value is greater than guardright-guardregionsize (see illustration), then the right sub-domain                           *
         *    is affected too.                                                                                                                 *
         *      |_ _ |_ _ _ x | _ _|                                                                                                           *
         *               |_ _ |_ _ _ _|                                                                                                        *
         *                                                                                                                                     *
         *    Only add if the solution region of the  periodic image of the sub-domain is not already in the range between lStart,lEnd         *
         ***************************************************************************************************************************************/



        sd2=&subdomains_comb.at((j2+Nsubdomains)%Nsubdomains);
        sd1=&subdomains_comb.at((j1+Nsubdomains)%Nsubdomains);

        //the complicated expression for the modulo operation is necessary since % is machine dependent, for some machines
        // -12%10 will return -2, on others 8, in this case the number has to be >=0.
        int tempind=((j1-1)+Nsubdomains)%Nsubdomains;
        if(x1< guardsize && neighbors && (!((subdomains_comb.at(tempind).guardleft<(fmod(lEnd,odtP.domainLength)+mainLine->posf[0])&&
                            subdomains_comb.at(tempind).guardright>(fmod(lEnd,odtP.domainLength)+mainLine->posf[0]))))) {
            j1--;
        }

        if(x2> (solregsize-guardsize) && neighbors){
            tempind=(((j2+1)+Nsubdomains)%Nsubdomains);
            if(!((subdomains_comb.at(tempind).guardleft<lStart && subdomains_comb.at(tempind).guardright>lStart))){
                j2++;
            }
        }

        for(int i=j1;i<=j2;i++){
            domains_affected.push_back((i+Nsubdomains)%Nsubdomains);
        }

        sd1=0;sd2=0;
        //D        for(int i =0; i< updategroup.size(); i++)
        //D            delete [] updategroup.at(i);
    }

    else{
        //case of non-periodic boundary conditions: search linearly for the first sub-domain
        //and then go until the sub-domain range is no longer overlapping with the map
        //        ptr_vector<subdomain>::iterator it2,it3,it4;
        int i,j1,j2;       


        //find the first affected subdomain
        for(i=0; !(subdomains_comb.at(i).leftedge < lStart && subdomains_comb.at(i).rightedge > lStart);i++) { 
        }
        //find the last or second last affected sub-domain
        for(j1=i; !(subdomains_comb.at(j1).leftedge < lEnd && subdomains_comb.at(j1).rightedge > lEnd); j1++){
        }
        //              if(it2==it3){
        //                  it4=it2;
        //                  cout << " sdit4 "<<it4->sdnumber;
        //                  domains_comb.push_back(&(*it4));
        //              }
        j1++;

        for(j2=i;j2!=j1;){
            domains_affected.push_back(j2); 
            i++;
            if(i!=j1) j2++; else break;
        }


        //it might be that there is another sub-domain affected but the break criteria in the second for loop
        //above already holds for the one before

        if(j2!=subdomains_comb.size()-1){            
            j2++;
            if(subdomains_comb.at(j2).leftedge < lEnd){
                domains_affected.push_back(j2); 
            }
        }

    }

}






void odtSolver::gather_sol(odtline &odtl,double lStart, double lEnd, vector<int> sdindex){
    /*******************************************************************************************************
     *     get all the solution regions from lStart to lEnd together in a new odtline.                      *
     *     -------------------------------------------------------------------------------------------------*
     *     ----------------------------------------Input----------------------------------------------------*
     *     the start/end coordinate of the desired range, a vector of all the affected sub-domains          *
     *     the domains have to be ordered, that is the first sub-domain in domains is the leftmost.         *
     *     lStart is taken to be within mainLine->posf[0]/mainLine->posf[0]+odtP.domainlength                 *
     *     and has to lie within the solution region of the first sub-domain                                *
     *     only subdomains of which solution regions are affected should be in domains.                     *
     *     ----------------------------------------Output---------------------------------------------------*
     *     odtl will contain the the parts of the odtlines of the sub-domains such that                     *
     *     odtl.posf[0]==lStart, odtl.posf[odtl.ngrd]==lEnd -> odtl is shifted relative from the subdomain  *
     *     odtl's in the periodic cases                                                                     *
     *                                                                                                      *
     *    -----------------------------------------------------------------------------------------------   *
     *     in case lEnd extends over the edge of the domain (jump periodic case),                           *
     *     the odtlines of the subdomains which have positions at the beginning of the domain are           *
     *     adapted in adding L to all the positions.                                                        *
     *                                                                                                      *
     *******************************************************************************************************/


    if((lEnd-odtP.domainLength) > 1e-10){
        cout<< endl << "(lEnd-odtP.domainLength) > 1e-10 "<< "  lEnd: " <<lEnd; cout.flush();
        lEnd=odtP.domainLength;
    }

    double left,right;
    double leftlocal,rightlocal;

    subdomain *sd;
    int i1,i2;

/*    odtline templine(0,0, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn);*/
    double length = lEnd-lStart;
    odtl = odtline(10*sdindex.size(),length, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn);
    vector<int>::iterator it;
    it=sdindex.begin();
    //in case only one part of one domain is needed when gather_sol is called,
    //extract segment and return immediatly. This case is rare but can happen
    //for wrap map when gathering the range of the domains not affected by the map
    //and the map is so big that the rest is within a domain.
    if(sdindex.size()==1){
        sd=&subdomains_comb.at(*it);
        i1=(sd->odtl).linePositionToIndex(lStart,true);
        i2=(sd->odtl).linePositionToIndex(lEnd,false);
        odtl = odtline(&sd->odtl,i1,i2, i1>i2,false);
        return;
    }

    sd=&subdomains_comb.at(*it);
    left=lStart;
    right=sd->guardright;
    odtl.addL(left);

    
    i1=(sd->odtl).linePositionToIndex(left,true);
    i2=(sd->odtl).linePositionToIndex(right,false);

    odtline templine= odtline(&sd->odtl,i1,i2, i1>i2,false);
    
     i1=0;
    i2=odtl.linePositionToIndex(right,false);
    odtl.insertEddy(templine,i1,i2,left,right);
    
     
   
    //reassign left/right such that the left and right in the for loop
    //have to values corresponding to the solution region boundaries
    //if(odtP.bcType==1 || sd->sdnumber!=0){
    left=sd->guardleft;
    right=sd->guardright;


    for(it++;it<sdindex.end()-1;it++){
        left+=solregsize;
        right+=solregsize;
	
        sd=&subdomains_comb.at(*it);
	
	

        leftlocal=sd->guardleft;
        rightlocal=sd->guardright;
	
        for(int sdn=sdindex[0]; sdn<sd->sdnumber; sdn++){
            for(int i=0;i<(sd->odtl).ngrd;i++){
                (sd->odtl).pos[i]+=subdomains_comb[sdn].sddiffuser.dx_exp;
                (sd->odtl).posf[i]+=subdomains_comb[sdn].sddiffuser.dx_exp;
            }
            (sd->odtl).posf[(sd->odtl).ngrd]+=subdomains_comb[sdn].sddiffuser.dx_exp;
        }
        i1=sd->odtl.linePositionToIndex(leftlocal,true);
        i2=sd->odtl.linePositionToIndex(rightlocal,false);
        templine= odtline(&sd->odtl,i1,i2, i1>i2,false);
        if(abs(leftlocal-left) > 1e-08){
            cout << endl << "abs(leftlocal-left) > 1e-08";cout.flush();
            templine.addL(left-leftlocal);
        }
        i1=odtl.linePositionToIndex(left,true);
        //i1--;
        i2=odtl.linePositionToIndex(right,false);
        odtl.insertEddy(templine,i1,i2,left,right);
	  
   

    }



    it=sdindex.end()-1;
    sd=&subdomains_comb.at(*it);
    left+=solregsize;
    right=lEnd;
    leftlocal=sd->guardleft;
    if(lEnd-(mainLine->posf[mainLine->ngrd])>1e-08 ){
        rightlocal=lEnd-odtP.domainLength;
    }
    else{
        rightlocal=lEnd;
    }
	
    for(int sdn=sdindex[0]; sdn<sd->sdnumber; sdn++){
        for(int i=0;i<(sd->odtl).ngrd; i++){
            (sd->odtl).pos[i]+=subdomains_comb[sdn].sddiffuser.dx_exp;
            (sd->odtl).posf[i]+=subdomains_comb[sdn].sddiffuser.dx_exp;
        }
        (sd->odtl).posf[(sd->odtl).ngrd]+=subdomains_comb[sdn].sddiffuser.dx_exp;
    }
    i1=sd->odtl.linePositionToIndex(leftlocal,true);
    i2=sd->odtl.linePositionToIndex(rightlocal,false)-1;
    templine=odtline(&sd->odtl,i1,i2, i1>i2,false);
    
     
    //      if(leftlocal!=left){
    //          templine.addL(left-leftlocal,odtP);
    //      } 
    i1=odtl.linePositionToIndex(left,true);
    i2=odtl.ngrd-1;
    odtl.insertEddy(templine,i1,i2,left,right);

      
   

    double dx_exp_sum=0.0;
    for(it=sdindex.begin();it!=sdindex.end()-1;it++){
      subdomains_comb[*it].setTimeLimit(dx_exp_sum, odtP);
       dx_exp_sum +=subdomains_comb[*it].sddiffuser.dx_exp;
        //sd=&subdomains_comb.at(*it); // stupid
        //subdomains_comb[sd->sdnumber].sddiffuser.dx_exp=0.0;
        subdomains_comb[*it].sddiffuser.dx_exp=0.0;

    }



    //find data that is nan and delete  <-- CHECK: reason for nan has to be found


    for(int i=odtl.ngrd-1;i>=0;i--){
        bool isnan = odtP.Lrxn > 0 ?  odtl.yspc[3][i]!=odtl.yspc[3][i] : odtl.rho[i]!=odtl.rho[i];
        if(isnan) {
            odtl.pos.erase( odtl.pos.begin()+i, odtl.pos.begin()+i+1 );
            odtl.posf.erase( odtl.posf.begin()+i, odtl.posf.begin()+i+1 );
            odtl.rho.erase( odtl.rho.begin()+i, odtl.rho.begin()+i+1 );
            odtl.molec.erase( odtl.molec.begin()+i, odtl.molec.begin()+i+1 );
	    odtl.lambda.erase( odtl.lambda.begin()+i, odtl.lambda.begin()+i+1 );
            odtl.phase.erase( odtl.phase.begin()+i, odtl.phase.begin()+i+1 );
            for(int k=0; k<odtl.nprops; k++) {
                (*odtl.props[k]).erase( (*odtl.props[k]).begin() + i, (*odtl.props[k]).begin() + i+1 );
            }
            odtl.ngrd--;
            odtl.ngrdf--;
        }
    }
    sd=0;
}




void odtSolver::gather_sol(odtline &odtl,vector<int> sdindex){    
    /******************************************************************************************
     *---------------------------Input---------------------------------------------------------*
     * a vector with pointers to all the sub-domains, a reference to the odtline---------------*
     *--------------------------Output-------------------------------------------------------- *
     * The solution regions of all the sub-domains in domains as well as the leftmost/rightmost*
     * guard region will be put toghether as a new odtline -                                   *
     * After execution, odtl.posf[0]== domain.begin().leftedge                                 *
     * and  odtl.posf[odtl.ngrd] == (domains.end()-1).rightedge                                *
     *                                                                                         *
     ******************************************************************************************/
    double left,right;
    double leftlocal,rightlocal;
    subdomain *sd;
    int i1,i2;

/*    odtline templine(0,0, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn);*/
    double length;
    //    if((domains.begin())->sdnumber==0 || ((domains.end()-1))->sdnumber==(Nsubdomains-1))
    
    if(sdindex[0]==0 || sdindex[sdindex.size()-1]==(Nsubdomains-1))
        length = (sdindex.size())*solregsize+1*guardsize;
    else 
        length = (sdindex.size())*solregsize+2*guardsize;
    
    odtl = odtline(2*sdindex.size(),length, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn);
    //ptr_vector<subdomain>::iterator it;
    vector<int>::iterator it;
    it=sdindex.begin();

    sd=&subdomains_comb.at(*it);
    left=sd->leftedge;
    right=sd->guardright;
    odtl.addL(left);


    //      if(sd->sdnumber > 0){
    //        cout << endl << "dx_exp " << subdomains_comb[sd->sdnumber-1].sddiffuser.dx_exp << endl;
    //        for(int i=0;i<(sd->odtl).ngrd;i++){
    // 	 (sd->odtl).pos[i]+=subdomains_comb[sd->sdnumber-1].sddiffuser.dx_exp;
    // 	 (sd->odtl).posf[i]+=subdomains_comb[sd->sdnumber-1].sddiffuser.dx_exp;
    //        }
    //        (sd->odtl).posf[(sd->odtl).ngrd]+=subdomains_comb[sd->sdnumber-1].sddiffuser.dx_exp;
    //      }

    i1=0;
    i2=sd->odtl.linePositionToIndex(right,false);
    odtline templine= odtline(&sd->odtl,i1,i2, i1>i2,false);

    i1=odtl.linePositionToIndex(left,true);
    i2=odtl.linePositionToIndex(right,false);
    odtl.insertEddy(templine,i1,i2,left,right);


    if(odtP.bcType==1 || sd->sdnumber !=0){
        left=sd->guardleft;
        right=sd->guardright;
    }
    else{
        left=sd->leftedge;//+guardsize;
        right=sd->guardright;
    }

    for(it++;it<sdindex.end()-1;it++){
        left+=solregsize;
        right+=solregsize;
        sd=&subdomains_comb.at(*it);
        leftlocal=sd->guardleft;
        rightlocal=sd->guardright;
        for(int sdn=*sdindex.begin(); sdn<sd->sdnumber; sdn++){
            for(int i=0;i<(sd->odtl).ngrd;i++){
                (sd->odtl).pos[i]+=subdomains_comb[sdn].sddiffuser.dx_exp;
                (sd->odtl).posf[i]+=subdomains_comb[sdn].sddiffuser.dx_exp;
            }
            (sd->odtl).posf[(sd->odtl).ngrd]+=subdomains_comb[sdn].sddiffuser.dx_exp;
        }
        i1=sd->odtl.linePositionToIndex(leftlocal,true);
        i2=sd->odtl.linePositionToIndex(rightlocal,false);
        templine = odtline(&sd->odtl,i1,i2, i1>i2,false);
        //          if(leftlocal!=left){
        //              templine.addL(left-leftlocal,odtP);
        //          }
        i1=odtl.linePositionToIndex(left,true);
        i2=odtl.linePositionToIndex(right,false);
        odtl.insertEddy(templine,i1,i2,left,right);

    }
    it=sdindex.end()-1;
    sd=&subdomains_comb.at(*it);
    left+=solregsize;
    if(sdindex[sdindex.size()-1]==(Nsubdomains-1))
      right+=solregsize;//+guardsize;?? why +guardsize
    else
      right+=solregsize+guardsize;
	
    leftlocal=sd->guardleft;
    rightlocal=sd->rightedge;
    for(int sdn=*sdindex.begin(); sdn<sd->sdnumber; sdn++){
        cout << endl << "dx_exp " << subdomains_comb[sdn].sddiffuser.dx_exp << endl;
        for(int i=0;i<(sd->odtl).ngrd; i++){
            (sd->odtl).pos[i]+=subdomains_comb[sdn].sddiffuser.dx_exp;
            (sd->odtl).posf[i]+=subdomains_comb[sdn].sddiffuser.dx_exp;
        }
        (sd->odtl).posf[(sd->odtl).ngrd]+=subdomains_comb[sdn].sddiffuser.dx_exp;
    }
    i1=sd->odtl.linePositionToIndex(leftlocal,true);
    i2=sd->odtl.linePositionToIndex(rightlocal,false);//sd->odtl.ngrd-1;

    templine=odtline(&sd->odtl,i1,i2, i1>i2,false);



    //      if(leftlocal!=left){
    //          templine.addL(left-leftlocal,odtP);
    // 	 cout << endl <<" " <<left-leftlocal;
    //      }
        vector<double> yi(templine.nspc);


    i1=odtl.linePositionToIndex(left,true);
    i2=odtl.ngrd-1;
      
    odtl.insertEddy(templine,i1,i2,left,right);

   
    
    
    
    if((odtl.posf[i2+1]-right)>1.0e-9) { //TO DO: why odtl.posf[i2+1]-right?
      odtl.rho[odtl.ngrd-1]=odtl.rho[odtl.ngrd-2];
      odtl.molec[odtl.ngrd-1]=odtl.molec[odtl.ngrd-2];
      odtl.lambda[odtl.ngrd-1]=odtl.lambda[odtl.ngrd-2];
      odtl.phase[odtl.ngrd-1]=odtl.phase[odtl.ngrd-2];
      for(int k=0; k<odtl.nprops; k++) {
            (*odtl.props[k]).at(odtl.ngrd-1)=(*odtl.props[k]).at(odtl.ngrd-2);
            
      }
    }
    
      
    for(it=sdindex.begin();it!=sdindex.end()-1;it++){
        //sd=&subdomains_comb.at(*it); 
        //subdomains_comb[sd->sdnumber].sddiffuser.dx_exp=0.0;
        subdomains_comb[*it].sddiffuser.dx_exp=0.0;
    }

    for(int i=odtl.ngrd-1;i>=0;i--){

        bool isnan = odtP.Lrxn > 0 ?  odtl.yspc[3][i]!=odtl.yspc[3][i] : odtl.rho[i]!=odtl.rho[i];
        if(isnan) {
            odtl.pos.erase( odtl.pos.begin()+i, odtl.pos.begin()+i+1 );
            odtl.posf.erase( odtl.posf.begin()+i, odtl.posf.begin()+i+1 );
            odtl.rho.erase( odtl.rho.begin()+i, odtl.rho.begin()+i+1 );
            odtl.molec.erase( odtl.molec.begin()+i, odtl.molec.begin()+i+1 );
	    odtl.lambda.erase( odtl.lambda.begin()+i, odtl.lambda.begin()+i+1 );
            odtl.phase.erase( odtl.phase.begin()+i, odtl.phase.begin()+i+1 );
            for(int k=0; k<odtl.nprops; k++) {
                (*odtl.props[k]).erase( (*odtl.props[k]).begin() + i, (*odtl.props[k]).begin() + i+1 );
            }
            odtl.ngrd--;
            odtl.ngrdf--;
        }
    } 

    sd=0;
}

void odtSolver::distribute_sol(sdgroup &group, bool wrap){

    /*******************************************************************************************************************
     *distribute the odtline in group and the DA vector in group to the sub-domain.                                     *
     *------------------------------------------------------------------------------------------------------------------*
     *-------------------------------------------------Input------------------------------------------------------------*
     *A reference to the group. I case warp=true, the groups odtl should be in the same range as mainLine,              *
     *i.e. odtl.posf[0]=mainLine->posf[0] && odtl.posf[odtl.ngrd]==odtl.posf[odtl.ngrd], and odtl should only consist    *
     *of solution regions                                                                                               *
     *------------------------------------------------Output------------------------------------------------------------*
     *No output                                                                                                         *
     *-----------------------------------------------Modifies-----------------------------------------------------------*
     *creates the new sub-domains and copies them into the vector subdomains                                            *
     *------------------------------------------------------------------------------------------------------------------*
     *This function has to be modified for non-periodic boundary conditions because the indexing does not work          *
     *******************************************************************************************************************/

    odtline templine_comb(1,1, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn),templine_comb2(1,1, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn)/*,odtl(0,0, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn)*/;
    odtline odtl=group.odtl;

    int i1,i2;
    double left,right,guardleft,guardright,left2,right2,temptsav,oleft,oright;
    int sdnumber;
    vector<int>::iterator sdit;

    int counter=0;


    sdit=group.subdomains.begin();
    sdnumber=*sdit;
    if(sdnumber==0 && odtP.bcType!=1){
        templine_comb2 = odtline(5,subdomainsize-guardsize, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn);
    }
    else{
        templine_comb2 = odtline(10,subdomainsize, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn);
    }

    if(wrap){

        templine_comb2.addL(mainLine->posf[0]-guardsize);
        //wrap case (periodic boundary conditions)
        //first get the solution of the right side of the domain:
        left=(mainLine->posf[mainLine->ngrd]-guardsize);
        i1=odtl.linePositionToIndex(left,true);
        i2=odtl.ngrd-1;
        templine_comb = odtline(&odtl,i1,i2, i1>i2,false);

        //insert the line from the end of the domain in the beginning of templine2
        //corresponding to the beginning of the domain
        i1=0;
        i2=templine_comb2.linePositionToIndex(mainLine->posf[0],false);
        left=mainLine->posf[0]-guardsize;
        right=mainLine->posf[0];
        //shift the piece down
        templine_comb.addL(-odtP.domainLength);
        templine_comb2.insertEddy(templine_comb,i1,i2,left,right,i1>i2);

        //get the rest from the beginning of the domain, insert it after the piece in templine2:
        left=mainLine->posf[0];
        right=mainLine->posf[0]+solregsize+guardsize;

        i1=odtl.linePositionToIndex(left,true);
        i2=odtl.linePositionToIndex(right,false);
        templine_comb=odtline(&odtl,i1,i2, i1>i2,false);
        i1=templine_comb2.linePositionToIndex(left,true);
        i2=templine_comb2.linePositionToIndex(right,false);
        templine_comb2.insertEddy(templine_comb,i1,i2,left,right,i1>i2);

        left=mainLine->posf[0]-guardsize;
        right=left+subdomainsize;
        oleft=left;
        oright=right;
    }
    else{
        //calculate the sd-number. the point lStart+1.1*guardsize is a)  >= 0 and falls into the solution region of the sub-domain
        //from which the number is to be calculated

        left=odtl.posf[0];

        if(sdnumber==0 && odtP.bcType!=1)
            right=left+subdomainsize-guardsize;
        else
            right=left+subdomainsize;
        left2=subdomains_comb.at(sdnumber).leftedge;
        right2=subdomains_comb.at(sdnumber).rightedge;

        templine_comb2.addL(left2);
        i1 =0;
        i2 = odtl.linePositionToIndex(right,false);
        templine_comb = odtline(&odtl,i1,i2, i1>i2,false);
        if(left!=left2){
            templine_comb.addL(left-left2);
        }
 
 
 
        i1=templine_comb2.linePositionToIndex(left2,true);
        i2=templine_comb2.linePositionToIndex(right2,false);
        templine_comb2.insertEddy(templine_comb,i1,i2,left2,right2,i1>i2);

        oleft=left;
        oright=right;
    }
    counter++;

    left2=subdomains_comb.at(sdnumber).leftedge;
    right2=subdomains_comb.at(sdnumber).rightedge;

    guardleft=subdomains_comb.at(sdnumber).guardleft;
    guardright=subdomains_comb.at(sdnumber).guardright;

        
    subdomains_comb[sdnumber].set_subdomain(left2,right2,guardleft,guardright,group.curtime,sdnumber,&odtP,this, mainLine->gas, mainLine->tran, mainLine->strm);
    subdomains_comb[sdnumber].set_odtline(templine_comb2);

 
  

    if(sdnumber==0 && odtP.bcType!=1)
        left=oleft+counter*solregsize-guardsize;
    else
        left=oleft+counter*solregsize;
    right=oright+counter*solregsize;


    for(sdit=(group.subdomains.begin()+1);sdit<(group.subdomains.end()-1);sdit++){
        sdnumber = *sdit;

        left2=subdomains_comb.at(sdnumber).leftedge;
        right2=subdomains_comb.at(sdnumber).rightedge;


        i1 = odtl.linePositionToIndex(left,true);
        i2 = odtl.linePositionToIndex(right,false);

        templine_comb = odtline(&odtl,i1,i2, i1>i2,false);


        templine_comb2 = odtline(10,subdomainsize, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn);
        templine_comb2.addL(left2);


        if(left2-left > 1.0e-9){
            templine_comb.addL(left2-left);
        }

        i1=templine_comb2.linePositionToIndex(left2,true);
        i2=templine_comb2.linePositionToIndex(right2,false);


        templine_comb2.insertEddy(templine_comb,i1,i2,left2,right2,i1>i2);

        guardleft=left2+guardsize;
        guardright=right2-guardsize;

        temptsav=subdomains_comb[sdnumber].tsaverage;
        subdomains_comb[sdnumber].set_subdomain(left2,right2,guardleft,guardright,group.curtime,sdnumber,&odtP,this, mainLine->gas, mainLine->tran, mainLine->strm);
        subdomains_comb[sdnumber].tsaverage=temptsav;



        subdomains_comb[sdnumber].set_odtline(templine_comb2);
   
	
        counter++;
        left=oleft+counter*solregsize;
        right=oright+counter*solregsize;
    }

    //last domain: get part from the beginning of the domain when it was a wrapped map
    sdnumber=*sdit;
    left2=subdomains_comb.at(sdnumber).leftedge;
    right2=subdomains_comb.at(sdnumber).rightedge;
    if(sdnumber==(Nsubdomains-1) && odtP.bcType!=1){
        templine_comb2 = odtline(5,subdomainsize-guardsize, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn);
    }
    else{
        templine_comb2 = odtline(10,subdomainsize, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn);
    }
    templine_comb2.addL(left2);

    //remark: the left/right are taken from mainLine since errors in the numerics otherwise hard to detect might show up in inconsitensies between mainLine and
    //templine2/odtl, resulting in an "Error in Line Position To Index" from anyline when calculating the indexes of the lines
    //in templine2/odtl
    if (wrap){
        //get the solution from the beginning of the domain
        right=mainLine->posf[0]+guardsize;
        left=mainLine->posf[0];
        i1=0;
        i2=odtl.linePositionToIndex(right,false);
        templine_comb = odtline(&odtl,i1,i2, i1>i2,false);

        //this left and right are where the piece of the odtl is put in the domain
        left=mainLine->posf[mainLine->ngrd];
        right=mainLine->posf[mainLine->ngrd]+guardsize;
        i1=templine_comb2.linePositionToIndex(left,true);
        i2=templine_comb2.ngrd-1;
        templine_comb.addL(odtP.domainLength);
        templine_comb2.insertEddy(templine_comb,i1,i2,left,right,i1>i2);

        //get the rest from the end of the domain:
        left=mainLine->posf[mainLine->ngrd]-solregsize-guardsize;
        right=mainLine->posf[mainLine->ngrd];

        i1=odtl.linePositionToIndex(left,true);
        i2=odtl.ngrd-1;
        templine_comb=odtline(&odtl,i1,i2, i1>i2,false);
        i1=templine_comb2.linePositionToIndex(left,true);
        i2=templine_comb2.linePositionToIndex(right,false);
       	
	templine_comb2.insertEddy(templine_comb,i1,i2,left,right,i1>i2);

        left=mainLine->posf[mainLine->ngrd]-solregsize-guardsize;
        right=mainLine->posf[mainLine->ngrd]+guardsize;
        guardleft=left+guardsize;
        guardright=right-guardsize;

    }
    else{
        //left and right already assigned

	
        i1=odtl.linePositionToIndex(left,true);
	if(sdnumber== Nsubdomains-1 && odtP.bcType!=1)
	  right-=guardsize;
        i2 = odtl.linePositionToIndex(right,false);//i2=odtl.ngrd-1;

        templine_comb = odtline(&odtl,i1,i2, i1>i2,false);

        if(left2-left > 1.0e-9){
            templine_comb.addL(left2-left);
        }

        i1=templine_comb2.linePositionToIndex(left2,true);
        i2=templine_comb2.linePositionToIndex(right2,false);
	

	templine_comb2.insertEddy(templine_comb,i1,i2,left2,right2,i1>i2);

      
        if(sdnumber==(Nsubdomains-1) && odtP.bcType!=1){
            guardleft=left2+guardsize;
            guardright=right2;
        }
        else{
            guardleft=left2+guardsize;
            guardright=right2-guardsize;
        }
    }
    temptsav=subdomains_comb[sdnumber].tsaverage;
    subdomains_comb[sdnumber].set_subdomain(left2,right2,guardleft,guardright,group.curtime,sdnumber,&odtP,this, mainLine->gas, mainLine->tran, mainLine->strm);
    subdomains_comb[sdnumber].tsaverage=temptsav;
    

    subdomains_comb[sdnumber].set_odtline(templine_comb2);
   



}

bool odtSolver::check_guard_region(priority_queue<point,vector<point>, point::greater> &points, double lStart, double lEnd){

    /*********************************************************************************************************************************
     *checks if in the heap "points", after inserting lStart and lEnd, there is a segment > than limit (=0.5 guardregion size e.g.)   *
     *(a segment is between a start and an endpoint)                                                                                  *
     *--------------------------------------------------------------------------------------------------------------------------------*
     *---------------------------------------------------Input------------------------------------------------------------------------*
     *the coordinate lStart of the first point of the new segment, the coordinate lEnd of the last point of the segment               *
     *a reference to the heap where the previously inserted points are stored                                                         *
     *--------------------------------------------------Output------------------------------------------------------------------------*
     *a boolean which is true if the segment was inserted such that no segment is > the limit. If the result is true, also a new heap *
     *is generated which contains the new points and does not contain all the points which are covered by the new segment. The new    *
     *heap is copied onto points                                                                                                      *
     *--------------------------------------------------------------------------------------------------------------------------------*
     *A point contains a coordinate and a boolean flag if the point is at the beginning of the segment                                *
     *                                                                                                                                *
     *********************************************************************************************************************************/


    int active_segments=0;
    priority_queue<point,vector<point>, point::greater > newpoints;
    point p1(lStart,true);
    point p2(lEnd,false);
    point toppoint,firstpoint;

    points.push(p1);
    points.push(p2);

    while (!points.empty()){
        toppoint=points.top();
        points.pop();
        if(active_segments==0){
            firstpoint=toppoint;
        }
        active_segments=(toppoint.isbeginning)?active_segments+1:active_segments-1;

        if(active_segments==0){
            if((toppoint.pos-firstpoint.pos)>limit){
                return false;
            }
            else{
                newpoints.push(firstpoint);
                newpoints.push(toppoint);
            }
        }
    }

    points=newpoints;
    return true;
}

/** 
 * when called from a sub-domain, this funcion gets the event with index stopnumber from the vector, assembles
 * the sub-domain group, and then (if needed) applies the triplet map, and redistributes the updated odtline.
 * -------------------------------------------------------------------------------------------------------------------
 * -------------------------------------------Input-------------------------------------------------------------------
 * The index of the event
 * -------------------------------------------Modifies----------------------------------------------------------------
 * May overwrite the mainLine when a map over the edge is execute
 * ---------------------------------------------------------------------------------------------------------------
 * there are two main cases
 * dependent maps and update of the guard regions are essentially the same (the map is not executed when the
 * update of the guard region is due to the violation of the time-contstraint
 * then there is the case when a map exends over the edge of the domain.
 * the needed range looks like this (the numbered cells are the ones affected by the map.
 * they are gathered from the corresponding subdomains)
 * 
 *  |4|5|6|x|x|x|x|x|x|1|2|3|
 * 
 *   make a new line from these cells, the one at the beginning of the domain are shifted with adding the domain
 *   length to their pos and posf vectors and adding length/Ldomain*jump to their scalar vectors (e.g. uvel), jump
 * is the change of value in the scalar the domainlength to the right.
 * to this new line, the tripletmap is applied. then the mainline is reinitialized and from there, the new subdomains
 * are created.
 */
void odtSolver::apply_event(int stopnumber){

  
      
    t1=wt.walltime_return(&clockzero);

    //ptr_vector<subdomain>::iterator it;
    vector<int>::iterator it;

    int i1,i2;
    double l1,l2;
    bool wrapmap=false;
    //      vector<subdomain> domains, domains2;
    //ptr_vector<subdomain> domains_comb, domains2_comb;
    vector<int> sdindex;
    vector<int> sdindex2;

    odtParam sdgparam(odtP);


    //      odtline eddyline,templine, templine2;
    odtline eddyline_comb(1,1, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn),templine_comb(1,1, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn), templine_comb2(1,1, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn);

    stopevent se=stopevents.at(stopnumber);

    sdgroup tempgroup(0,0,0,0,0,sdgparam,this,se.subdomains, mainLine->gas, mainLine->tran, mainLine->strm);


    //    domains_comb = get_locations(se.subdomains);  // method is deprecated !!
    //    WORK TO DO HERE !!!!
    //    // stopevents -> ptr_vector !?!?!??!
    //    se.subdomains: vector<int> of sdindex
    int i=0;
    sdindex = se.subdomains;
    for (it=sdindex.begin();it<sdindex.end();it++)
        subdomains_comb.at(*it).stopnumber=-1;
	

    /*  //ptr_vector<sd> domains_comb replaced by vector<int> !!
        for(it=domains_comb.begin();it<domains_comb.end();it++){
        sdnumbers[i]=(&(*it))->sdnumber;
        (&(*it))->stopnumber=-1;
        i++;
        }
        */
    //          for(int i=1;i<domains_comb.size();i++){
    // 	   if (sdnumbers[i]!=sdnumbers[i-1]+1 && sdnumbers[i]!=sdnumbers[i-1]-1) return;} //???

    if(se.donotgather==true){
        return;
    }

    if(se.lEnd>(mainLine->posf[0]+odtP.domainLength)){
        wrapmap=true;
	cout << endl << " wrapmap true";cout.flush();
    }
    if(wrapmap){

        //if the map extends over the edges (cells with * is the map)
        // |*|*|*| | | |*|*|
        // first get all * cells from the subdomains, a new odtline with posf[0]=lStart,
        //posf[ngrd]=lEnd
        //then get the rest of the cells. odtline with range posf[0]=lEnd-odtP.domainLength, lStart
        //insert both this lines into mainLine

        //get eddyline and apply triplet map
        return_subdomains(se.lStart,se.lEnd,sdindex2,false);


        gather_sol(eddyline_comb,se.lStart,se.lEnd,sdindex2);

        //get the rest of the cells (also of sub-domains not affected by the event) insert into mainLine
        l1=se.lEnd-odtP.domainLength;
        l2=se.lStart;
        return_subdomains(l1,l2,sdindex,false);
        //        subdomain *sd2_comb;
        int sd2_comb;        



        sd2_comb=*sdindex.begin();
        if(l1<subdomains_comb.at(sd2_comb).leftedge){
            sdindex.erase(sdindex.begin());
            sdindex.push_back(sd2_comb);
        }
        gather_sol(templine_comb2,l1,l2,sdindex);

        templine_comb =odtline(10,odtP.domainLength, mainLine->gas, mainLine->tran, mainLine->strm, &odtP, "", LhasVel, odtP.Lrxn);
        templine_comb.addL(se.lEnd-odtP.domainLength);

        l1=se.lEnd-odtP.domainLength;
        l2=se.lStart;
        i1=0;


        i2=templine_comb.linePositionToIndex(se.lStart,false);
        templine_comb.insertEddy(templine_comb2,i1,i2,l1,l2,i1>i2);



        //insert eddyline
        l1=se.lStart;
        l2=se.lEnd;
        i1=templine_comb.linePositionToIndex(l1,true);
        i2=templine_comb.ngrd-1;

        templine_comb.insertEddy(eddyline_comb,i1,i2,l1,l2,i1>i2);

        tempgroup = sdgroup(templine_comb.posf[0],templine_comb.posf[templine_comb.ngrd],se.time,se.lStart,se.lEnd,sdgparam,this,se.subdomains, mainLine->gas, mainLine->tran, mainLine->strm);


        tempgroup.set_odtline(templine_comb);

        sdindex.clear();
        for(int i=0;i<Nsubdomains;i++)
            sdindex.push_back(i);


        tempgroup.apply_maps();

        tempgroup.shift_group(mainLine->posf[0]-tempgroup.odtl.posf[0]);

        distribute_sol(tempgroup,true);
        mainLine=&tempgroup.odtl;

        //        sd2_comb=0;
    }
    else{

      
        gather_sol(templine_comb,sdindex);
      

        stopnumber=-1;

        tempgroup = sdgroup(templine_comb.posf[0],templine_comb.posf[templine_comb.ngrd],se.time,se.lStart,se.lEnd,sdgparam,this,se.subdomains, mainLine->gas, mainLine->tran, mainLine->strm);
        tempgroup.set_odtline(templine_comb);
        if(se.applymap){
            tempgroup.apply_maps();
        }
    
	stopnumber=-1;
        distribute_sol(tempgroup);
      
	
    }

    if(!wrapmap){
        //domains_comb = get_locations(se.subdomains);
        sdindex = se.subdomains;
        for(it=sdindex.begin();it!=sdindex.end();it++){
            distribute_maps(subdomains_comb.at(*it));
        }
    }
    if(wrapmap){
        timeincontrol+=wt.walltime_return(&t1);
    }
    else{
        stopeventtime+=wt.walltime_return(&t1);
    }



}


/**********************************************************************************************************
 * copies the maps in the vector tripletmaps into the list in the subdomains                               *
 * --------------------------------------------------------------------------------------------------------*
 * ---------------------------------------------Input------------------------------------------------------*
 * a reference to a sub-domains                                                                            *
 * --------------------------------------------Modifies----------------------------------------------------*
 * deletes maps from the heap corresponding to the subdomain in the vector subdomains                      *
 * --------------------------------------------------------------------------------------------------------*
 **********************************************************************************************************/

void odtSolver::distribute_maps(subdomain &sd){

    int sdn=sd.sdnumber;
    tripletmap map;
    while(!(tripletmaps.at(sdn).empty())){
        distributed_maps++;
        map=tripletmaps.at(sdn).top();
        tripletmaps.at(sdn).pop();
        sd.maps.push_back(map);
        if(map.stopnumber>=0){
            break;
        }
    }

}



//ptr_vector<subdomain> odtSolver::get_locations(vector<int> &sdindexes){
/**********************************************************************************************
 * --------------------------------------------------------------------------------------------*
 * ----------------------------------Input-----------------------------------------------------*
 * a vector with the indexes of some sub-domains                                               *
 * ----------------------------------Output----------------------------------------------------*
 *  vector with the pointers to the sub-domains int the same order as in sdindexes             *
 * --------------------------------------------------------------------------------------------*
 **********************************************************************************************/


//ptr_vector<subdomain> domains;
//vector<int>::iterator it;

//for(it=sdindexes.begin();it!=sdindexes.end();it++){
//domains.push_back(&subdomains_comb[*it]);
//}

//return domains;
//[>
//for(int i =0; i< domains.size(); i++){
//domains.at(i)=0;}
//*/
//}


double odtSolver::modL(const double &inL){
    /********************************************************************************************
     *-------------------------------------------------------------------------------------------*
     *-------------------------------------Input-------------------------------------------------*
     *some number in double format                                                               *
     *------------------------------------Output-------------------------------------------------*
     *the the input number modulo the domainlength but within mainLine->posf[0]                   *
     *and mainLine->posf[mainLine->ngrd]                                                           *
     *-------------------------------------------------------------------------------------------*
     ********************************************************************************************/

    double templ=inL;

    if(inL<mainLine->posf[0]){
        templ+=odtP.domainLength;
    }
    else if(inL>mainLine->posf[mainLine->ngrd]){
        templ-=odtP.domainLength;
    }

    return templ;

}

bool odtSolver::is_in_guardregion(const double &inp, const subdomain *sd, bool left){
    /************************************************************************************************************************************
     *checks if the point inp is in the left guardregion of sd (if left=true), else if inp is in the right guard region of sd            *
     *---------------------------------------------------------------------------------------------------------------------------------- *
     *---------------------------------------------------------Input-------------------------------------------------------------------- *
     * the coordinate of the point, a pointer to the sub-domain, a boolean indicating if the right or left guardregion should be checked *
     *---------------------------------------------------------Output------------------------------------------------------------------- *
     *a boolean (true if inp is in the guardregion                                                                                       *
     *---------------------------------------------------------------------------------------------------------------------------------  *
     *left guard region: in case the point is between [domainlength-guardsize, domainlength], it falls into the guardregion of the left  *
     *subdomain too                                                                                                                      *
     * since modL returns L for input L and 0 for input 0, special treatment is needed.                                                  *
     ************************************************************************************************************************************/


    if(left){

        if((inp>=modL(sd->leftedge) && inp<=modL(sd->guardleft)) || ((inp>=modL(sd->leftedge)&& (sd->guardleft==0) && inp<=(sd->guardleft+odtP.domainLength)) && odtP.Lperiodic)) {
            return true;
        }

    }
    //add domain length to inp rather than modL(sd->righedge) since modL(sd->guardright) would need a additional comparision as above (since
    //modL returns domainLength with input domainLength)
    //but numerical accuracy when comparing sd->guardright==mainLine->posf[mainLine->ngrd] cannot be guaranteed.
    else if((inp>=sd->guardright && inp< sd->rightedge) || ((inp+odtP.domainLength)>=sd->guardright && (inp+odtP.domainLength) < sd->rightedge  && odtP.Lperiodic)){
        return true;
    }

    return false;

}


void odtSolver::vector_out(vector<double> &v){
    //prints out a vector
    for(int i=0;i<v.size();i++){
        cout << i << " " << v[i] << endl;
    }
    cout << endl;

}


bool odtSolver::check_time_contraint(int sdindex, double time, bool &leftok, bool &rightok){
    subdomain *sd;
    sd = &subdomains_comb.at(sdindex);
    //check if the guardregions of sd need an update due to the time contraint. leftok is true if
    //the left guardregion does not need an update at time "time", the same holds for rightok
    // the return value is just a logical and of leftok and rightok
    leftok=true;
    rightok=true;

    if((odtP.Lperiodic || sd->sdnumber!= 0) && (sd->timelimit-(time-sd->t0left))<1e-10){
        leftok=false;
    }
    else if(!odtP.Lperiodic && sd->sdnumber== 0) leftok=true;


    if((odtP.Lperiodic || sd->sdnumber!= Nsubdomains-1) &&time-sd->t0right>sd->timelimit){
        rightok=false;
    }
    else if(!odtP.Lperiodic && sd->sdnumber== Nsubdomains-1) rightok=true;

    sd = 0;
    return leftok*rightok;

}



/**
 * makes sure the subdomain "sd" fullfills the time-contraints for the guardregions
 * ----------------------------------------------------------------------------------------------------------------
 * ---------------------------------------------Input--------------------------------------------------------------
 * a reference to the subdomain, a double for the time when the constraints needs to be fullfilled,
 * two booleans guardleftok, guardrightok which tell the function which sides are NOT ok at the time "time"
 * -------------------------------------------Modifies-------------------------------------------------------------
 * calls insert_stopevent
 * ----------------------------------------------------------------------------------------------------------------
 * The following assertion hold for this functions:
 * at the beginning, the guardregion indicated by guardleftok/guardleftok violate the time contstraint.
 * at function exit, this sub-domain will not violate the time contraint.
 * 
 * when this function is called, the subdomain at time "time" violates the time-constraint for the sides indicated
 * by the booleans guardleftok/guardrightok (false if the constraint is violated).
 * At the end of the function the function
 * will insert a new stopevent at the time (t0left|t0right)+limit (from now on called "updatetime"),
 * which will update the corresponding guardregion  once the sub-domains begin to apply
 * the maps. Before inserting this stopevent, this function must be sure that
 * a) the time between the inserted stopevent and the time "time" is smaller than  the timelimit, and
 * b) that the neighbors do not violate the timeconstraint at the time "updatetime".
 * This in turn applies of course to every sub-domain checked. Thus, the function calls itself recursively.
 *
 */
void odtSolver::updt_guardregion(int sdindex, double time,bool guardleftok, bool guardrightok) {

    bool isleft;
    bool isright;
    bool isok;
    int sdnleft,sdnright;
    double newtime;
    //ptr_vector<subdomain> updategroup;
    vector<int> updategroup;

    subdomain *sd;
    sd = &subdomains_comb.at(sdindex);

    //check if the domains sd fullfills the time contraints at time "newtime"
    if(sd->sdnumber!=0 && sd->sdnumber!=Nsubdomains-1) newtime=min(sd->t0left,sd->t0right)+sd->timelimit; //C.Schroedinger
    else if(sd->sdnumber==0) newtime=sd->t0right+sd->timelimit;
    else newtime=sd->t0left+sd->timelimit;
    sd->t0left=(guardleftok)?sd->t0left:newtime;
    sd->t0right=(guardrightok)?sd->t0right:newtime;



    //now make sure both the left and the right neighbor are in a state at time "time" where they do not violate the time constraint


    /////////////////get the neighboring domains which are needed to perform an update///////////////////////////////


    //find the left neighbor of the domain. Depending on the type of bc, the left/right most sub-domain might have
    //have neighbors or not.
    if(odtP.bcType==1){
        sdnleft=(sd->sdnumber!=0)?sd->sdnumber-1:Nsubdomains-1;
    }
    //in case non-periodic boundary conditions are applied, subtract one and catch the exeption for
    //the boundaries later
    else{
        sdnleft=sd->sdnumber-1;
    }

    //find the right neighbor of the domain.
    if(odtP.bcType==1){
        sdnright=(sd->sdnumber+1)%Nsubdomains;
    }
    else{
        sdnright=sd->sdnumber+1;
    }


    //////////check the left neighbor//////////////////////////////////////////////////////////////////
    //in case non-periodic boundary conditions are applied, don't check if the sub-domain is the left/rightmost

    if(!( odtP.bcType!=1 && sdnleft ==-1) ||  (odtP.bcType!=1 && sdnleft !=-1) ){  //C.Schroedinger
        isok=check_time_contraint(sdnleft,newtime,isleft,isright);

        //the left neighbor might need to be updated. pass which side needs an update
        if(!isok){
            updt_guardregion(sdnleft,newtime,isleft,isright);
        }

    }

    ////////check the right neighbor //////////////////////////////////////////////////////////////////////

    if(!((sdnright==Nsubdomains) && odtP.bcType!=1)||  (odtP.bcType!=1 && sdnright!=Nsubdomains)){  //C.Schroedinger

        isok=check_time_contraint(sdnright,newtime,isleft,isright);
        if(!isok){
            updt_guardregion(sdnright,newtime,isleft,isright);
        }

    }


    //needs left neighbor for update

    if(!guardleftok && sdnleft !=-1){
        updategroup.push_back(subdomains_comb.at(sdnleft).sdnumber); // = sdnleft
    }
    updategroup.push_back(subdomains_comb.at(sd->sdnumber).sdnumber); // = sd->sdnumber == sdindex 

    //needs right neighbor for update
    if(!guardrightok && sdnright!=Nsubdomains){
        updategroup.push_back(subdomains_comb.at(sdnright).sdnumber); // sdnright
    }

    tripletmap newmap(0,0,newtime);
    newmap.applymap=false;
    insert_stop_event(newmap,updategroup,false,2);

    isok=check_time_contraint(sdindex,time,isleft,isright);
    if(!isok){
        updt_guardregion(sdindex,time,isleft,isright);
    }



}
#endif


///shifting the domain for premixed cases and inflow on one side; TO DO

void odtSolver::shiftFlame(/*double x_0, double t_0, double t, std::ofstream &sL_file, int iEtrials*/) {

 /*   if(odtP.constShiftSpeed){
        double dx_shift = odtP.ShiftSpeed* (t-t_0);

        int i;
        double phi_in = odtP.phi, temp_in = odtP.tempIn;
        if(odtP.timeInput){
            for(i=0; i< lengthTimeData; i++){
                if(timedata_time[i]> t) break;
            }
        }
        if(odtP.timeInput && i==0) phi_in = timedata_phi[0];
        else if(odtP.timeInput && i==lengthTimeData-1) phi_in = timedata_phi[lengthTimeData-1];
        else if(odtP.timeInput){
            phi_in = timedata_phi[i-1] + (timedata_phi[i] -timedata_phi[i-1])/(timedata_time[i] -timedata_time[i-1]) * (t-timedata_time[i-1]);
        }
        else if(odtP.SinusTimeInput && odtP.timeInpSwitch){
            if (odtP.timeInpSwitch == 1)  // phi
                phi_in = odtP.phi + odtP.SinusAmpl*odtP.phi * sin(2*PI*odtP.sinusFreq*t);
            else if (odtP.timeInpSwitch == 2) 
                temp_in = odtP.tempIn + odtP.SinusAmpl*odtP.tempIn * sin(2*PI*odtP.sinusFreq*t);
        }
        cout << endl << phi_in << " " << temp_in << endl;
        combLine.expand1(dx_shift, 0.0, phi_in, temp_in);
        lineDiffuser.chopOutflowGrid(0.0, odtP.domainLength);
        if(neddies % (odtP.modDisp) == 0)  *proc.ostrm << endl << "dx_shift = " << dx_shift;
    }
    else if(!odtP.sL_spec){            // calculate dx_shift finding difference of position of mean integrated temperature and Ldomain/2

        double x;

        double Sum, Prod;

        int CH4 = combLine.gas->speciesIndex("CH4");

        Sum = 0.0;
        Prod = combLine.yspc[CH4][0] * combLine.rho[0];//*(odtS.combLine.posf[1]-odtS.combLine.posf[0]);
        //Prod = combLine.temp[combLine.ngrd-1]-combLine.temp[0];
        for(int i=0; i<combLine.ngrd; i++) {
            Sum += combLine.yspc[CH4][i] * combLine.rho[i]* (combLine.posf[i+1]-combLine.posf[i]);
            //Sum += (combLine.temp[i]-combLine.temp[0]) *(combLine.posf[i+1]-combLine.posf[i]);
        }
        x = Sum/Prod;// *odtS.combLine.Ldomain/odtS.combLine.ngrd;combLine.Ldomain -

        double dx_shift = x_0-x;
        if(neddies % (odtP.modDisp) == 0) *proc.ostrm << endl << "dx_shift = " <<dx_shift;
        double sL =   (x_0-x)/(t-t_0);
        if(neddies % (odtP.modDisp) == 0) *proc.ostrm << endl <<"sL: "<< sL << endl;

        if(dx_shift > 0.0){
            { 
                int i;
                double phi_in = odtP.phi, temp_in = odtP.tempIn;
                if(odtP.timeInput){
                    for(i=0; i< lengthTimeData; i++){
                        if(timedata_time[i]> t) break;
                    }
                }
                if(odtP.timeInput && i==0 ) phi_in = timedata_phi[0];
                //else if(i==lengthTimeData-1) phi_in = timedata_phi[lengthTimeData-1];
                else if(odtP.timeInput){  
                    phi_in = timedata_phi[i-1] + (timedata_phi[i] -timedata_phi[i-1])/(timedata_time[i] -timedata_time[i-1]) * (t-timedata_time[i-1]);}
                else if(odtP.SinusTimeInput && odtP.timeInpSwitch){
                    if (odtP.timeInpSwitch == 1)  // phi
                        phi_in = odtP.phi + odtP.SinusAmpl*odtP.phi * sin(2*PI*odtP.sinusFreq*t);
                    else if (odtP.timeInpSwitch == 2) 
                        temp_in = odtP.tempIn + odtP.SinusAmpl*odtP.tempIn * sin(2*PI*odtP.sinusFreq*t);
                }
                else  phi_in =odtP.phi;

                combLine.expand1(dx_shift, 0.0, phi_in, temp_in);}

                lineDiffuser.chopOutflowGrid(0.0, odtP.domainLength);}


                if(iEtrials%50==0){
                    sL_file << scientific;
                    sL_file << setprecision(10);
                    sL_file << endl;
                    sL_file << setw(19) << t  << setw(19) << sL;

                }
    }
    else{              //calculate dx_shift using flame speed of species
        int iCO2 = combLine.gas->speciesIndex("CO2");
        int iH2O = combLine.gas->speciesIndex("H2O");
        int iCH4 = combLine.gas->speciesIndex("CH4");
        int iO2  = combLine.gas->speciesIndex("O2");
        int iNO  = combLine.gas->speciesIndex("NO");
        combLine.calculateFlameSpeeds();

        if(neddies % (odtP.modDisp) == 0) *proc.ostrm << endl <<"sL_CH4: "<< combLine.sL[iCH4] <<" sL_O2: "<< combLine.sL[iO2] <<" sL_CO2: "<< combLine.sL[iCO2] <<" sL_H2O: "<< combLine.sL[iH2O]<<" sL_NO: "<< combLine.sL[iNO];

        //        double dx_shift = 1.0*(combLine.sL[iCO2] +combLine.sL[iH2O]+combLine.sL[iO2]+combLine.sL[iCH4]+combLine.sL[iNO])/5.0*(t-t_0);
        double dx_shift = 1.17*(combLine.sL[iCO2]+combLine.sL[iH2O]+combLine.sL[iCH4]+combLine.sL[iO2])/4.*(t-t_0);

        cout <<endl<< "sL_mean:  " << 1.0*(combLine.sL[iCO2] +combLine.sL[iH2O]+combLine.sL[iO2]+combLine.sL[iCH4]+combLine.sL[iNO])/5.0; cout.flush();

        if(dx_shift < combLine.Ldomain/4){ // !!

            int i;
            double phi_in = odtP.phi, temp_in = odtP.tempIn;
            if(odtP.timeInput){
                for(i=0; i< lengthTimeData; i++){
                    if(timedata_time[i]> t) break;
                }
            }
            if(odtP.timeInput && i==0) phi_in = timedata_phi[0];
            else if(odtP.timeInput && i==lengthTimeData-1) phi_in = timedata_phi[lengthTimeData-1];
            else if(odtP.timeInput){
                phi_in = timedata_phi[i-1] + (timedata_phi[i] -timedata_phi[i-1])/(timedata_time[i] -timedata_time[i-1]) * (t-timedata_time[i-1]);}
            else if(odtP.SinusTimeInput && odtP.timeInpSwitch){
                if (odtP.timeInpSwitch == 1)  // phi
                    phi_in = odtP.phi + odtP.SinusAmpl*odtP.phi * sin(2*PI*odtP.sinusFreq*t);
                else if (odtP.timeInpSwitch == 2) 
                    temp_in = odtP.tempIn + odtP.SinusAmpl*odtP.tempIn * sin(2*PI*odtP.sinusFreq*t);
            }
            else  phi_in =odtP.phi;    
            combLine.expand1(dx_shift, 0.0, phi_in, temp_in);


            lineDiffuser.chopOutflowGrid(0.0, odtP.domainLength);
            if(neddies % (odtP.modDisp) == 0)  *proc.ostrm << endl << "dx_shift = " << dx_shift;
        }


        if(iEtrials%50==0 && !odtP.binary){


            sL_file << scientific;
            sL_file << setprecision(10);
            sL_file << endl;
            sL_file << setw(19) << t 
                << setw(19) << combLine.sL[iCH4]
                << setw(19) << combLine.sL[iO2]
                << setw(19) << combLine.sL[iCO2]
                << setw(19) << combLine.sL[iH2O]
                << setw(19) << combLine.sL[iNO];

        }

    }          */       

}




