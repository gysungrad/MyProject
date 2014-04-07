#ifdef BOOST

#include "subdomain.h"
#include <iomanip>
#include <cmath>          // fabs
#include <sstream>
#include "odtSolver.h"



using namespace std;
using namespace boost;

subdomain::subdomain(double inleftedge, double inrightedge, double inguardleft, double inguardright, double intime0, int insdnumber,
        odtParam &inPar, odtSolver *incontrol, IdealGasMix *cantIG, Transport   *cantTran, streams     *strm_p ):
    linegroup(inleftedge, inrightedge,intime0,inPar,incontrol, cantIG, cantTran, strm_p),
    sddiffuser(/*&odtl,*/&odtl,inPar, 0)

{


    guardleft=inguardleft; guardright=inguardright;

    stopnumber=-1;
    sdnumber=insdnumber;
    t0left=intime0;
    t0right=intime0;
    printcounter=0;
    tsaverage=0;



    //set the maximum time the subdomain can run independently.
    //it's the diffusion time for half the length of the guard regions.
    //make sure the guardsize is the actual guard size in case of non-periodic
    //boundary conditions
    double guardsize=max(guardleft-leftedge,rightedge-guardright);
    double D, fs;
    
    double domainsize;
    
    domainsize= rightedge-leftedge;
    
    D=10*inPar.Dt*min(1.,pow(domainsize/inPar.Lmax, 4/3)); // max diffusivity within the subdomain
    if(odtl.probType==6) 
	fs =10; // max flame speed
    else 
	fs = 0;
        
    timelimit = 0.5* pow(guardsize,2)/(D+fs*guardsize);  // Guardsize² / (diffusivity-Flame speed*Guardsize); 
    
    control=incontrol;

    sddiffuser.subLeft_bc=0;
    sddiffuser.subRight_bc=0;
    
    switch(inPar.bcType){
      
      case 2:
	if(sdnumber==0)
	  sddiffuser.subLeft_bc=1;
	if(sdnumber==incontrol->Nsubdomains-1)
	  sddiffuser.subRight_bc=1;
      case 4:
	if(sdnumber==0)
	  sddiffuser.subLeft_bc=1;
      case 5:
	if(sdnumber==0)
	  sddiffuser.subLeft_bc=1;
    }
    
    //     double guardsize=max(guardleft-leftedge,rightedge-guardright);
    //     timelimit=pow(guardsize,2)/inPar.visc_0/2;
    //     control=incontrol;

}


void subdomain::set_subdomain(double inleftedge, double inrightedge, double inguardleft, double inguardright, double intime0, int insdnumber, odtParam *inPar, odtSolver *incontrol, IdealGasMix *cantIG, Transport   *cantTran, streams     *strm_p ) {

    //linegroup(inleftedge, inrightedge,intime0,inPar,incontrol, cantIG, cantTran, strm_p);

    leftedge = inleftedge;      // init in linegroup constructor
    rightedge = inrightedge;    // init in linegroup constructor
    curtime = intime0;          // init in linegroup constructor
    

    
    odtl=odtline(static_cast<int>(1.1*(inrightedge-inleftedge)/inPar->dxmax),inrightedge-inleftedge, cantIG, cantTran, strm_p, inPar, "", incontrol->LhasVel, odtP.Lrxn);    /*,lastDA(inodtP.sLastDA,0),*/
    odtl.addL(inleftedge);    

//     meshAdapter= adaptMesh(&odtl,&odtl,&inPar,&odtl.temp);

//    sddiffuser = diffuser(&odtl,inPar); // check this ! D.N.
    sddiffuser.set_diffuser(&odtl); // doesn't seem to have any effect! D.N.

    guardleft = inguardleft;
    guardright = inguardright;



    stopnumber=-1;
    sdnumber = insdnumber;
    t0left = intime0;
    t0right = intime0;
    printcounter=0;
    tsaverage=0;
    //set the maximum time the subdomain can run independently.
    //it's the diffusion time for half the length of the guard regions.
    //make sure the guardsize is the actual guard size in case of non-periodic
    //boundary conditions
    double guardsize=max(guardleft-leftedge,rightedge-guardright);
    double D, fs;
    double domainsize;
    
    domainsize= rightedge-leftedge;
    
    D=10*inPar->Dt*min(1.,pow(domainsize/inPar->Lmax, 4/3)); // max diffusivity within the subdomain
    if(odtl.probType==6) 
	fs =10; // max flame speed
    else 
	fs = 0;
        
    timelimit = pow(guardsize,2)/(D+fs*guardsize);  // Guardsize² / (diffusivity-Flame speed*Guardsize); 
    control=incontrol;

}


//subdomain::subdomain(){}

subdomain::subdomain(const subdomain &insd):linegroup(insd),tsaverage(0),sddiffuser(insd.sddiffuser){
    //void subdomain::operator=(const subdomain &insd): linegroup(insd),    tsaverage(0),    sddiffuser(insd.sddiffuser) {

    maps=insd.maps;
    guardleft=insd.guardleft;
    guardright=insd.guardright;
    timelimit=insd.timelimit;
    t0left=insd.t0left;
    t0right=insd.t0right;
    //sddiffuser=insd.sddiffuser;
    stopnumber=insd.stopnumber;
    sdnumber = insd.sdnumber;
    printcounter=insd.printcounter;



}



bool subdomain::apply_maps(){
    /******************************************************************************************************************************************
     * apply tripletmaps until hitting a stop event                                                                                            *
     * ----------------------------------------------------------------------------------------------------------------------------------------*
     * -------------------------------------------------Output---------------------------------------------------------------------------------*
     * returns false if apply_maps stops due to a stop event, true if it runs out of maps                                                      *
     * ------------------------------------------------Modifies--------------------------------------------------------------------------------*
     * changes odtl, deletes maps from list                                                                                                    *
     * ----------------------------------------------------------------------------------------------------------------------------------------*
     * Do simple maps until reaching a stop event. Upon hitting one, check if all the other sub-domain involved in the event                   *
     * are ready. If so, call the controlers function to update deal with the event, e.g. apply the dependent map or update the guard regions. *
     *                                                                                                                                         *
     ******************************************************************************************************************************************/

    
    
    if(stopnumber>=0){

        return false;
    }

    vector<int> sdindex;
    vector<int>::iterator it;
    list<double>::iterator lit;
    bool notallready=false;
    bool found;
    nmapsindependent=0;
    tripletmap tempmap;
    int mapstart,mapend;
    double t1,t2,clockzero;
    t2=0;
    clockzero=0;


    dtCUmax  = computeDtCUmax();


    
   bool breakcon=false;
    while(!maps.empty()){

        nmapsindependent++;
        tempmap=maps.front();
        maps.pop_front();


	

        t1=wt.walltime_return(&clockzero);

        if(tempmap.stopnumber>=0){
#pragma omp critical (crit1)
            {
                found=false;
                for(lit=control->acquired_stopnumbers.begin();lit!=(control->acquired_stopnumbers.end());lit++){
                    if((*lit)==tempmap.stopnumber){
                        maps.push_front(tempmap);
                        found=true;
                    }
                }
                if(!found){
                    control->acquired_stopnumbers.push_back(tempmap.stopnumber);
                }
            }
            if(!found){
                stopnumber=tempmap.stopnumber;
                sdindex=tempmap.subdomains;

                for(it=sdindex.begin(); it!= sdindex.end(); it++){
                    if(!(control->subdomains_comb.at(*it).stopnumber==stopnumber)){
                        notallready=true;
                        break;
                    }
                }
         

	        

                if(!notallready){
#pragma omp critical (apply_event)
                    {
                        control->apply_event(stopnumber);
                    }
   cout << endl << "  apply maps "; 

                }
#pragma omp critical (crit1)
                {
                    control->acquired_stopnumbers.remove(tempmap.stopnumber);

                }
		
                odtl.meshAdapter.adaptGrid(0, odtl.ngrd-1); 
		

            }
            control->timeinsd+=t2;
            cout << endl << "curtime: " <<curtime; cout.flush();

            return false;
        }
       
	
	
        t1=wt.walltime_return(&clockzero);
        control->eddycount++;
	 		
        mapstart=odtl.linePositionToIndex(tempmap.lStart,true);

        mapend=odtl.linePositionToIndex(tempmap.lEnd,false);

        odtline segment(&odtl, mapstart,mapend, (mapstart>mapend), false);
	
		
        eddyf.apply_single_map(segment,tempmap.lStart,tempmap.lEnd);
        
	
        odtl.insertEddy(segment, mapstart,mapend,tempmap.lStart,tempmap.lEnd, (mapstart > mapend));

	
        adaptEddyRegionOfMesh(mapstart, mapend, dtCUmax,curtime);

	    

        t2+=wt.walltime_return(&t1);
        // 	meshAdapter.adaptGrid(0, odtl.ngrd-1); 

        //mapstart and mapend need not to be defined, they are not used but overwritten
        //in adaptAfterSufficientDiffTime.
        //LEM, this is useless anyway since dtCUmax is set to be so high, this function actualy does
        //nothing

					    

        diffuse(tempmap.time);
        // adaptAfterSufficientDiffTime(curtime, tLastDA, cLastDA,mapstart, mapend, dtCUmax);
        //t2+=wt.walltime_return(&t1);

	
        curtime=tempmap.time;
       

    }


    control->timeinsd+=t2;

    cout << endl << "curtime: " <<curtime; cout.flush();
    //if(sdnumber == 5) print_subdomain();
    odtl.meshAdapter.adaptGrid(0, odtl.ngrd-1); 



    return true;

}


// 
void subdomain::diffuse(double intime){
    /***********************************************************************************************************************
     * diffues odtl until time intime                                                                                       *
     * ---------------------------------------------------------------------------------------------------------------------*
     * ---------------------------------------------------Modifies----------------------------------------------------------*
     * calls the diffuser which will modify the odtl                                                                        *
     * ---------------------------------------------------------------------------------------------------------------------*
     *                                                                                                                      *
     ***********************************************************************************************************************/
    if(sdnumber < control->Nsubdomains-1) {
        sddiffuser.indNextSub_start = odtl.linePositionToIndex(sddiffuser.posNextSub_start,true);

        //     vector<double> icp(1);
        //     icp[0] = sddiffuser.posNextSub_start- odtl.posf[sddiffuser.indNextSub_start];
        //     odtl.splitCell(sddiffuser.indNextSub_start, 1,icp , odtP, false);
        //     sddiffuser.indNextSub_start++;
        //  
        // //    cout << endl << "1: dx1 " << odtl.posf[sddiffuser.indNextSub_start-1]-odtl.posf[sddiffuser.indNextSub_start-2]<< " dx2 " << odtl.posf[sddiffuser.indNextSub_start]-odtl.posf[sddiffuser.indNextSub_start-1]<< " dx2 " << odtl.posf[sddiffuser.indNextSub_start+1]-odtl.posf[sddiffuser.indNextSub_start];
        //       odtl.setTempVec();
        // //       if(  (fabs(odtl.temp[sddiffuser.indNextSub_start-1]-odtl.temp[sddiffuser.indNextSub_start-2])- odtP.largeGradFrac*(1500) < 1.0e-8) && ((odtl.posf[sddiffuser.indNextSub_start]-odtl.posf[sddiffuser.indNextSub_start-2])- odtP.dxmax < 1.0e-8)){
        // 	odtl.merge2cells(sddiffuser.indNextSub_start-2);
        // 	sddiffuser.indNextSub_start--;
        // 
        // //      }
        // //       if(  ((fabs(odtl.temp[sddiffuser.indNextSub_start+1]-odtl.temp[sddiffuser.indNextSub_start])- odtP.largeGradFrac*(1500)) < 1.0e-8) && (((odtl.posf[sddiffuser.indNextSub_start+2]-odtl.posf[sddiffuser.indNextSub_start])- odtP.dxmax) < 1.0e-8))
        // 	odtl.merge2cells(sddiffuser.indNextSub_start);
        // //    cout << endl << "2: dx1 " << odtl.posf[sddiffuser.indNextSub_start-1]-odtl.posf[sddiffuser.indNextSub_start-2]<< " dx2 " << odtl.posf[sddiffuser.indNextSub_start]-odtl.posf[sddiffuser.indNextSub_start-1]<< " dx2 " << odtl.posf[sddiffuser.indNextSub_start+1]-odtl.posf[sddiffuser.indNextSub_start];*/
        // 	
        //      }
        //       if((odtl.posf[sddiffuser.indNextSub_start]-odtl.posf[sddiffuser.indNextSub_start-1]-(odtl.posf[sddiffuser.indNextSub_start+1]-odtl.posf[sddiffuser.indNextSub_start])< 1.0e-8) && ((odtl.posf[sddiffuser.indNextSub_start]-odtl.posf[sddiffuser.indNextSub_start-2]-odtP.dxmax)< 1.0e-8)  && (fabs(odtl.temp[sddiffuser.indNextSub_start-1]-odtl.temp[sddiffuser.indNextSub_start-2])- odtP.largeGradFrac*(1500) < 1.0e-8)){
        // 	odtl.merge2cells(sddiffuser.indNextSub_start-2);
        // 	sddiffuser.indNextSub_start--;
        // 
        //       }
        //       else if((odtl.posf[sddiffuser.indNextSub_start]-odtl.posf[sddiffuser.indNextSub_start-1]-(odtl.posf[sddiffuser.indNextSub_start+1]-odtl.posf[sddiffuser.indNextSub_start])> 1.0e-8) && ((odtl.posf[sddiffuser.indNextSub_start+2]-odtl.posf[sddiffuser.indNextSub_start] - odtP.dxmax) < 1.0e-8) && (fabs(odtl.temp[sddiffuser.indNextSub_start+1]-odtl.temp[sddiffuser.indNextSub_start])- odtP.largeGradFrac*(1500) < 1.0e-8)){
        // 	odtl.merge2cells(sddiffuser.indNextSub_start);
        // 
        //       }
        //       else if(((odtl.posf[sddiffuser.indNextSub_start]-odtl.posf[sddiffuser.indNextSub_start-1]-odtP.dxmin)< 1.0e-8) ){// || (fabs(odtl.temp[sddiffuser.indNextSub_start-1]-odtl.temp[sddiffuser.indNextSub_start-2])- odtP.smallGradFrac*(1500) < 1.0e-8)){
        // 	odtl.merge2cells(sddiffuser.indNextSub_start-2);
        // 	sddiffuser.indNextSub_start--;
        // 
        //       }
        //       else if(((odtl.posf[sddiffuser.indNextSub_start+1]-odtl.posf[sddiffuser.indNextSub_start]-odtP.dxmin)< 1.0e-8) ){//  || (fabs(odtl.temp[sddiffuser.indNextSub_start+1]-odtl.temp[sddiffuser.indNextSub_start])- odtP.smallGradFrac*(1500) < 1.0e-8)){
        // 	odtl.merge2cells(sddiffuser.indNextSub_start);
        //       
        //       }
        //       else if((fabs(odtl.temp[sddiffuser.indNextSub_start-1]-odtl.temp[sddiffuser.indNextSub_start-2])- odtP.smallGradFrac*(1500) < 1.0e-8)){
        // 	odtl.merge2cells(sddiffuser.indNextSub_start-2);
        // 	sddiffuser.indNextSub_start--;
        //       }
        //       else if((fabs(odtl.temp[sddiffuser.indNextSub_start+1]-odtl.temp[sddiffuser.indNextSub_start])- odtP.smallGradFrac*(1500) < 1.0e-8)){
        // 	odtl.merge2cells(sddiffuser.indNextSub_start);
        //       
        //       }
    }
    else sddiffuser.indNextSub_start = odtl.ngrd -1;

    double dtdiffuse=intime-curtime;

    if(dtdiffuse > 0) {
        sddiffuser.diffuseOdtLine(dtdiffuse,curtime);
    }
    tsaverage+=sddiffuser.numberofsteps;

}




void subdomain::print_subdomain(){
    stringstream ss1;  string s1;
    ss1 << "sd" << sdnumber << "_" << printcounter; ss1 >> s1;
    odtl.outputProperties("../data/"+s1);
    printcounter++;


}



void subdomain::set_odtline(odtline &inodtline){

    linegroup::set_odtline(inodtline);

    sddiffuser.odtl = &odtl;
    sddiffuser.brxr.odtl = &odtl;
}

void subdomain::operator=(const subdomain &insd){
    if(this != &insd){
        linegroup::operator=(insd);
        maps=insd.maps;
        guardleft=insd.guardleft;
        guardright=insd.guardright;
        timelimit=insd.timelimit;
        t0left=insd.t0left;
        t0right=insd.t0right;
        sddiffuser=insd.sddiffuser;
        stopnumber=insd.stopnumber;
        sdnumber = insd.sdnumber;
        printcounter=insd.printcounter;
    }
    // return *this;


} 
/*
   subdomain& subdomain::operator=(const subdomain &insd){
   cout << endl << "subdomain copy constr 2" << endl; cout.flush();
   if(this != &insd){
   linegroup::operator=(insd);
   maps=insd.maps;
   guardleft=insd.guardleft;
   guardright=insd.guardright;
   timelimit=insd.timelimit;
   t0left=insd.t0left;
   t0right=insd.t0right;
   sddiffuser=insd.sddiffuser;
   stopnumber=insd.stopnumber;
   sdnumber = insd.sdnumber;
   printcounter=insd.printcounter;
   }
   return *this;


   } */

void subdomain::setTimeLimit(double dx_exp, odtParam &inPar){
  
  
  double guardsize=max(guardleft-leftedge,rightedge-guardright);
  double D, fs;
    
  double domainsize;
  
  domainsize= rightedge-leftedge;
  
  D=10*inPar.Dt*min(1.,pow(domainsize/inPar.Lmax, 4/3)); // max diffusivity within the subdomain
  if(odtl.probType==6) 
      fs =10; // max flame speed
  else 
      fs = 0;
      
  timelimit = 0.5* pow(guardsize,2)/(D+fs*guardsize);
  
  timelimit = min(timelimit, guardsize/dx_exp * inPar.dtGatherSubdomains * 0.03);
  cout << endl << "timelimit " << timelimit;
}

#endif
























