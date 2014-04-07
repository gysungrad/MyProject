#include "linegroup.h"
#include "odtSolver.h"

using namespace std;

linegroup::linegroup(double inleft,double inright, double intime0, odtParam inodtP, odtSolver *incontrol, IdealGasMix *cantIG, Transport   *cantTran, streams     *strm_p ):
    leftedge(inleft),
    rightedge(inright),
    curtime(intime0), /*,tLastDA(0),cLastDA(0)*//*,odtl(1.1*(inright-inleft)/inodtP.dxmax,inright-inleft),*/
    odtl(*(incontrol->mainLine)),    /*,lastDA(inodtP.sLastDA,0),*/
    odtP(inodtP)    /*,odtStats(inodtP,1,100)*/
{

    odtl.addL(inleft);

//     meshAdapter= adaptMesh(&odtl,&odtP,&odtl.mesherPhi);

    control=incontrol;
}




linegroup::linegroup(const linegroup &ingroup):odtl(ingroup.odtl),/*odtl(ingroup.odtl),*/odtP(ingroup.odtP){
    /*******************************************************************************************************************
     *copy contructor. the address of the odtline changes  is line group specific which is why the copy of              *
     *the mesh adapter needs the new address of the new copied odtline                                                  *
     *******************************************************************************************************************/

    leftedge=ingroup.leftedge;
    rightedge=ingroup.rightedge;
    dtCUmax = ingroup.dtCUmax;
    curtime = ingroup.curtime;

    eddyf=ingroup.eddyf;
    ed=ingroup.ed;


    control = ingroup.control;

//     meshAdapter = ingroup.meshAdapter;
    odtl=ingroup.odtl;
//     meshAdapter.anyl=&odtl;
//     meshAdapter.phi=&odtl.temp;

}







// ///////////////////////////////////////////////////////////////////////////////
void linegroup::adaptEddyRegionOfMesh(int iStart, int iEnd, double &dtCUmax,
        const double &time) {


    double posLower;
    double posUpper;
    int    iStart2, iEnd2;                // index positions for odtl
    double left,right;
    if(iStart > 0)             iStart--;


    posLower = odtl.posf[iStart];   // to update iStart, iEnd below
    posUpper = odtl.posf[iEnd+1];
    left = odtl.posf[iStart];   // to update iStart, iEnd below
    right = odtl.posf[iEnd];
    cout << endl << "1 posLower, posUpper, iStart, iEnd " <<posLower <<" " <<posUpper << " " << odtl.posf[iStart] <<" " <<odtl.posf[iEnd]; cout.flush();

    iStart2 = iStart;
    iEnd2 = iEnd;
    odtl.setTempVec();
    cout << endl << "1 posLower, posUpper, iStart, iEnd " <<posLower <<" " <<posUpper << " " << odtl.posf[iStart2] <<" " <<odtl.posf[iEnd2]; cout.flush();

    odtl.meshAdapter.adaptGrid(iStart2, iEnd2);      
    dtCUmax = computeDtCUmax();

    cout << endl << "2 posLower, posUpper " <<posLower <<" " <<posUpper; cout.flush();
    iStart = odtl.linePositionToIndex(posLower, true); 
    iEnd   = odtl.linePositionToIndex(posUpper, false); 


}

void linegroup::set_odtline(odtline &inputline){

    /***************************************************************************************************
     * when the odtline of the linegroup changes, the mesh adapter needs an update of its new address   *
     ***************************************************************************************************/

    odtl=inputline;
//     //meshAdapter.resetPointers(&odtl, &odtl.temp, &odtl.phase);
//     meshAdapter.anyl=&odtl;
//     if(odtP.Lcombustion) meshAdapter.odtl=&odtl;
//     meshAdapter.phi=&odtl.temp;
//     meshAdapter.bdy=&odtl.phase;
//     meshAdapter.ngrd=odtl.ngrd;
//     meshAdapter.ngrdf=odtl.ngrdf;
}
// // 
double linegroup::computeDtCUmax(){
    //make sure DtCUmax is the highest number possible
    return numeric_limits<double>::max();

}

linegroup& linegroup::operator=(const linegroup &ingroup){

    if(this!=&ingroup){
        leftedge=ingroup.leftedge;
        rightedge=ingroup.rightedge;
        dtCUmax = ingroup.dtCUmax;
        curtime = ingroup.curtime;
        eddyf=ingroup.eddyf;
        ed=ingroup.ed;
        odtP = ingroup.odtP;
        control = ingroup.control;
//         meshAdapter = ingroup.meshAdapter;

        odtl=ingroup.odtl;
//         meshAdapter.anyl=&odtl;
//         meshAdapter.phi=&odtl.temp;
//         meshAdapter.odtP=&odtP;
//         meshAdapter.bdy=&odtl.phase;

    }
    return *this;

}


