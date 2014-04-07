#include "sdgroup.h"
#include "odtSolver.h"

using namespace std;

sdgroup::sdgroup(double inleft,double inright, double intime, double inmapleft, double inmapright, odtParam &inodtP, odtSolver *incontrol, vector<int> insd, IdealGasMix *cantIG, Transport   *cantTran, streams     *strm_p  ):linegroup(inleft,inright,intime,inodtP,incontrol,cantIG,cantTran,strm_p),mapleft(inmapleft),mapright(inmapright),subdomains(insd){}

// sdgroup::sdgroup(){}

sdgroup::sdgroup(const sdgroup &insdg):linegroup(insdg){
mapleft=insdg.mapleft;
mapright=insdg.mapright;
subdomains=insdg.subdomains;
}


bool sdgroup::apply_maps(){


    int i1,i2;


   if(mapleft<odtl.posf[0]){
        mapleft+=odtP.domainLength;
        mapright+=odtP.domainLength;
    }

    control->eddycount++;
    i1 = odtl.linePositionToIndex(mapleft,true);
    i2 = odtl.linePositionToIndex(mapright,false);

    odtline eddyline(&odtl,i1,i2,i1>i2,false);

    eddyf.apply_single_map(eddyline,mapleft,mapright);

    odtl.insertEddy(eddyline,i1,i2,mapleft,mapright,false);

    adaptEddyRegionOfMesh(i1,i2,dtCUmax,curtime/*,tLastDA,cLastDA*/);


    i1 = odtl.linePositionToIndex(mapleft,true);
    i2 = odtl.linePositionToIndex(mapright,false);

    return true;


}




void sdgroup::shift_group(double inL){
 

    odtl.addL(inL);
    leftedge+=inL;
    rightedge+=inL;
    mapleft+=inL;
    mapright+=inL;

}

sdgroup& sdgroup::operator=(const sdgroup &insdg){
    if(this!=&insdg){
        linegroup::operator=(insdg);
        mapleft=insdg.mapleft;
        mapright=insdg.mapright;
        subdomains=insdg.subdomains;
    }
    return *this;
}
