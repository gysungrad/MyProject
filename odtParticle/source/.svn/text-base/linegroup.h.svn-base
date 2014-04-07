
/***************************************************************************************
The class line group is the parent class for the classes sub-domain                    *
and sd group and provides all the common function and especially                       *
for the functions for the mesh adaption                                                *
                                                                                       *
***************************************************************************************/
#ifndef LINEGROUP_H
#define LINEGROUP_H
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <limits>
#include "anyline.h"
#include <list>
#include "tripletmap.h"
#include "eddy_functions.h"
#include "eddy.h"
#include "adaptMesh.h"
#include "diffuser.h"
#include "odtParam.h"

class odtSolver;

/**
 * The class line group is the parent class for the classes sub-domain                    
 * and sd group and provides all the common function and especially                       
 * for the functions for the mesh adaption                                                
**/
class linegroup {

public:
    double leftedge,rightedge;
    double dtCUmax;
    double curtime;

    odtline odtl;
    anyline anyl;


    eddy_functions eddyf; 
    eddy        ed;                                  //functions like applying the triplet maps are in the class eddy_functions

    // std::vector<double> lastDA;
								    //time is the time during indepenent running of the subdomaindiffuser sddiffus
    odtParam odtP;
    odtSolver *control;
//     adaptMesh meshAdapter;


    bool apply_maps();                                              //apply all the maps until a dependent maps shows up
    void set_odtline(odtline &inputline);

    //////////////diffusion adaption////////////////////////////////////////////////////////////////////
    void   adaptEddyRegionOfMesh(int iStart, int iEnd, double &dtCUmax,const double &time);
    double computeDtCUmax();
    linegroup& operator=(const linegroup &ingroup);

    ////////////constructors///////////////////////

    linegroup(double inleft, double inright, double intime, odtParam inodtP, odtSolver *incontrol, IdealGasMix *cantIG, Transport   *cantTran, streams     *strm_p);
    linegroup(const linegroup &ingroup);
    ~linegroup(){ control=0;};

};

#endif
