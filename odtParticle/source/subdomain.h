/*****************************************************************************************************************
The class sub-domain represents the sub-domain, it inherits from odtline,                                        *
represents a subdomain and has a list of maps which are executed by the subdomain without any regard             *
for other subdomains. The subdomain has a solution region, from guardleft to guardright, all                     *
the solution regions together are the solution.                                                                  *
The subdomain has global positions, i.e. subdomain.pos[i] gives the global position on the mainline.             *
update the odtline only through the function provided for it                                                     *
*****************************************************************************************************************/

#ifndef SUBDOMAIN_H
#define SUBDOMAIN_H
#include <iostream>
#include <vector>
#ifdef BOOST
#include <boost/ptr_container/ptr_vector.hpp>
#endif
#include <fstream>
#include <string>
#include <limits>
#include "anyline.h"
#include <list>
#include "tripletmap.h"
#include "eddy_functions.h"
#include "eddy.h"
// #include "adaptMesh.h"
#include "diffuser.h"
#include "linegroup.h"
#include <omp.h>
#include "walltime.h"
//#include "odtSolver.h"

class odtSolver;


/**
* The class sub-domain represents the sub-domain, it inherits from odtline,                                        
* represents a subdomain and has a list of maps which are executed by the subdomain without any regard            
* for other subdomains. The subdomain has a solution region, from guardleft to guardright, all                   
* the solution regions together are the solution. \n
* The subdomain has global positions, i.e. subdomain.pos[i] gives the global position on the mainline. \n
* update the odtline only through the function provided for it                                                
*/

class subdomain : public linegroup {

public:
//////////Constructors/////////////////////////////////////////////////

subdomain(double inleftedge, double inrightedge, double inguardleft, double inguardright, double intime0, int insdnumber,odtParam &inPar, odtSolver *incontrol, IdealGasMix *cantIG, Transport   *cantTran, streams     *strm_p );
//subdomain(double inleftedge, double inrightedge);
subdomain(const subdomain &insd);
~subdomain(){
    //cout << "sd D  &sd : " << this<<endl;cout.flush();
};

std::list<tripletmap>    maps;



double guardleft,guardright;                                         //the borders of the subdomain in real space
double timelimit,t0left,t0right;                                    //time0 is the time when the subdomain began running indepentendly, i.e. time of last synchronisation
int printcounter;                                                   //time is the time during indepenent running of the subdomaindiffuser sddiffus
int stopnumber;
int sdnumber;
int nmapsindependent;

diffuser sddiffuser;

walltime wt;

double tsaverage;

bool apply_maps();

void diffuse(double intime);
void set_odtline(odtline &inodtline);

void set_subdomain(double inleftedge, double inrightedge, double inguardleft, double inguardright, double intime0,int insdnumber,odtParam *inPar, odtSolver *incontrol, IdealGasMix *cantIG, Transport   *cantTran, streams     *strm_p);


void print_subdomain();

void setTimeLimit(double dx_exp, odtParam &inPar);

//subdomain& operator=(const subdomain &insd);
void operator=(const subdomain &insd);



};


#endif
