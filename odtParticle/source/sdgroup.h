/*************************************************************************************************************
The class sdgroup is used when multiple sub-domains are merged together. It has it's own                     *
function to applye a triplet map to it and it's own vector which parts of the sd-group is adapted last.      *
This vector is compiled from the vector of the sub-domains.                                                  *
                                                                                                             *
*************************************************************************************************************/
#ifndef SDGROUP_H
#define SDGROUP_H

#include "linegroup.h"
#include <vector>
#include "subdomain.h"
#include <cmath>



class sdgroup: public linegroup{
   public:

    double mapleft,mapright;

    bool apply_maps();

    sdgroup(double inleft,double inright, double intime, double inmapleft, double inmapright, odtParam &inodtP, odtSolver *incontrol,std::vector<int> insd, IdealGasMix *cantIG , Transport   *cantTran, streams     *strm_p );
    sdgroup(const sdgroup &ingroup);
    ~sdgroup(){};

    std::vector<double> return_lastDArange(double lStart, double lEnd,bool wrap=false);
    void shift_group(double inL);

    sdgroup& operator=(const sdgroup &insdg);
    std::vector<int> subdomains;

};

#endif
