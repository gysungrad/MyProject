/****************************************************************************************************************************
This class saves the start and end point as also which type of stop-event it is (if it is an event with or                  *
without the triplet map. For certain events, it is not necessary that the solution regions are gathered                     *
(for example if it is just a marker event to return after a certain time has passed, e.g. for the gather of statistics      *
                                                                                                                            *
*****************************************************************************************************************************/

#ifndef STOPEVENT_H
#define STOPEVENT_H

#include <vector>
#include "subdomain.h"

class stopevent{
 public:
 
    ~stopevent(){ }
    
    
    stopevent(double inlStart, double inlEnd, double intime, std::vector<int> insubdomains);

    double lStart,lEnd;
    int type;                                                  // 1 for update event, 2 for dependent map event
    std::vector<int> subdomains;
    double time;
    bool applymap;
    bool donotgather;



};

#endif
