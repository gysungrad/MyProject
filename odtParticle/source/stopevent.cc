
#include "stopevent.h"

stopevent::stopevent(double inlStart, double inlEnd, double intime, std::vector<int> insubdomains):
    lStart(inlStart), 
    lEnd(inlEnd),
    subdomains(insubdomains),
    donotgather(false),
    time(intime) {};

