#ifndef TRIPLETMAP_H
#define TRIPLETMAP_H
#include "odtline.h"
#include <cstdlib>
#include <functional>

class subdomain;

/**
 * This class saves the location and the time of a triplet map. 
 * Also, it has a function to give order to a triplet map. 
 * The criteria is the time, if the time is equal the second criteria 
 * is the stop number it has.   \n 
 * triplet map with different stop numbers but the same time should never
 * be equal in terms of the opererator greater   \n
 * Do not include the sub-domain header file in here, this will lead to a circular dependency. 
 */
class tripletmap{

public:

tripletmap(double lStart, double lEnd, double intime, odtline &inodtline);
tripletmap(double lStart, double lEnd, double intime);

tripletmap();
~tripletmap(){};

double lStart,lEnd;

double time;

bool applymap;

int stopnumber;

//debug
int type; //0 undefined/simplemap, 1 stopevent, 2 time update 3 wrap map


 struct greater:public  std::binary_function<tripletmap,tripletmap,bool> {

  bool operator()(const tripletmap &a1, const tripletmap &a2){
      if(a1.time!=a2.time){
    return a1.time>=a2.time;
      }
      else{
       return a1.stopnumber>a2.stopnumber;
      }
  }
 };


std::vector<int> subdomains;




};

#endif
