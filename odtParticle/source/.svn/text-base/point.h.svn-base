/***************************************************************************************
*A point is a point for the list kept for every guard regions to check if              *
*triplet maps have already covered most of the guard region. It's only                 *
*function is a stuct which make a binary comparision between it's values.              *
*It also has a boolean which is true if the point is at the beginning of a segment     *
*                                                                                      *
***************************************************************************************/
#ifndef POINT_H
#define POINT_H
#include<iostream>
#include <cstdlib>
#include <functional>

class point{

public:

 point (double pos, bool inisbeginning);
 point ();
 ~point(){};
 
 double pos;

 bool isbeginning;

 struct greater:public  std::binary_function<point,point,bool> {

  bool operator()(const point &a1, const point &a2){
  return a1.pos>a2.pos;
  }

};


};


#endif
