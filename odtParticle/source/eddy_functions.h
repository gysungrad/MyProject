
/***************************************************************************************************************************
This is a reduced version of the class eddy.cc it has only a the function apply_single map, which applies some sort of     *
triplet maps to it.                                                                                                        *
***************************************************************************************************************************/
#ifndef EDDYFUNCTIONS_H
#define EDDYFUNCTIONS_H

#include <iostream>
#include <vector>
#include <fstream>
#include "anyline.h"
#include "odtParam.h"
#include "randomGenerator.h"

class eddy_functions {
public:
void apply_single_map(anyline &anyl,double leftEdge, double rightEdge);
~eddy_functions(){}
};

#endif
