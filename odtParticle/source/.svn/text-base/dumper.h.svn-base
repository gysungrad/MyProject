#ifndef DUMPER_H
#define DUMPER_H

#include "anyline.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

extern string inputFileDir;

///////////////////////////////////////////////////////////////////////////////

/**
 * @class dumper
 * 
 * Class managing output of line properties at specified times.
 *  
 *  @author Devin Rappleye
 */

class dumper {
    
    public:

        ////////////////////////// DATA MEMBERS

        vector<double> dumpTimes;    ///< File of times for data dumps
        int n_time;                  ///< Vector size indicator for dumptimes
        double endTime;              ///< endTime for the simulation
        string currentTimeFile;      ///< Names the file that will be dumped

        int iNextTime;               ///< Index of next time in list for output 
       
        ////////////////////////// MEMBER FUNCTIONS 
        
        void readTimeFile(string fname);
        void makeFilename(int i, double currentTime);
        bool checkDumpTime(double currentTime);

        ////////////////////////// CONSTRUCTORS

        dumper(double endTime, string fname="../input/" + inputFileDir + "dumpTimes.inp");
        vector<double> linspace(double x1,
            double x2,
            int n);
       
};

#endif

