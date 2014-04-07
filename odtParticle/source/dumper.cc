#include "dumper.h"
#include "processor.h"
#include <sstream>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include <cmath>

using namespace std;

extern processor proc;
extern string inputFileDir;

///////////////////////////////////////////////////////////////////////////////////

/** Constructor function
 * @param fname \input name of file contained times to output data at.
 *  fname defaults to "../input/" + inputFileDir + "/dumpTimes.inp"
 */

dumper::dumper(double endTime, std::string fname) {

    iNextTime = 0;
    this->endTime = endTime;

    ifstream ifile(fname.c_str());
    if(!ifile){
        *proc.ostrm << "\n Unable to locate file. \n Where is " << fname << "?" << endl;
        exit(EXIT_FAILURE);
    }

    string s1;
    //double d1; //  !!!!!   unused variable
    
    ifile >> n_time;     getline(ifile, s1);
    dumpTimes.resize(n_time);
    //doldb dumpTimes = linspace(0.0,endTime,n_time);
    //doldb return;
    for(int i=0; i<n_time; i++) {
        ifile >> dumpTimes[i];
        if(ifile.eof()) {
            *proc.ostrm << "\nNot enough times in ../input/" + inputFileDir + "/dumpTimes.inp" << endl;
            exit(EXIT_FAILURE);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////

/** Function returns true as a signal to write the next line property file.
 *  Function queries the list of times.
 * @param currentTime \input interrogation time, from another advancement routine, like diffuser.
 */

bool dumper::checkDumpTime(double currentTime) {

    if(n_time==0)
        return false;

    if(iNextTime < n_time && currentTime >= dumpTimes[iNextTime] ) {
        makeFilename(iNextTime, currentTime);
        *proc.ostrm << "\n# Dumping files at currentTime, dumpTime "
                            << currentTime << " " << dumpTimes[iNextTime] << endl;
        iNextTime++;
        return true;
    }
    else
        return false;
}

///////////////////////////////////////////////////////////////////////////////////

/** Called by checkDumpTime, writes 
 *  Function queries the list of times.
 * @param i           \input file index
 * @param currentTime \input interrogation time, from another advancement routine, like diffuser.
 */

void dumper::makeFilename(int i, double currentTime) {


    stringstream ss1;
    ss1 << scientific;
    ss1 << setprecision(2);
    ss1 << i+1 << ".dat";
    ss1 >> currentTimeFile;

}
////////////////////////////////////////////////////////////////////////////
/**
 * generates n row vector of n linearly equally spaced points between x1 and x2
 * @param x1 start point
 * @param x2 end point
 * @param n points between x1 and x2
 * @return a row vector of n linearly equally spaced points between x1 and x2
 */
vector<double> dumper::linspace(double x1, double x2, int n) {
    vector<double> linSpacedVector(floor(n), 1);
    if (floor(n) <= 0) {
        return linSpacedVector;
    }
    if (n == 1) {
        linSpacedVector[0] = x2;
    } else {
        linSpacedVector[0] = x1;
        for (int i = 1; i < floor(n); i++) {
            double dx = (x2 - x1) / (double(floor(n)) - 1);
            linSpacedVector[i] = x1 + dx*i;
        }
    }
    return linSpacedVector;

}

///////////////////////////////////////////////////////////////////////////////////
