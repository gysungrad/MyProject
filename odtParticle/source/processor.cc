/**
 * @file processor.cc
 * Source file for class processor
 */

#include "processor.h"
#include <sys/stat.h>             // for mkdir
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

int processor::nInst;

///////////////////////////////////////////////////////////////////////////////

/** Constructor */

processor::processor() {

#ifdef DOMPI

    //----------- set MPI Stuff (if on)

    if((++nInst)==1) {                 // Only ever call MPI_* once
        int fake_argc = 0;
        char** fake_argv;
        MPI_Init(&fake_argc, &fake_argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &myid);
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    }
    if(nInst > 1) 
        cout << endl << "****** WARNING, more than one processor class object" << endl;

    if(myid==0) 
        cout << "\n# MPI IS ON" << "; Nprocs = " << nproc << endl;
    MPI_Barrier(MPI_COMM_WORLD);

    //----------- set the data directory

    stringstream ss1;  
    string s1;         
    ss1.clear(); ss1 << myid; ss1 >> s1;
    dataDir = "../data/data_" + s1 + "/";   // e.g., "../data_01", etc.

    //----------- create the data directory

    int iflag = mkdir(dataDir.c_str(), 0755);
    if(iflag != 0)  {
        cout << "\n***** Error, process " << myid << " failed to create " 
             << dataDir << ", or it was already there" << endl;
        exit(0);
    }

    //----------- set the runtime file

    ss1.clear(); ss1 << myid; ss1 >> s1; 
    string fname = "./RUNTIME/runtime_" + s1;
    ostrm = new ofstream(fname.c_str());
    // write with *ostrm << "stuff" or *proc.ostrm << "stuff"

#else

    dataDir = "../data/";

    ostrm = &cout;

    myid  = 0;
    nproc = 1;

#endif

}


///////////////////////////////////////////////////////////////////////////////

/** Destructor */

processor::~processor() {

#ifdef DOMPI

    if((--nInst)==0)             // Only ever finalize mpi once
        MPI_Finalize();

    delete ostrm; 

#endif
}

///////////////////////////////////////////////////////////////////////////////





