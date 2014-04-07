/**
 * @file randomGenerator.h
 * Header file for classes randomGenerator 
 */

#ifndef DOMTRAND

#ifndef RRRNG_H
#define RRRNG_H

#include <cstdlib>
#include "processor.h"

extern processor proc;

////////////////////////////////////////////////////////////////////////////////////

/** A random number generator class.  There are two versions, one of which is 
 *  selected at compile time.  This version just uses the std built in generator.
 *  The other version uses a custom generator (Mersenne twister).
 *  
 *  @author David O. Lignell
 */

class randomGenerator {

    private :

    int seed;                       ///< random number seed 
    double maximum_rand;            ///< used to normalize 0 to 1
    
    public :

    /** Returns the next random number
     *  @returns a normalized random number between 0 and 1
     */
    inline double getRand() {
        return static_cast<double>(rand()) / maximum_rand;
    }
    
    /**
     * Constructor that allows user to set the seed
     *
     * @param aseed \input seed to use. When aseed < 0, a random seed value will be used.
     */
    randomGenerator(int aseed) {
        seed = aseed;
        maximum_rand = static_cast<double> (RAND_MAX);
        if(seed < 0)
            srand(time(NULL)+proc.myid*1000);
        else
            srand(seed);
    }

    /**
     * Default constructor. Uses the seed '25001'.
     */
    randomGenerator() {
        seed = 25001;
        maximum_rand = static_cast<double> (RAND_MAX);
        srand(seed);
    }
};


#endif

#endif


///////////////////////////////////////////////////////////////////////////////////////

#ifdef DOMTRAND

#ifndef RRRNG_H
#define RRRNG_H

#include "Mersenne/MersenneTwister.h"

#include "processor.h"

extern processor proc;

/** A random number generator class.  There are two versions, one of which is 
 *  selected at compile time.  This version sets up and calls the Mersenne twister.
 *  The other version uses the standard c library generator.
 */

///////////////////////////////////////////////////////////////////////////////////////

class randomGenerator {

    private :

    MTRand mtwist;                   ///< Mersenne twister object

    public :

    inline double getRand() {
        return mtwist.rand();
    }
    
    randomGenerator(int aseed) : mtwist(aseed) {
        if(aseed < 0)               // randomize the seed
            mtwist.seed();
            //mtwist.seed(aseed+proc.myid+12800);    
            //mtwist.seed(proc.myid-aseed*1000);    
    }
    randomGenerator()          : mtwist() {}
};


#endif

#endif
