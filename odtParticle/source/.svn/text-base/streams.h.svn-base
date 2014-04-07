/**
 * @file streams.h
 * Header file for classes streams 
 */

#ifndef STREAMS_H
#define STREAMS_H

#ifdef DOCANTERA
#ifndef CANTERA21
#include "Cantera.h"
#endif
#include "IdealGasMix.h"
#endif

#include "cantera_shell_functions.h"
#include "inputFile.h"

#ifdef CANTERA18
using Cantera_CXX::IdealGasMix;
#else
using Cantera::IdealGasMix;
#endif

using namespace std;

extern string inputFileDir;
////////////////////////////////////////////////////////////////////////////////

/** Class implementing streams for use in mixing and or reaction problems.
 *  This is writting in terms of mixture fraction with streams defined in an
 *  input file.  The class can implement products of complete combustion, or
 *  equilibrium (through the Cantera IdealGasMix object, if desired).
 *  This class holds a pointer to a Cantera IdealGasMix object (defined up front
 *  in main) that computes thermodynamic, kinetic, and transport data.
 *  
 *  @author David O. Lignell
 */

class streams {

    public:

    //////////////////// DATA MEMBERS //////////////////////
        
        double              T0;              ///< stream mixf=0 temperature
        double              T1;              ///< stream mixf=1 temperature
        vector<double> y0;              ///< stream mixf=0 composition vector
        vector<double> y1;              ///< stream mixf=1 composition vector

        double              h0;              ///< stream mixf=0 enthalpy
        double              h1;              ///< stream mixf=1 enthalpy
        double              pres;            ///< stream pressure (system)

	    IdealGasMix        *gas;             ///< pointer to a (just one) gas object

        int                 nspc;            ///< number of chemical species
        double              mixfStoic;       ///< stoichiometric mixture fraction

                                             /// \fun{\text{mixf} = \frac{\beta-\beta_0}{\beta_1-\beta_0}}
        double              beta0;           //< mixf = (beta-beta0) / (beta1-beta0)
                                             /// \fun{\text{mixf} = \frac{\beta-\beta_0}{\beta_1-\beta_0}}
        double              beta1;           //< mixf = (beta-beta0) / (beta1-beta0)

	int		    comp2;	    ///< for premixed combustion to distinguish between different input possibilities for the mixture


// Members for new input file stuff

        vector<string>  y0Labels;  ///< Composition vector labels (chemical sepcies names), stream 0/1
        vector<string>  y1Labels;  ///< Composition vector labels (chemical sepcies names), stream 0/1

        inputFile           streamInput_;    ///< inputFile object for stream input file


    //////////////////// MEMBER FUNCTIONS /////////////////

        void readStreams(string fname="../input/" + inputFileDir + "streams.inp");
        void getProdOfCompleteComb(double mixf, vector<double> &ypcc, 
                                   double &hpcc, double &Tpcc, int probType);
        double getMixtureFraction(double *y, bool doBeta01=false);

    private:

        void readStreamComp(ifstream &ifile, vector<double> &xy);
        void setStreamComp(vector<double> &y, vector<string> &yLabel);
        void setStoicMixf();
        vector<double> setElementMassFracs(double *y);
        vector<double> setElementMoleFracs(double *y);
        vector<double> getElementMoles(double *x, double &nOnotFromO2,
                                                       double &nHnotFromH2O,
                                                       double &nCnotFromCO2);

    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        streams() {}
        streams(IdealGasMix *cantIG);

        ~streams(){ gas=0; }
        



};



#endif
