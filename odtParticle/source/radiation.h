/**
 * @file radiation.h
 * Header file for class radiation
 */

#ifndef RADIATION_H
#define RADIATION_H

#include "odtParam.h"
#include "particles.h"
#include <vector>
#include <string>

#ifdef DOCANTERA
#ifndef CANTERA21
#include "Cantera.h"
#endif
#include "IdealGasMix.h"
#endif

#include "cantera_shell_functions.h"

#ifdef CANTERA18
using Cantera_CXX::IdealGasMix;
#else
using Cantera::IdealGasMix;
#endif
using namespace std;

///////////////////////////////////////////////////////////////////////////////

/** Class implementing radiation models: optically thin or two flux.          
 *  Assumes gray gases                                 
 *  
 *  @author David O. Lignell
 */

class radiation {

  ////////////////////// DATA MEMBERS /////////////////////

    public:

        double                       nRadSp;      ///< number of radiating species
        vector<vector<double> >      radCoefs;    ///< Radiation Coefficients [spc][coef]
        int                          Imode;       ///< 1=opthin, 2=twoflux
        vector<int>                  iRadIndx;    ///< radiation species indicies: ch4 co2 h2o co: negative if not present

        double                       sigmaSB;     ///< Stefan Boltzman const

        double                       TloBC;       ///< T on lower Boundary
        double                       ThiBC;       ///< T on upper Boundary;

        double                       sootFactor;  ///< Ksoot = 1863 * fvsoot * T


    ////////////////////// MEMBER FUNCTIONS  /////////////////////

        void getRadHeatSource(const vector<vector<double> > &xMoleSp, 
                              const vector<double>          &temp, 
                              const double                  pressure,
                              const vector<double>          &xPosf,
                              vector<double>                &radSource_G, 
                              vector<double>                &radSource_P, 
                              vector<double>                &ka,
                              particles                     *part, 
                              const vector<double>          &fvSoot = vector<double>(0,0.0));
    private:

        void opthinRadHeatSource(const vector<vector<double> > &xMoleSp, 
                                 const vector<double> &temp, 
                                         const double pressure,
                                 vector<double> &radSource_G, 
                                 vector<double> &ka, 
                                 const vector<double> &fvSoot=vector<double>(0,0.0));

        double getGasAbsorptionCoefficient(const vector<double> &xMole, const double &T,
                                         const double &pressure, const double &fvSoot=0);

        void twoFluxRadHeatSource(const vector<vector<double> > &xMoleSp, 
                                  const vector<double> &temp, 
                                         const double pressure,
                                  const vector<double> &xPosf,
                                        vector<double> &radSource_G, 
                                        vector<double> &radSource_P, 
                                        vector<double> &ka, 
                                        particles *part,
                                  const vector<double> &fvSoot=vector<double>(0,0.0));

    ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:

        radiation(odtParam &odtpp, IdealGasMix *gas=0);   // constructor
        radiation(){}   // constructor



};

#endif
