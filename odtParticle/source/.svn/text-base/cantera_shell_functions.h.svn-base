/**
 * @file cantera_shell_functions.h
 * Classes for Cantera shell functions
 */

#ifndef DOCANTERA
#ifndef CANTERA_SHELL_FUNCTIONS_H
#define CANTERA_SHELL_FUNCTIONS_H

// This file is used to store dummy function definitions when cantera is not available

#include <string>
#include <vector>

using std::string;
using std::vector;

namespace Cantera {

//-------------------------------------------------------------------------------

/** Cantera shell functions when not compiling with Cantera.  
 *  When cantera is not needed, compile with these
 *  routines that have the same headers as those of Cantera.  Add new routines
 *  when new cantera functions are called.  This also provides a template for
 *  implementing user-defined functions instead of calling Cantera.  However,
 *  the routines are not meant to be called, and return errors if they are.
 *  
 *  @author David O. Lignell
 */

class IdealGasMix {

    public :

        void           setState_PY(double p, double* y);
        void           setState_HP(double h, double p, double tol = 1.e-8);
        double         temperature();
        double         density();
        double         cp_mole();
        double         cp_mass();
        void           getEnthalpy_RT(double* hrt);
        void           getNetProductionRates(double* rr);
        double         molecularWeight(int k);
        int            speciesIndex(string name);
        int            nSpecies();
        int            nElements();
        double         meanMolecularWeight();
        void           setMoleFractions(double *x);
        void           setMassFractions(double *y);
        void           getMassFractions(double *y);
        void           getMoleFractions(double *x);
        void           setState_TPY( double T, double p, double *y);
        double         enthalpy_mass();
        int            elementIndex(string name);
        void           getAtoms(int k, double *Z);
        double         atomicWeight(int k);
        double         moleFraction(int k);
        double         massFraction(int k);
        double         pressure();
        void           setMassFractions_NoNorm(double *y);
        vector<string> speciesNames();
        string         speciesName(int x);

        IdealGasMix(string mechfile);

};

//-------------------------------------------------------------------------------


/** Cantera shell functions when not compiling with Cantera.  
 *  When cantera is not needed, compile with these
 *  routines that have the same headers as those of Cantera.  Add new routines
 *  when new cantera functions are called.  This also provides a template for
 *  implementing user-defined functions instead of calling Cantera.  However,
 *  the routines are not meant to be called, and return errors if they are.
 *  
 *  @author David O. Lignell
 */

class Transport {

    public :

        double viscosity();
        double thermalConductivity();
        void getMixDiffCoeffs(double *D);

};

//-------------------------------------------------------------------------------

Transport* newTransportMgr(string s1, IdealGasMix *gas);

} // end namespace Cantera

#endif
#endif

