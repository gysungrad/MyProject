/**
 * @file cantera_shell_functions.cc
 * Source file for class cantera_shell_functions
 */

#ifndef DOCANTERA

// This file is used to store dummy function definitions when cantera is not available
// Just for compilation

#include "cantera_shell_functions.h"
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace Cantera;

//-------------------------------------------------------------------------------

void errorMessage(string name) {
    cout << "\nERROR, FUNCTIONS NOT IMPLEMENTED IN CANTERA_SHELL_FUNCTIONS " 
         << name << endl;
    exit(0);
}
void warningMessage(string name) {
    cout << "\nWarning, FUNCTIONS NOT IMPLEMENTED IN CANTERA_SHELL_FUNCTIONS " 
         << name << endl;
}

//-------------------------------------------------------------------------------

void   IdealGasMix::setState_PY(double p, double* y)             { errorMessage("setState_PY"); }
void   IdealGasMix::setState_HP(double h, double p, double tol)  { errorMessage("setState_HP"); }
double IdealGasMix::temperature()                                { errorMessage("temperature"); return 0; }
double IdealGasMix::density()                                    { errorMessage("density"); return 0; }
double IdealGasMix::cp_mole()                                    { errorMessage("cp_mole"); return 0; }
double IdealGasMix::cp_mass()                                    { errorMessage("cp_mass"); return 0; }
void   IdealGasMix::getEnthalpy_RT(double* hrt)                  { errorMessage("getEnthalpy_RT"); }
void   IdealGasMix::getNetProductionRates(double* rr)            { errorMessage("getNetProductionRates"); }
double IdealGasMix::molecularWeight(int k)                       { errorMessage("molecularWeight"); return 0; }
int    IdealGasMix::speciesIndex(string name)                    { errorMessage("speciesIndex"); return 0; }
int    IdealGasMix::nSpecies()                                   { errorMessage("nSpecies"); return 0; }
int    IdealGasMix::nElements()                                  { errorMessage("nElements"); return 0; }
double IdealGasMix::meanMolecularWeight()                        { errorMessage("meanMolecularWeight"); return 0; }
void   IdealGasMix::setMoleFractions(double *x)                  { errorMessage("setMoleFractions"); }
void   IdealGasMix::setMassFractions(double *y)                  { errorMessage("setMassFractions"); }
void   IdealGasMix::getMassFractions(double *y)                  { errorMessage("getMassFractions"); }
void   IdealGasMix::getMoleFractions(double *x)                  { errorMessage("getMoleFractions"); }
void   IdealGasMix::setState_TPY( double T, double p, double *y) { errorMessage("setState_TPY"); }
double IdealGasMix::enthalpy_mass()                              { errorMessage("enthalpy_mass"); return 0; }
int    IdealGasMix::elementIndex(string name)                    { errorMessage("elementIndex"); return 0; }
void   IdealGasMix::getAtoms(int k, double *Z)                   { errorMessage("getAtoms"); }
double IdealGasMix::atomicWeight(int k)                          { errorMessage("atomicWeight"); return 0; }
double IdealGasMix::moleFraction(int k)                          { errorMessage("moleFraction"); return 0; }
double IdealGasMix::massFraction(int k)                          { errorMessage("massFraction"); return 0; }
double IdealGasMix::pressure()                                   { errorMessage("pressure"); return 0; }
void   IdealGasMix::setMassFractions_NoNorm(double *y)           { errorMessage("setMassFractions_NoNorm"); }
vector<string> IdealGasMix::speciesNames()                       { errorMessage("speciesNames"); return vector<string>(); }
string IdealGasMix::speciesName( int x)                          { errorMessage("speciesName"); return string(); }
IdealGasMix::IdealGasMix(string mechfile)                        { warningMessage("IdealGasMix"); }

//-------------------------------------------------------------------------------

double Transport::viscosity()                                    { errorMessage("viscosity"); return 0; }
double Transport::thermalConductivity()                          { errorMessage("thermalCond"); return 0; }
void   Transport::getMixDiffCoeffs(double *D)                    { errorMessage("getMixDiffC"); }

//-------------------------------------------------------------------------------

Transport* Cantera::newTransportMgr(string s1, IdealGasMix *gas)          { warningMessage("newTransportMgr"); return 0; }

#endif    //(*proc.ostrm).flush();


