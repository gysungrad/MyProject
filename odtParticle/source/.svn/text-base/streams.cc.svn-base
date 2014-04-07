/**
 * @file streams.cc
 * Header file for class streams
 */

#include "streams.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "processor.h"

using namespace std;
#ifdef CANTERA18
using Cantera_CXX::IdealGasMix;
#else
using Cantera::IdealGasMix;
#endif

using namespace std;

extern processor proc;

///////////////////////////////////////////////////////////////////////////////

/** Constructor.
 *
 * @param cantIG \input ideal gas object to set pointer gas
 */
streams::streams(IdealGasMix *cantIG) {


#ifdef DOCANTERA

    gas = cantIG;

    nspc = gas->nSpecies();
    y0   = vector<double>(nspc, 0.0);
    y1   = vector<double>(nspc, 0.0);
    T0   = 300.0;
    T1   = 300.0;

    y0Labels = vector<string>(nspc);
    y1Labels = vector<string>(nspc);

    pres = 101325.0;                      // overwritten in readStreams below

    
    readStreams();
    
    if(comp2==0) setStoicMixf(); //fix this should be 7 which doesn't exist anymore

    getMixtureFraction(&y0[0], true);    // set beta0 and beta1

#endif   // DOCANTERA


}
///////////////////////////////////////////////////////////////////////////////

/** Read streams file.
 *
 *  @param fname \input streams input file name.
 */
void streams::readStreams(string fname) {

    *proc.ostrm << endl << "# ************* Reading streams.inp input file" << endl;

    streamInput_.setFile(fname);

    // first: Lmole?
    bool Lmole = false;

    string tString;
    streamInput_.getParameter("moleOrMass", &tString);

    if (tString=="MOLE" 
     || tString=="mole"
     || tString=="Mole") {
      Lmole = true;
    }
    
    // next is the pressure
    streamInput_.getParameter ("pres", &pres);

    /// Questions:
    ///
    /// This is the streams class, not the stream class
    /// But why fix the number of streams at 2?
    /// 
    /// What about more complicated mixing?
    /// (Like, a primary inlet, a secondary inlet...
    ///  two gas mixture fractions, three "streams")
    /// - coal gas mf
    /// - primary/secondary mf
    
    /// Need to provide some way of grouping multiple stream variables
    /// e.g. T and composition vector
    ///
    /// Once you have that you're good
    ///
    /// For now, just assume 2 streams, fuel and oxidizer
    ///
    /// Maybe need a "subclass" for Stream type
   
#ifdef DOCANTERA
    /// Find T and Y
    /// USES CANTERA
    /// 1. grab values
    /// 2. check to make sure they are in mechanism file
    /// 3. normalize mass fractions, 0 to 1
    
    streamInput_.getParameter("composition2", &comp2, 0);
    if(comp2==0)  {  
      streamInput_.getParameter("T0", &T0);
      streamInput_.getVector("composition0", &y0, &y0Labels);
      setStreamComp(y0,y0Labels);
  

      // Find second stream's T and Y
      // USES CANTERA
      streamInput_.getParameter("T1", &T1);
      streamInput_.getVector("composition1", &y1, &y1Labels);     
      
      setStreamComp(y1,y1Labels);
      // convert from mole -> mass if necessary
      // USES CANTERA
      gas->setMoleFractions( &y0[0] );
      gas->getMassFractions( &y0[0] );
      gas->setMoleFractions( &y1[0] );
      gas->getMassFractions( &y1[0] );
      // compute H0, H1 (independent of pressure)
      // USES CANTERA
      gas->setState_TPY( T0, pres, &y0[0] );
      h0 = gas->enthalpy_mass();
      gas->setState_TPY( T1, pres, &y1[0] );
      h1 = gas->enthalpy_mass();
      
      *proc.ostrm << endl << "# ************* Stream_0 properties " << endl;
      *proc.ostrm << endl << "T0: "<<T0 << endl;
      for(int i=0; i < y0.size(); i++){	      *proc.ostrm << "Species index: "<< gas->speciesIndex(y0Labels[i])<<"	" <<y0Labels[i]<<"	" <<y0[i]<< endl;	}
      
      *proc.ostrm << endl << "# ************* Stream_1 properties " << endl;
      *proc.ostrm << endl << "T1: "<<T1 << endl;
      for(int i=0; i < y1.size(); i++){	      *proc.ostrm << "Species index: "<< gas->speciesIndex(y1Labels[i])<<"	" <<y1Labels[i]<<"	" <<y1[i]<< endl;	}

   }
      
  #else
  // This section is not called since this function is called inside a DOCANTERA block in the streams constructor.
  // Fix this or remove it
  // --dol
  // With your change putting in the curly bracket this part is as it was (we didn't change this). Do what you please with this. Alan and Christina
      *proc.ostrm  << endl;
      *proc.ostrm  << "WARNING: streams:readStreams(): a routine is required to access "
	    << "the enthalpy, temperature, mass, and mole fractions of the gas, "
	    << "not to mention a bunch of other stuff to come later on." 
	    << endl;
      *proc.ostrm  << endl;
      *proc.ostrm  << "You should define these routines for your particular case elsewhere "
	    << "(for example, cantera_shell_functions.cc/.h or my_scalar_type.cc/.h)."
	    << endl;
      *proc.ostrm  << endl;
  #endif
}

///////////////////////////////////////////////////////////////////////////////

/** For a given simulation, if a Cantera input file is given, the species for that simulation
  * are specified in the Cantera input file.  This means that the mole fractions or
  * other composition data must be put in the "right" order - that is, an order that matches
  * Cantera's list of species.
  *
  * This also means that species not listed in the input file must have their mole fractions set
  * to 0.
  *
  * @param y \input vector of mole fractions.
  * @param yLabels \input  vector of species labels (same size and same order as y).
  */
void streams::setStreamComp(vector<double> &y, vector<string> &yLabels) {

  int    z;             // counter
  double sum_molefracs; // cumulative sum of mole fractions

  //---------- create a temporary vector to hold the (correct) information

  vector<double> y_temp = vector<double>(nspc, 0.0);

  //---------- for each species (that is, for each species label), 
  //           find the corresponding cantera index

  z = 0;
  for (vector<string>::iterator iLabels = yLabels.begin(); iLabels != yLabels.end(); ++iLabels) {
    int isp = gas->speciesIndex(*iLabels);
    if (isp == -1) {
      // THROW EXCEPTION: SPECIES SPECIFIED IN INPUT FILE COULD NOT BE FOUND IN MECHANSIM FILE
    }
    y_temp[isp] = y[z];
    ++z;
  }

  //---------- now normalize if the species fractions don't add up to 1

  sum_molefracs = 0;
  for (vector<double>::iterator iY = y_temp.begin(); iY != y_temp.end(); ++iY) {
    sum_molefracs += (*iY);        // comulative sum of all mole fractions
  }
  if (sum_molefracs == 0.0) {
    // THROW EXCPTION: NULL COMPOSITION VECTOR!!! SPECIFY COMPOSITIONS IN STREAMS.INP (or elsewhere)
  } else if (sum_molefracs != 1.0) {
    // Normalize
    for (vector<double>::iterator iY = y_temp.begin(); iY != y_temp.end(); ++iY) {
      (*iY) = (*iY)/sum_molefracs;
    }
  }

  //---------- Put the correct y vector (in y_temp) into the original y vector (y)
  //           (Note: species not specified in the input file are initialized to 0.0,
  //           so they don't need to be addressed here explicitly)

  z = 0;
  y.clear();
  yLabels.clear();							//ZJ
  
  for (vector<double>::iterator iY = y_temp.begin(); iY != y_temp.end(); ++iY) {
    y.push_back(*iY);
    yLabels.push_back(gas->speciesName(z));				//ZJ
    ++z;
  }

}


///////////////////////////////////////////////////////////////////////////////

/** Read stream composition from file. 
 *
 * @param ifile \input file name to read.
 * @param xy \input species concentrations.
 */
void streams::readStreamComp(ifstream &ifile, vector<double> &xy){

  *proc.ostrm << "******* WARNING: streams::readStreamComp(): obsolete method: use streams::setStreamComp() with inputFile class." << endl;

   string s1;
   double d1;
   int    isp;
   
   for(;;) {
       ifile >> s1;
       if(s1 == "END" || s1=="end" || s1=="End" || ifile.eof())
           break;
       ifile >> d1;
       isp = gas->speciesIndex(s1);
       if(isp==-1) {
           *proc.ostrm << endl << "ERROR in readStreamComp: species " << s1 
                        << " not found (compare to mechanism file)" << endl;
           exit(0);
       }
       xy[isp] = d1;
   }
   
   d1 = 0.0;
   for(int i=0; i<(int)xy.size(); i++)
       d1 += xy[i];
   if(d1 == 0.0) {
       *proc.ostrm << endl << "ERROR in readStreamComp: empty composition vector" << endl;
       exit(0);
   }
   if(d1 != 1.0) 
       for(int i=0; i<(int)xy.size(); i++)
           xy[i] /= d1;

}

///////////////////////////////////////////////////////////////////////////////

/** Computes the temperature and composition of complete combustion at the given mixf.
 *  For nonpremixed flames (don't do anything funny, like have oxygen in the fuel stream)
 *
 *  @param mixf \input mixture fraction, defines elemental composition.
 *  @param ypcc \output mass fractions of products of complete combustion.
 *  @param hpcc \output enthalpy of products of complete combustion.
 *  @param Tpcc \output temperature of products of complete combustion.
 */
void streams::getProdOfCompleteComb(double mixf, vector<double> &ypcc, 
                                    double &hpcc, double &Tpcc, int probType) {

  if(probType==6) return;	
   //---------- Compute the mixing mass fractions and enthalpy
   
   ypcc.resize(nspc);
   double d1 = 1.0-mixf;
   for(int k=0; k<nspc; k++) {
       ypcc[k] = y1[k]*mixf + y0[k]*d1;
   }
   hpcc = h1*mixf + h0*d1;
   
  if(probType!=4)					// for opposedJets, don't start with combustion
  {
   //--------- Set gas and element indicicies

   int iC = gas->elementIndex("C");
   int iH = gas->elementIndex("H");
   //int iO = gas->elementIndex("O"); // !!!!!  currently unusedvariable
   int iN = gas->elementIndex("N");
   int iCO2 = gas->speciesIndex("CO2");
   int iH2O = gas->speciesIndex("H2O");
   int iN2  = gas->speciesIndex("N2");
   int iO2  = gas->speciesIndex("O2");

   //---------- Set ypcc as the mixing mole fractions: Take a basis of one mole:
   // now we are working in moles
   // elemM are moles of each element
   // when stoic: 
   // CxHyNz   + (x+y/4)O2  ==>  (z/2)N2 + (x)CO2 + (y/2)H2O
   // otherwise:
   // CxHyNz   + (beta)O2   ==>  (z/2)N2 + (x)CO2 + (y/2)H2O
   // Note this allows fuels with nitrogen and oxygen

   gas->setMassFractions( &ypcc[0] );
   gas->getMoleFractions( &ypcc[0] );

   double nOnotFromO2  = 0.0;
   double nHnotFromH2O = 0.0;
   double nCnotFromCO2 = 0.0;
   vector<double> elemM = getElementMoles( &ypcc[0], nOnotFromO2, 
                                                     nHnotFromH2O, nCnotFromCO2 );          

   double x    = elemM[iC];
   double y    = elemM[iH];
   double z    = elemM[iN];
   double beta = elemM[gas->elementIndex("O")] * 0.5;        // moles of O as O2

   //--------------------------------------------------------------------------

   if(mixf < mixfStoic) {                        // lean: burn all fuel, leftover O2

       for(int k=0; k<nspc; k++)
           ypcc[k] = 0.0;
       
       if(iCO2 > 0)                              // CO2 is not in the H2 mechs
           ypcc[iCO2] = x;         
       ypcc[iH2O] = y*0.5;
       ypcc[iN2]  = z*0.5;
       ypcc[iO2]  = beta - (x+y/4.0);

   }
   else{                                         // rich: burn all O2, leftover fuel

       //double eta = beta/(x+y/4.0); // extent of reaction
       double eta = (beta-nOnotFromO2/2)/(x+y/4-nOnotFromO2/2); // extent of reaction
       if(eta > 1.0)
           *proc.ostrm << endl << "eta > 1.0" << endl;
       d1 = 1.0-eta;                            // fraction of fuel unburnt

       for(int k=0; k<nspc; k++)       
           ypcc[k] *= d1;                       // initialize everything then correct
       if(iCO2 > 0)                             // CO2 is not in the H2 mechs
           ypcc[iCO2] = (x-nCnotFromCO2)     + nCnotFromCO2*eta;       // init + formed
       ypcc[iH2O] = (y-nHnotFromH2O)*0.5 + nHnotFromH2O*0.5*eta;   // init + formed
       ypcc[iN2]  = z*0.5;
       ypcc[iO2]  = 0.0;

   }

   //--------------------------------------------------------------------------

   double sum = 0.0;                       // normalize moles
   for(int k=0; k<nspc; k++)
       sum += ypcc[k];
   assert(sum != 0.0);
   for(int k=0; k<nspc; k++)
       ypcc[k] /= sum;
   gas->setMoleFractions( &ypcc[0] );      // set mole fractions
   gas->getMassFractions( &ypcc[0] );      // set ypcc as mass fractions


   //--------------------------------------------------------------------------

   gas->setState_HP(hpcc, pres, 1.E-10);    // get temperature as Tadiabatic
   //gas->setState_HP(hpcc, pres);    // get temperature as Tadiabatic
   Tpcc = gas->temperature();
  }
  
  
}

///////////////////////////////////////////////////////////////////////////////

/** Set the stoichiometric mixture fraction using Bilger's definition */

void streams::setStoicMixf() {

    vector<double> x_c(nspc,0.0);

    double mc0, mc1, mo0, mo1, mh0, mh1;

    vector<double> elemMassFrac0 = setElementMassFracs(&y0[0]);
    vector<double> elemMassFrac1 = setElementMassFracs(&y1[0]);

    mc0 = elemMassFrac0[gas->elementIndex("C")]/
          gas->atomicWeight(gas->elementIndex("C"));
    mc1 = elemMassFrac1[gas->elementIndex("C")]/
          gas->atomicWeight(gas->elementIndex("C"));
    mh0 = elemMassFrac0[gas->elementIndex("H")]/
          gas->atomicWeight(gas->elementIndex("H"));
    mh1 = elemMassFrac1[gas->elementIndex("H")]/
          gas->atomicWeight(gas->elementIndex("H"));
    mo0 = elemMassFrac0[gas->elementIndex("O")]/
          gas->atomicWeight(gas->elementIndex("O"));
    mo1 = elemMassFrac1[gas->elementIndex("O")]/
          gas->atomicWeight(gas->elementIndex("O"));

    mixfStoic = (2.0*mc0 + 0.5*mh0 - mo0) / 
                (mo1-mo0 + 2.0*(mc0-mc1) + 0.5*(mh0-mh1));

    *proc.ostrm << endl << "mixfStoic = m_fuel/(m_fuel+m_air) = " << mixfStoic << endl;

}

///////////////////////////////////////////////////////////////////////////////

/** Sets the elements to have the correct Mass Fractions based on the specified array.
 *  @param y \input mass fraction array to use to get corresponding element fractions.
 *  @return vector of element mass fractions.
 */

vector<double> streams::setElementMassFracs(double *y) {


    vector<double> atomArr(gas->nElements());
    vector<double> elemMassFrac(gas->nElements(), 0.0);
    double sum = 0.0;

    gas->setMassFractions( &y[0] );

    for(int k=0; k<nspc; k++) {
        sum=0.0;
        gas->getAtoms(k, &atomArr[0]);           // [nelements] in sp k
        for(int m=0; m<(int)atomArr.size(); m++)
            sum += atomArr[m]*gas->atomicWeight(m);
        for(int m=0; m<(int)atomArr.size(); m++) 
            elemMassFrac[m] += y[k] * atomArr[m]/sum*gas->atomicWeight(m);
                              // is * mass frac of elem in sp
    }

    return elemMassFrac;

}

///////////////////////////////////////////////////////////////////////////////

/** Sets the elements to have the correct Mole Fractions based on the specified array.
 *  @param y \input mass fraction array to use to get corresponding element fractions.
 *  @return vector of element mole fractions.
 */

vector<double> streams::setElementMoleFracs(double *y) {


    vector<double> atomArr(gas->nElements());
    vector<double> elemMoleFrac(gas->nElements(), 0.0);

    gas->setMassFractions( &y[0] );

    for(int k=0; k<nspc; k++) {
        gas->getAtoms(k, &atomArr[0]);           // [nelements] in sp k
        for(int m=0; m<(int)atomArr.size(); m++) 
            elemMoleFrac[m] += gas->moleFraction(k) * atomArr[m];
    }
    double sum = 0.0;
    for(int m=0; m<(int)atomArr.size(); m++)
        sum += elemMoleFrac[m];
    assert(sum != 0.0);
    for(int m=0; m<(int)atomArr.size(); m++)
        elemMoleFrac[m] /= sum;

    return elemMoleFrac;

}

///////////////////////////////////////////////////////////////////////////////

/** Get amount of moles for each element.
 *  @param x \input pointer to vector of species mole fractions.
 *  @param nOnotFromO2 \input number of moles of oxygen not from O2 (oxygen in the base fuel).
 *  @param nHnotFromH2O \input number of moles of hydrogen not from H2O.
 *  @param nCnotFromCO2 \input number of moles of carbon not from CO2.
 *  @return vector of element moles.
 */

vector<double> streams::getElementMoles(double *x, double &nOnotFromO2,
                                                   double &nHnotFromH2O,
                                                   double &nCnotFromCO2) {


    vector<double> atomArr(gas->nElements());
    vector<double> elemM(gas->nElements(), 0.0);
    int iO2  = gas->speciesIndex("O2");
    int iO   = gas->elementIndex("O");
    int iCO2 = gas->speciesIndex("CO2");
    int iC   = gas->elementIndex("C");
    int iH2O = gas->speciesIndex("H2O");
    int iH   = gas->elementIndex("H");

    for(int k=0; k<nspc; k++) {
        gas->getAtoms(k, &atomArr[0]);           // [nelements] in sp k
        for(int m=0; m<(int)atomArr.size(); m++) 
            elemM[m] += x[k] * atomArr[m];
        if(k != iO2)  nOnotFromO2  += atomArr[iO] * x[k];
        if(k != iCO2) nCnotFromCO2 += atomArr[iC] * x[k];
        if(k != iH2O) nHnotFromH2O += atomArr[iH] * x[k];
    }
    return elemM;

}

///////////////////////////////////////////////////////////////////////////////

/**Compute the mixture fraction from the mass fractions using Bilger's mixf.
 * Set doBeta01=true on first call to initialize members beta0, beta1.
 * Later calls of this function only use the first parameter.
 *
 * @param y \input vector of species mass fractions.
 * @param doBeta01 \input flag=true on first call to set members beta0, beta1.
 */
double streams::getMixtureFraction(double *y, bool doBeta01) {


    double gO = -1.0/gas->atomicWeight(gas->elementIndex("O"));
    double gC =  2.0/gas->atomicWeight(gas->elementIndex("C"));
    double gH =  0.5/gas->atomicWeight(gas->elementIndex("H"));
    double gN =  0.0;        

    vector<double> elemMF;
    double         beta;

    elemMF = setElementMassFracs(y);
    beta   = gC * elemMF[gas->elementIndex("C")] + 
             gH * elemMF[gas->elementIndex("H")] + 
             gO * elemMF[gas->elementIndex("O")] + 
             gN * elemMF[gas->elementIndex("N")];

    if(doBeta01) {
        elemMF = setElementMassFracs(&y0[0]);
        beta0  = gC * elemMF[gas->elementIndex("C")] + 
            gH * elemMF[gas->elementIndex("H")] + 
            gO * elemMF[gas->elementIndex("O")] + 
            gN * elemMF[gas->elementIndex("N")];

        elemMF = setElementMassFracs(&y1[0]);
        beta1  = gC * elemMF[gas->elementIndex("C")] + 
             gH * elemMF[gas->elementIndex("H")] + 
             gO * elemMF[gas->elementIndex("O")] + 
             gN * elemMF[gas->elementIndex("N")];
    }

    return (beta - beta0)/(beta1 - beta0);

}

