/**
 * @file radiation.cc
 * Source file for class radiation
 */

#include "radiation.h"
#include "processor.h"
#include <cstdlib>
#include <cmath>
#include "diffuser.h"

using namespace std;
#ifdef CANTERA18
using Cantera_CXX::IdealGasMix;
#else
using Cantera::IdealGasMix;
#endif

extern processor proc;

///////////////////////////////////////////////////////////////////////////////
/** Constructor 
 *  @param odtpp \input odtParam object to set things up.
 *  @param gas \input pointer to the gas object
 */

radiation::radiation(odtParam &odtpp, IdealGasMix *gas) {

    if(!odtpp.Lrxn && !odtpp.ItableLookup) {
        cout << "\n#********** WARNING radiation is on, but Lrxn and ItableLookup are false" << endl;
        return;
    }

    nRadSp = 4;                // CH4, CO2, H2O, CO

    radCoefs.resize(nRadSp,vector<double>(6,0));  // kp (=) 1/atm*m

    radCoefs[0][0] =  1.017015E+1;    // ch4; kp=2.798 at 1150K
    radCoefs[0][1] = -7.947312E-03;
    radCoefs[0][2] =  4.342446E-7;
    radCoefs[0][3] =  1.048611E-9;
    radCoefs[0][4] = -2.287861E-13;
    radCoefs[0][5] =  0.000000E+0;
    radCoefs[1][0] =  3.24442E+1;     // co2; kp=29.197 at 925K
    radCoefs[1][1] =  7.537513E-02;
    radCoefs[1][2] = -1.535140E-04;
    radCoefs[1][3] =  9.48794E-8;
    radCoefs[1][4] = -2.509259E-11;
    radCoefs[1][5] =  2.447995E-15;
    radCoefs[2][0] =  6.86948E+1;     // h2o; kp=4.474  at 1119K
    radCoefs[2][1] = -1.52349E-01;
    radCoefs[2][2] =  1.417848E-04;
    radCoefs[2][3] = -6.620996E-8;
    radCoefs[2][4] =  1.52415E-11;
    radCoefs[2][5] = -1.373456E-15;
    radCoefs[3][0] =  1.56536E+0;    // co; kp=2.501 at 1007 K
    radCoefs[3][1] =  1.483914E-02;
    radCoefs[3][2] = -2.656035E-05;
    radCoefs[3][3] =  1.68798E-8;
    radCoefs[3][4] = -4.674473E-12;
    radCoefs[3][5] =  4.767887E-16;

    Imode = odtpp.Iradiation; 

    sigmaSB = 5.670E-8;       // W/m2*K4

    sootFactor = 1863.0;     // 1220

    TloBC  =  odtpp.inletTemp;   
    ThiBC  =  odtpp.freeStreamTemp;   

   if(Imode==1 && TloBC != ThiBC) {
       cout << "\n#********** WARNING in opthinRad TloBC != ThiBC, using ThiBC" << endl;
   }

   //-------------- set iRadIndx: rad species not in the mechanism are -1


   if(odtpp.Lrxn) {

       bool fmissing = false;
       iRadIndx.resize(nRadSp);
       int isp;

       isp = gas->speciesIndex("CH4");      
       isp = (isp > 0) ? isp : gas->speciesIndex("ch4");
       iRadIndx[0] = isp;
       if(isp < 0) fmissing = true;

       isp = gas->speciesIndex("CO2");      
       isp = (isp > 0) ? isp : gas->speciesIndex("co2");
       iRadIndx[1] = isp;
       if(isp < 0) fmissing = true;

       isp = gas->speciesIndex("H2O");      
       isp = (isp > 0) ? isp : gas->speciesIndex("h2o");
       iRadIndx[2] = isp;
       if(isp < 0) fmissing = true;

       isp = gas->speciesIndex("CO");      
       isp = (isp > 0) ? isp : gas->speciesIndex("co");
       iRadIndx[3] = isp;
       if(isp < 0) fmissing = true;

       if(fmissing) 
           cout << endl << "Warning one or more radiating species missing from mechanism" << endl;
   }


}

///////////////////////////////////////////////////////////////////////////////
/** Function computes radiatative heat source (W/m3).
 *  Use optically thin for Imode == 1, twoflux for Imod == 2.
 *  @param xMoleSp   \input xMoleSp[grid point][species] -  Species are all species in mechanism
 *  @param temp      \input Input temperatures (K)
 *  @param pressure  \input (Pa)
 *  @param xPosf     \input cell face positions
 *  @param radSource_G \output gas radiation volumetric heat source (W/m3)
 *  @param fvSoot    \input optional soot volume fraction
 */
void radiation::getRadHeatSource(const vector<vector<double> > &xMoleSp, 
                      const vector<double> &temp, 
                      const double pressure,
                      const vector<double> &xPosf,
                      vector<double> &radSource_G, 
                      vector<double> &radSource_P, 
                      vector<double> &ka,
                      particles *part, 
                      const vector<double> &fvSoot) {

    if(Imode == 1) {
        if(fvSoot.size() > 0) 
            opthinRadHeatSource(xMoleSp, temp, pressure, radSource_G, ka, fvSoot);
        else
            opthinRadHeatSource(xMoleSp, temp, pressure, radSource_G, ka);
    }
    else if(Imode == 2) {
        if(fvSoot.size() > 0) 
            twoFluxRadHeatSource(xMoleSp, temp, pressure, xPosf, radSource_G, radSource_P, ka, part, fvSoot);
        else {
            twoFluxRadHeatSource(xMoleSp, temp, pressure, xPosf, radSource_G, radSource_P, ka, part);
        }
    }

}

///////////////////////////////////////////////////////////////////////////////
/** Function computes optically thin volumetric radiative heat source.
 *  @param xMoleSp   \input xMole[grid point][species] - Species are all species in mechanism
 *  @param temp      \input Input temperature (K)
 *  @param pressure  \input (Pa)
 *  @param radSource_G \output (W/m3)
 *  @param fvSoot    \input optional soot volume fraction
 */

void radiation::opthinRadHeatSource(const vector<vector<double> > &xMoleSp, 
                                    const vector<double> &temp, 
                                    const double pressure,
                                    vector<double> &radSource_G,
                                    vector<double> &ka,
                                    const vector<double> &fvSoot) {
   if(Imode != 1) {
       cout << endl << "Error opthinRadHeatSource called but Imode is not right" << endl;
       exit(0);
   }

   int npts = radSource_G.size();  
   double Kabsorption;

   bool Bhave_ka = (ka.size()) ? 1 : 0;

   for(int i=0; i<npts; i++) {

       if(Bhave_ka) 
           Kabsorption = ka[i];
       else {
           if(fvSoot.size() == 0)
               Kabsorption = getGasAbsorptionCoefficient( xMoleSp[i], temp[i], pressure );
           else
               Kabsorption = getGasAbsorptionCoefficient( xMoleSp[i], temp[i], pressure, fvSoot[i] );
       }

       radSource_G[i] = -4 * sigmaSB * Kabsorption * (pow(temp[i],4.0)-pow(ThiBC,4.0));


   }
}

///////////////////////////////////////////////////////////////////////////////
/** Function computes mean absorption coefficient of gas
 *  \cond
 *  Kp = sum x_k * P * K_k
 *  K_k = a5*T^5 + a4*T^4 + a3*T^3 + a2*T^2 + a1*T + a0
 *  \endcond
 *
 *  \f[
 *  K_p = \sum_k{x_k * P * K_k}
 *  \f]
 *  \f[
 *  K_k = a_5 T^5 + a_4 T^4 + a_3 T^3 + a_2 T^2 + a_1 T + a_0
 *  \f]
 *
 *  @param xMole \input xMole[species] - Species are all species in mechanism
 *  @param T \input Input temperature (K)
 *  @param pressure \input (Pa)
 *  @param fvSoot \input optional soot volume fraction
 *  @return Plank mean absorption coeficient (1/m)
 */

double radiation::getGasAbsorptionCoefficient(const vector<double> &xMole, const double &T,
                                const double &pressure,
                                const double &fvSoot) {

    double P = pressure/101325.;         // Pa --> atm

    double Kabs;
    double Ktot = 0.0;

    for(int k=0; k<nRadSp; k++) {
        if(iRadIndx[k] < 0)          // this radiation species k is not in the mechanism
            continue; 
        Kabs = radCoefs[k][5];
        for(int j=4; j>=0; j--)
            Kabs = Kabs * T + radCoefs[k][j];

        Ktot += xMole[iRadIndx[k]]*P*Kabs;
    }

    if(fvSoot > 0.0)
        Ktot += sootFactor * fvSoot * T;

    return Ktot;
}

///////////////////////////////////////////////////////////////////////////////
/** Function computes radiative source (W/m3) at grid points using the two-flux model
*
 *
 *
 *  \cond
 *  Solving (13.30a) and (13.30b) in Modest 1993 Radiative Heat Transfer P. 492:
 *  dIp/dx =  2*k*Ib - 2*k*Ip,    bc x=0: Ip = sig*Tbclo^4 / pi
 *  dIm/dx = -2*k*Kb + 2*k*Im,    bc x=0: Im = sig*Tbchi^4 / pi
 *
 *  with del dot q = rad source = k*(4*pi*Ib - G) (eqn 8.54) = 2*pik*(2*Ib-Ip-Im) using (eqn 13.32)
 *
 *  Here, written in terms of q+ and q- not I+ and I- (Ip, Im).
 *  With particles: 
 *     dq+/ds =  2kg*sigma*Tg^4 + 2*sigma*sum_j(kp_j * Tp_j^4) - 2kg*q+ - 2*q+*sum_j(kp_j)
 *     dq-/ds = -2kg*sigma*Tg^4 - 2*sigma*sum_j(kp_j * Tp_j^4) + 2kg*q- + 2*q-*sum_j(kp_j)
 *  or: qout = qin + gas emission + particle emission - gas absorption - particle absorption
 * 
 *  The solution grid is dx between cell centers:
 *  odt grid:      | * |        *        |                       *                     |
 *  solution grid: * *          *                                *                     *
 *  with dx between stars in the solution grid
 *
 *  Verified against Example 13.4.  Also computed del dot q directly using (13.33) (verif).
 *  This was done for constant dx, and constant k, and uniform T.
 *  The two ODEs are integrated at odtline cell center points using a finite difference grid
 *    using implicit euler.
 *  That is: q+_i = [q+_{i-1} + ds*(2Kg*sig*Tg^4 + 2*sig*sum_j(Kp_j*Tp_j^4)) ]_i / (1+ds*(2Kg 2sum_j(Kp_j)))_i
 *  The terms in this equation are computed and stored in arrays.
 *  First compute q+, q-, then Gas source (W/m3) is dq+/ds - dq-/ds.
 *  The particle source W/m3 / #/m3 --> W/particle.  
 *  
 *  Particles in a given cell are assumed to live at the grid center for purposes of computing radiative sources
 *  \endcond
 *
 *  @param xMoleSp \input xMole[grid point][species];  Species are all species in mechanism
 *  @param temp \input Input temperature (K)
 *  @param pressure \input (Pa)
 *  @param xPosf \input cell face positions
 *  @param radSource \output (W/m3)
 *  @param fvSoot \input optional soot volume fraction
 */

void radiation::twoFluxRadHeatSource(const vector<vector<double> > &xMoleSp, 
                                     const vector<double> &temp, 
                                     const double pressure,
                                     const vector<double> &xPosf,
                                     vector<double> &radSource_G, 
                                     vector<double> &radSource_P, 
                                     vector<double> &ka, 
                                     particles *part,
                                     const vector<double> &fvSoot) {

   if(Imode != 2) {
       cout << endl << "Error twoFluxRad called but Imode is not right" << endl;
       exit(0);
   }

   //-------------

   int npts = radSource_G.size()+2;       // add in the two boundaries  
   vector<double> qp(npts);
   vector<double> qm(npts);
   vector<double> Kabs(npts);
   vector<double> dx(npts-1);
   vector<double> gasEmmTerm(npts);       // 2*kg*sigma*Tg^4

   //------------- Get the dx array

   dx[0] = 0.5*(xPosf[1]-xPosf[0]);
   for(int i=1; i<npts-2; i++)
       dx[i] = 0.5*(xPosf[i+1]-xPosf[i-1]);
   dx[npts-2] = 0.5*(xPosf[xPosf.size()-1]-xPosf[xPosf.size()-2]);

   //------------- set particle number densities in the dx array

   vector<double> nDens(0,0.0);
   if(part) {
       part->set_iyPos();         // may already be set --> redundant ?
       nDens.resize(part->nPart,0.0);
       for(int k=0; k<part->nPart; k++) {
           if(!part->pActive[k])
                   continue;
           int idxPart = (part->yPos[k] <= 0.5*(xPosf[part->iyPos[k]]+xPosf[part->iyPos[k]+1])) 
                        ? part->iyPos[k] : part->iyPos[k] + 1;
           nDens[k] = part->nInPseudoPart[k] / dx[idxPart];
       }
   }

   //------------- Get the gas absorption coefficient

   bool Lhave_ka   = (ka.size() > 0)     ? true : false;
   bool Lhave_soot = (fvSoot.size() > 0) ? true : false;

   for(int i=0; i<npts-2; i++) {

       if(Lhave_ka) {
           Kabs[i+1] = ka[i];
           if(Lhave_soot){
               Kabs[i+1] += sootFactor * fvSoot[i] * temp[i];
           }
           
       }
       else {
           if(Lhave_soot)
               Kabs[i+1] = getGasAbsorptionCoefficient( xMoleSp[i], temp[i], pressure, fvSoot[i] );
           else
               Kabs[i+1] = getGasAbsorptionCoefficient( xMoleSp[i], temp[i], pressure );
       }
   }

   Kabs[0] = Kabs[1];
   Kabs[npts-1] = Kabs[npts-2];

   //------------- Get the gas emmision term 

   for(int i=0; i<npts-2; i++) 
       gasEmmTerm[i+1] = 2.0*Kabs[i+1] * sigmaSB*pow(temp[i],4.0); 
   gasEmmTerm[0]      = gasEmmTerm[1];
   gasEmmTerm[npts-1] = gasEmmTerm[npts-2];

   //------------- for particles

   vector<double> partEmmTerm;        // 2*sigma*sum_j(kp_j*Tp_j^4) in each cell
   vector<double> Kabs_p;             // sum_j(kp_j) in each cell

   if(part) {
       partEmmTerm.resize(npts,0.0);    
       Kabs_p.resize(npts,0.0);        

       for(int k=0; k<part->nPart; k++) {
           if(!part->pActive[k]) continue;
           double kk; 
           if(part->pShape==1)
               kk = part->emiss*2.0/M_PI*(2*part->pRadi[k]*part->pLength)*nDens[k];       // absorption coefficient of particle k
           else if(part->pShape==2)
               kk = part->emiss*M_PI*part->pRadi[k]*part->pRadi[k]*nDens[k];       // absorption coefficient of particle k
           else {
               cout << "\nERROR: in radiation only pShape=1 or 2 allowed" << endl;
               exit(0);
           }
           Kabs_p[ part->iyPos[k]+1 ]      += kk;
           double partSurfTemp = part->getSurfTemp(k);
           partEmmTerm[ part->iyPos[k]+1 ] += 2.0*sigmaSB*kk*pow(partSurfTemp,4.0);
       }
       Kabs_p[0]           = Kabs_p[1];
       Kabs_p[npts-1]      = Kabs_p[npts-2];
       partEmmTerm[0]      = partEmmTerm[1];
       partEmmTerm[npts-1] = partEmmTerm[npts-2];
   }

   //------------- Get radiative fluxes: qp, qm
   // Marching qp low to high, qm high to low

   double qpBClo = sigmaSB * pow(TloBC,4.0);
   double qmBChi = sigmaSB * pow(ThiBC,4.0);

   qp[0] = qpBClo;
       for(int i=1; i<npts; i++)
       qp[i] = ( qp[i-1] + dx[i-1]*gasEmmTerm[i] + ((!part) ? 0 : dx[i-1]*partEmmTerm[i]) )  /
               ( 1. + 2.*Kabs[i]*dx[i-1]         + ((!part) ? 0 : 2.0*dx[i-1]*Kabs_p[i]) );

   qm[npts-1] = qmBChi;
       for(int i=npts-2; i>=0; i--)
           qm[i] = ( qm[i+1] + dx[i]*gasEmmTerm[i] + ((!part) ? 0 : dx[i]*partEmmTerm[i]) )  /
                   ( 1. + 2.*Kabs[i]*dx[i]         + ((!part) ? 0 : 2.0*dx[i]*Kabs_p[i]) );

   //------------- Get gas radiative source: (div q) (=) W/m3 (additive source to E-bal)

    for(int i=0, ip=1; i<npts-2; i++, ip++){
       radSource_G[i] = -2.*Kabs[ip] * (2.*sigmaSB*pow(temp[i],4.0) - qp[ip] - qm[ip]);
    }

   //------------- Get particle radiative sources: J/s (additive) per individual particle (not collection)

   if(part) {
       for(int k=0; k<part->nPart; k++) {
           if(!part->pActive[k]) continue;
           double kk; 
           if(part->pShape==1)
               kk = part->emiss*2.0/M_PI*(2*part->pRadi[k]*part->pLength)*nDens[k];       // absorption coefficient of particle k
           else if(part->pShape==2)
               kk = part->emiss*M_PI*part->pRadi[k]*part->pRadi[k]*nDens[k];       // absorption coefficient of particle k
           else {
               cout << "\nERROR: in radiation only pShape=1 or 2 allowed" << endl;
               exit(0);
           }
           int    i  = part->iyPos[k];
           double partSurfTemp = part->getSurfTemp(k);
           radSource_P[k] = -2.*kk * (2.*sigmaSB*pow(partSurfTemp,4.0) -qp[i]-qm[i]) / nDens[k]; // W/particle
       }
   }

}

////////////////////////////////DOXYGEN DOCUMENTATION//////////////////////////////////

/*! \fn void radiation::twoFluxRadHeatSource(const vector<vector<double> > &xMoleSp, 
                                     const vector<double> &temp, 
                                     const double pressure,
                                     const vector<double> &xPosf,
                                     vector<double> &radSource_G, 
                                     vector<double> &radSource_P, 
                                     vector<double> &ka, 
                                     particles *part,
                                     const vector<double> &fvSoot)
 *
 *  Solving (13.30a) and (13.30b) in Modest 1993 Radiative Heat Transfer P. 492:
 *  \f[
 *      \frac{d I_p}{dx} = 2 k I_b - 2 k I_p,    bc: x=0: Ip = \sigma \frac{{Tbc_{lo}}^4 }{ \pi}
 *  \f]
 *  \f[
 *      \frac{d I_m}{dx} = -2 k K_b + 2 k I_m,    bc: x=0: Im = \sigma \frac{{Tbc_{hi}}^4}{ \pi}
 *  \f]
 *
 *  with
 *  \f[\nabla \cdot q = \text{rad source} = k (4 \pi I_b - G) \f]
 *  \f[\text{(eqn 8.54)} = 2 \pi k (2 I_b-I_p-I_m) \text{ using (eqn 13.32)} \f]
 *
 *  Here, written in terms of \fun{q+} and \fun{q-} not \fun{I+} and \fun{I-} (\fun{Ip}, \fun{Im}).
 *  With particles: 
 *     \f[
 *         \frac{dq+}{ds} =  2kg * \sigma * {Tg}^4 + 2*\sigma*\sum_j{kp_j * {{Tp}_j}^4} - 2kg*{q+} - 2*{q+}*\sum_j{kp_j}
 *     \f]
 *     \f[
 *         \frac{d{q-}}{ds} = -2kg*\sigma*{Tg}^4 - 2*\sigma*\sum_j{kp_j * {{Tp}_j}^4} + 2kg*{q-} + 2*{q-}*\sum_j{kp_j}
 *		 \f]
 *  or: 
 *     \f[
 *          q_{out} = q_{in} + \text{ gas emission } + \text{ particle emission } - \text{ gas absorption } - \text{ particle absorption }
 *     \f]
 * 
 *  The solution grid is dx between cell centers:\vc{
 *  odt grid:      | * |        *        |                       *                     |
 *  solution grid: * *          *                                *                     *
 *  }
 *  with dx between stars in the solution grid
 *
 *  Verified against Example 13.4.  Also computed \fun{\nabla \cdot q} directly using (13.33) (verif).
 *  This was done for constant \fun{dx}, constant \fun{k}, and uniform \fun{T}.
 *  The two ODEs are integrated at odtline cell center points using a finite difference grid
 *    using implicit euler.
 *  That is: 
 *  \f[
 *      {q+}_i = \frac{[{q+}_{i-1} + ds*(2Kg*\sigma*{Tg}^4 + 2*\sigma*\sum_j{Kp_j*{{Tp}_j}^4}) ]_i }{ [1+ds*(2Kg 2\sum_j{Kp_j})]_i}
 *  \f]
 *  The terms in this equation are computed and stored in arrays.
 *  First compute \fun{q+}, \fun{q-}, then Gas source (\fun{W/m^3}) is \fun{dq+/ds - dq-/ds}.
 *  The particle source \fun{(W/m^3) / (\#/m^3) \to W/\text{particle}}.  
 *  
 *  Particles in a given cell are assumed to live at the grid center for purposes of computing radiative sources.
 */
