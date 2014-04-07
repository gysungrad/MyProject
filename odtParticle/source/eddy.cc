/**
 * @file eddy.cc
 * Source file for class eddy
 */

#include "eddy.h"
#include "odtParam.h"
#include <cmath>   
#include <iomanip>   
#include <stdexcept>
#include <algorithm>

using namespace std;

///////////////////////////////////////////////////////////////////////////////

/** Default constructor function. */
eddy::eddy() { 

    leftEdge = 0.0;
    rightEdge = 1.0;
    invTauEddy = 1.0;
    cCoef = vector<double> (3,0.0); 
    bCoef = vector<double> (3,0.0); 
    partMomSrcInEddy = vector<double> (3,0.0);
    partErgSrcInEddy = vector<double> (3,0.0);
    LperiodicEddy = false;
    K = vector<double> (1,0.0); //AKeddy arbitrary initialization

}

///////////////////////////////////////////////////////////////////////////////

/** Copy constructor function  .
 *
 *  @param ed \input eddy object to copy.
 */
eddy::eddy(const eddy &ed) {

    leftEdge      = ed.leftEdge;
    rightEdge     = ed.rightEdge;
    eddySize      = ed.eddySize;
    Pa            = ed.Pa;
    invTauEddy    = ed.invTauEddy;
    cCoef         = ed.cCoef;
    bCoef         = ed.bCoef;
    partMomSrcInEddy = ed.partMomSrcInEddy;
    partErgSrcInEddy = ed.partMomSrcInEddy;
    LperiodicEddy = ed.LperiodicEddy;
    K = vector<double> (1,0.0); //AKeddy arbitrary initialization
    eddyLife      = ed.eddyLife;

}

///////////////////////////////////////////////////////////////////////////////

/** Sample an eddy size from the eddy size distribution.
 * This could be as simple as a tophat, but this function is
 * more accurate, hence more efficient.
 * 
 * @param odtP \input parameters object.
 * @param rr   \input random number generator object.
 */
void eddy::sampleEddySize(const odtParam &odtP) {

    double rndd = 0.0;
    eddySize = 0.0;
    
    // I have added this while-loop due to the problem described in function
    // eddy::sampleEddyPosition.
    do
    {
        rndd = rnd->getRand();
        if(!odtP.Llem) eddySize = odtP.esdp1 / log( rndd * odtP.esdp2 + odtP.esdp3 );
        else eddySize = odtP.Lmin * pow(( 1.0 - ( 1.0 - pow((odtP.Lmin/odtP.Lmax), 5./3.)) * rndd ), -0.6);
    }
    while (odtP.domainLength-eddySize < 2*odtP.domainLength*odtP.dxmin*0.1);
}

///////////////////////////////////////////////////////////////////////////////

/** Uniformly sample the eddy position on the line.  For periodic domains the 
 * position can be anywhere.  For nonperiodic domains, the position 
 * is from 0 to the end where end is the domain size - the eddy size 
 * (since the eddy has to fit in the domain).  This also means that 
 * the eddy position is the left edge of the eddy.
 * For periodic domains, rightEdge is greater than leftEdge (even for eddies that 
 * wrap the domain (since the eddy class is autonomous)), and may be 
 * outside the range of the base odtline.
 * 
 * @param odtP \input parameters object.
 * @param line \input anyline to sample eddies from.
 * @param rr   \input random number generator object.
 */
void eddy::sampleEddyPosition(const odtParam &odtP, anyline &line) {


    if(!odtP.Lperiodic) {
        leftEdge  = rnd->getRand() * (odtP.domainLength-eddySize);
        rightEdge = leftEdge + eddySize;
        // I don't know why, but if one age is equal to a domain boundary, the
        // code gives an segmentation fault.
        while (odtP.domainLength-rightEdge < odtP.domainLength*odtP.dxmin*0.1
                || leftEdge < odtP.domainLength*odtP.dxmin*0.1){
            leftEdge  = rnd->getRand() * (odtP.domainLength-eddySize);
            rightEdge = leftEdge + eddySize;
        }
    }
    else {
        LperiodicEddy = false;             // reset the value
        leftEdge  = rnd->getRand() * odtP.domainLength + line.posf[0];
        rightEdge = leftEdge + eddySize;
        if(rightEdge > (odtP.domainLength + line.posf[0])) {
            LperiodicEddy = true;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

/** Eddy acceptance probability is computed after several previous steps. 
 * 
 * @param odtP \input parameters object.
 * @param dtSample \input eddy sample time step.
 **/
void eddy::computeEddyAcceptanceProb(const odtParam &odtP, double &dtSample) {


    double f, g;

    // f = odtP.esdp4/(eddySize*eddySize) * exp(odtP.esdp1/eddySize);
    // Pa = dtSample * invTauEddy * odtP.C_param / (eddySize * eddySize * f * g);
    // note cancellation of eddySize * eddySize

    f = odtP.esdp4 * exp(odtP.esdp1/eddySize);
    
    if(!odtP.Lperiodic)
        g = 1.0/(odtP.domainLength - eddySize);
    else
        g = 1.0/odtP.domainLength;

    Pa = dtSample * invTauEddy * odtP.C_param / (f * g);

}

///////////////////////////////////////////////////////////////////////////////

/** Apply the triplet map to the eddy.  This is continuous, not discrete.
 * That is, we take the "exact" profile to be the peicewise constant 
 * profiles within cells.  The eddy region is copied three times,
 * Then each copy is compressed by a factor of three, then the middle
 * copy is inserted as a mirror image:
 * \vc{
 * grid 1 :  | {     |     |     |     |     } |       ngrid = 5 (6 faces)
 *                   a     b     c     d
 * grid 2 :  | ( | | | | | | | | | | | | | | ) |       ngrid = 15 (16 faces)
 *               a b c d N d c b a N a b c d           N's line up with ()
 * }
 * Note that the cell is not split at the eddy boundary {, }.  However, later
 * when the eddy is inserted (if accepted), the cell is split there.
 *
 * @param anyl input/output: line object to apply triplet map to.
 * @see anyline::insertEddy
 */
 void eddy::tripMap(anyline &anyl) {
   
    int    ngrd  = anyl.ngrd;
    int    ngrd2 = ngrd*2;
    int    ngrd3 = ngrd*3;
    double cFac   = 1./3.;
    int    i, j, k;
    
    ///////////// make space for 2nd, 3rd segments
    
    anyl.pos.resize(ngrd3);
    anyl.posf.resize(ngrd3+1);
    anyl.rho.resize(ngrd3);
    anyl.molec.resize(ngrd3);
    anyl.lambda.resize(ngrd3);
    anyl.phase.resize(ngrd3);

    for(k=0; k<anyl.nprops; k++) 
        (*anyl.props[k]).resize(ngrd3);

#ifndef COMPSCI
    ///////////// fill second segment

    for(i=ngrd, j=ngrd-1; i<ngrd2; i++, j--) {
        anyl.rho[i]  = anyl.rho[j];
        anyl.molec[i] = anyl.molec[j];
        anyl.lambda[i] = anyl.lambda[j];
        anyl.phase[i] = anyl.phase[j];
        for( k=0; k<anyl.nprops; k++ )
            (*anyl.props[k])[i] = (*anyl.props[k])[j];
    }
    
    ///////////// fill third segment
    
    for(i=ngrd2, j=0; i<ngrd3; i++, j++) {
        anyl.rho[i]  = anyl.rho[j];
        anyl.molec[i] = anyl.molec[j];
        anyl.lambda[i] = anyl.lambda[j];
        anyl.phase[i] = anyl.phase[j];
        for( k=0; k<anyl.nprops; k++ )
            (*anyl.props[k])[i] = (*anyl.props[k])[j];
    }
#else
    //// first alternative
    //std::copy(&anyl.rho[0], &anyl.rho[ngrd], &anyl.rho[ngrd2]);
    //std::reverse_copy(&anyl.rho[0], &anyl.rho[ngrd], &anyl.rho[ngrd]);
    //std::copy(&anyl.molec[0], &anyl.molec[ngrd], &anyl.molec[ngrd2]);
    //std::reverse_copy(&anyl.molec[0], &anyl.molec[ngrd], &anyl.molec[ngrd]);
    //std::copy(&anyl.lambda[0], &anyl.lambda[ngrd], &anyl.lambda[ngrd2]);
    //std::reverse_copy(&anyl.lambda[0], &anyl.lambda[ngrd], &anyl.lambda[ngrd]);
    //std::copy(&anyl.phase[0], &anyl.phase[ngrd], &anyl.phase[ngrd2]);
    //std::reverse_copy(&anyl.phase[0], &anyl.phase[ngrd], &anyl.phase[ngrd]);
    //for( k=0; k<anyl.nprops; k++ ) {
    //    std::copy(&(*anyl.props[k])[0], &(*anyl.props[k])[ngrd], &(*anyl.props[k])[ngrd2]);
    //    std::reverse_copy(&(*anyl.props[k])[0], &(*anyl.props[k])[ngrd], &(*anyl.props[k])[ngrd]);
    //}
    
    // second alternative
    for(int i=0; i<ngrd; i++) { 
        anyl.rho[2*ngrd+i]      = anyl.rho[i];
        anyl.rho[2*ngrd-i-1]    = anyl.rho[i]; 
        anyl.molec[2*ngrd+i]    = anyl.molec[i];
        anyl.molec[2*ngrd-i-1]  = anyl.molec[i]; 
        anyl.lambda[2*ngrd+i]   = anyl.lambda[i];
        anyl.lambda[2*ngrd-i-1] = anyl.lambda[i];
        anyl.phase[2*ngrd+i]    = anyl.phase[i];
        anyl.phase[2*ngrd-i-1]  = anyl.phase[i]; 
    }
    for(int k=0; k<anyl.nprops; k++ ) {
        double* vec = &(*anyl.props[k])[0]; 
        for(int i=0; i<ngrd; i++) { 
            vec[2*ngrd+i]   = vec[i];
            vec[2*ngrd-i-1] = vec[i]; 
        }
    }
#endif
    
    //////////// write new cell and face positions
    
    //--------- cell face positions
    
    double newSegSize = (rightEdge-leftEdge)*cFac;
    
    anyl.posf[ngrd3] = anyl.posf[ngrd];       // last face doesn't change
    
    for(int i=1; i<ngrd; i++)       
        anyl.posf[i] = leftEdge + (anyl.posf[i]-leftEdge)*cFac;
    anyl.posf[ngrd] = leftEdge + newSegSize;               
    
    for(int i=ngrd+1, j=ngrd-1; i<ngrd2; i++, j--) 
        anyl.posf[i] = anyl.posf[ngrd] + (anyl.posf[ngrd]-anyl.posf[j]);
    anyl.posf[ngrd2] = anyl.posf[ngrd] + newSegSize;                

    for(int i=ngrd2+1, j=ngrd2-1; i<ngrd3; i++, j--)
        anyl.posf[i] = anyl.posf[ngrd2] + (anyl.posf[ngrd2]-anyl.posf[j]);


    //--------- cell center positions

    for(int i=0; i<ngrd3; i++)
        anyl.pos[i] = 0.5*(anyl.posf[i]+anyl.posf[i+1]);

    /////////// update ngrd

    anyl.ngrd  = ngrd3;
    anyl.ngrdf = anyl.ngrd+1;

}

///////////////////////////////////////////////////////////////////////////////

/** Fill velocity kernel K (used also for \fun{J=|K|})
 *
 *  @param line \input anyline line object, use for cell positions.
 */

void eddy::fillKernel(anyline &line) {

    ////////// Fill kernel

    int            i;
    int            nseg = line.ngrd / 3;

    K.resize(line.ngrd);

    //---------- 1st Segment

    K[0] = 2*leftEdge - (line.posf[1]+leftEdge);                    // 1st cell

    for(i=1; i<nseg; i++)                                           // other cells
        K[i] = 2.0*(leftEdge - line.pos[i]);

    //---------- Second Segment

    for(i=nseg; i<nseg*2; i++) 
        K[i] = 4.0*(line.pos[i] - leftEdge) - 2.0*eddySize;

    //---------- Third Segment

    for(i=nseg*2; i<line.ngrd-1; i++)                              // all but last
        K[i] = 2.0*(eddySize + leftEdge - line.pos[i]);

    int idmb = line.ngrd-1;                                            
    K[idmb] = 2.0*(eddySize+leftEdge)-(rightEdge+line.posf[idmb]); // last cell

}

///////////////////////////////////////////////////////////////////////////////

/** Apply kernels K and J to the velocity profile.
 * This is called after the kernel coefficients is computed in eddyTau
 *
 * @param line \inout odtline object velocities are updated.
 * @param odtP \input parameters object.
 */
void eddy::applyVelocityKernels(odtline &line, const odtParam &odtP) {

    int            i;
    if(odtP.Llem || odtP.Lspatial) // vanilla odt for spatial formulation
    return;

    ////////// update velocity profiles
    
    if(odtP.LconstProp)
        for(i=0; i<line.ngrd; i++) {
            line.uvel[i] += cCoef[0]*K[i];
            line.vvel[i] += cCoef[1]*K[i];
            line.wvel[i] += cCoef[2]*K[i];
        }
    else 
        for(i=0; i<line.ngrd; i++) {
            line.uvel[i] += cCoef[0]*K[i] + bCoef[0]*fabs(K[i]);
            line.vvel[i] += cCoef[1]*K[i] + bCoef[1]*fabs(K[i]);
            line.wvel[i] += cCoef[2]*K[i] + bCoef[2]*fabs(K[i]);
        }
}

///////////////////////////////////////////////////////////////////////////////

/** Compute the inverse of the eddy timescale. 
 *  Used to define the eddy acceptance probability.
 *  This is the variable property version.  See also eddyTauCP
 *
 *  @param line    \input odtline object used to compute eddy timescale.
 *  @param odtP    \input parameters object.
 *  @param Z_value \input large eddy suppression parameter.
 *  @return false for implausible eddy, true for the usual, proper eddy.
 */
bool eddy::eddyTau(odtline &line, const odtParam &odtP, double Z_value) {

    if(Z_value <= 0) Z_value = odtP.Z_param; //default is suppression coefficient

    double         KK=0;                                  // equivalent of the 4/27 fac
    double         rhoK=0, rhoJ=0, rhoKK=0, rhoJK=0;
    double         uRhoJ[3]={0,0,0}, uRhoK[3]={0,0,0};

    double         P[3];
    double         Q[3];
    double         S;
    double         A;
    double         eKinEddy;        
    double         eViscPenalty;
    double         rhoEddy=0.0, viscEddy=0.0;

    vector<double> dy(line.ngrd);    
    vector<double> intRhoKi(line.ngrd);
    vector<double> intRhoJi(line.ngrd); 

    int            i;
    double         dd;//, dd1, dd2; //  !!!!!  unused variables

    //////////// cell sizes

    dy[0] = line.posf[1]-leftEdge;
    for(i=1; i<line.ngrd-1; i++)
        dy[i] = line.posf[i+1]-line.posf[i];
    dy[line.ngrd-1] = rightEdge-line.posf[line.ngrd-1];


    ///////////// Fill in integral quantities
    
    //---------- rhoK, rhoJ, UrhoK, UrhoJ

    for(i=0; i<line.ngrd; i++) {
        intRhoKi[i]  = K[i]*line.rho[i]*dy[i];
        intRhoJi[i]  = fabs(intRhoKi[i]);
        KK          += K[i]*K[i]*dy[i];         
        if(!odtP.Lspatial) {
            rhoK        += intRhoKi[i];
            rhoJ        += fabs(intRhoKi[i]);
            rhoKK       += K[i]*intRhoKi[i];
            rhoJK       += K[i]*intRhoJi[i];
            uRhoK[0]    += intRhoKi[i]*line.uvel[i];
            uRhoK[1]    += intRhoKi[i]*line.vvel[i];
            uRhoK[2]    += intRhoKi[i]*line.wvel[i];
            uRhoJ[0]    += intRhoJi[i]*line.uvel[i];
            uRhoJ[1]    += intRhoJi[i]*line.vvel[i];
            uRhoJ[2]    += intRhoJi[i]*line.wvel[i];
        }
        else {       // not mass, mass flux
            rhoK        += line.uvel[i] * intRhoKi[i];
            rhoJ        += line.uvel[i] * fabs(intRhoKi[i]);
            rhoKK       += line.uvel[i] * K[i]*intRhoKi[i];
            rhoJK       += line.uvel[i] * K[i]*intRhoJi[i];
            uRhoK[0]    += line.uvel[i] * intRhoKi[i]*line.uvel[i];
            uRhoK[1]    += line.uvel[i] * intRhoKi[i]*line.vvel[i];
            uRhoK[2]    += line.uvel[i] * intRhoKi[i]*line.wvel[i];
            uRhoJ[0]    += line.uvel[i] * intRhoJi[i]*line.uvel[i];
            uRhoJ[1]    += line.uvel[i] * intRhoJi[i]*line.vvel[i];
            uRhoJ[2]    += line.uvel[i] * intRhoJi[i]*line.wvel[i];
        }
    }

    KK /= eddySize*eddySize*eddySize; // Exact model formulation gives KK=4/27

    //////////// Compute Eddy Energy

    A = rhoK/rhoJ;
    S = 0.5*(A*A+1.0)*rhoKK - A*rhoJK;
    for(i=0; i<3; i++) {
        P[i] = uRhoK[i] - A * uRhoJ[i];
        Q[i] = 0.25*P[i]*P[i]/S;
    }

    eKinEddy = KK * (Q[0]+Q[1]+Q[2]); 
    
    //////////// Compute Additional Surface Tension

    double eSurfTens = odtP.eSurfTens;

    if(odtP.eSurfTens != 0.0) {
        int iSurfTens = 0;
        LmultiPhaseEddy = false;
        
        for(i=0; i<line.ngrd-1; i++) {   
            if(line.phase[i] != line.phase[i+1])
                //return false; // temporary added by Falko
                iSurfTens++; //USE POINTER TO line.phase
        }
        
        eSurfTens = odtP.eSurfTens * KK * (iSurfTens * 2 * 2. / 3.); 
        // KK for eKinEddy scaling, 2 for isotropic surface/vol, 2/3 for new interfaces
        
        if (iSurfTens > 0){
            LmultiPhaseEddy = true;
            //cout << endl << "Falko: Verhaeltnis = " << eSurfTens / eKinEddy << endl;
        } 
    }

    //////////// Compute Viscous Energy Penalty

    for(i=0; i<line.ngrd; i++) { 
        rhoEddy +=line.rho[i]*dy[i];
        viscEddy += dy[i]/line.molec[i];
    }
    rhoEddy  /= eddySize;
    viscEddy  = eddySize/viscEddy;

    double Ufavre;

    eViscPenalty = Z_value*0.5 * viscEddy*viscEddy/rhoEddy/eddySize; 

    if(odtP.Lspatial)  {
        Ufavre = eddyFavreAvgVelocity(line);
        eViscPenalty *= Ufavre;
    }

    //////////// Compute invTauEddy
    //if (eSurfTens > 0){
    //    cout << endl << "Falko: Verhaeltnis = " << eSurfTens / eKinEddy << " " << eKinEddy - eSurfTens - eViscPenalty << " " << sqrt(2.0*KK/rhoKK * (eKinEddy - eSurfTens - eViscPenalty)) << endl;
    //}
    
    invTauEddy = 0.0;
    if(odtP.LheatedChannel){
        // FALKO debug: hier muss noch der Auftriebsterm hinzugefügt werden.
        double ePotEddy = 0.0;
        dd = eKinEddy - eSurfTens - ePotEddy - eViscPenalty;
    }
    else
        dd = eKinEddy - eSurfTens - eViscPenalty;

    if(dd < 0.0) 
        return false;
    invTauEddy = sqrt(2.0*KK/rhoKK * dd);

    if(odtP.Lspatial) 
        invTauEddy = invTauEddy/Ufavre; // 1/s --> 1/m


    //////////// Compute Kernel Coefficients (THIS IS WASTEFUL)

    if(odtP.Lspatial)
        return true;               // vanilla odt for spatial formulation

    cCoef[0] = 0.5/S * (-P[0] + (P[0]>0 ? 1.0 : -1.0) 
                 * sqrt( (1-odtP.A_param)*P[0]*P[0]  
                         //+ 0.5*odtP.A_param*(P[1]*P[1]+P[2]*P[2]) ));
                         + 0.5*odtP.A_param*(P[1]*P[1]+P[2]*P[2]) - S*eSurfTens/3.0));
    cCoef[1] = 0.5/S * (-P[1] + (P[1]>0 ? 1.0 : -1.0) 
                 * sqrt( (1-odtP.A_param)*P[1]*P[1]  
                         //+ 0.5*odtP.A_param*(P[0]*P[0]+P[2]*P[2]) ));
                         + 0.5*odtP.A_param*(P[0]*P[0]+P[2]*P[2]) - S*eSurfTens/3.0));
    cCoef[2] = 0.5/S * (-P[2] + (P[2]>0 ? 1.0 : -1.0) 
                 * sqrt( (1-odtP.A_param)*P[2]*P[2]  
                         //+ 0.5*odtP.A_param*(P[0]*P[0]+P[1]*P[1]) ));
                         + 0.5*odtP.A_param*(P[0]*P[0]+P[1]*P[1]) -S*eSurfTens/3.0));
                
    for(i=0; i<3; i++)
        bCoef[i] = -A*cCoef[i];

    return true;

}

///////////////////////////////////////////////////////////////////////////////

/** Compute the inverse of the eddy timescale. 
 *  Used to define the eddy acceptance probability.
 *  This is the variable property version.  See also eddyTauCP
 *
 *  @param line    \input odtline object used to compute eddy timescale.
 *  @param odtP    \input parameters object.
 *  @param Z_value \input large eddy suppression parameter.
 *  @return false for implausible eddy, true for the usual, proper eddy.
 */
bool eddy::eddyTauPartSrc(odtline &line, const odtParam &odtP, double Z_value) {

    if(Z_value <= 0) Z_value = odtP.Z_param; //default is suppression coefficient

    double         KK=0;                                  // equivalent of the 4/27 fac
    double         rhoK=0, rhoJ=0, rhoKK=0, rhoJK=0;
    double         uRhoJ[3]={0,0,0}, uRhoK[3]={0,0,0};

    double         P[3];
    double         Q[3];
    double         S;
    double         A;
    double         eKinEddy;        
    double         eViscPenalty;
    double         rhoEddy=0.0, viscEddy=0.0;

    vector<double> dy(line.ngrd);    
    vector<double> intRhoKi(line.ngrd);
    vector<double> intRhoJi(line.ngrd); 

    int            i;
    double         dd, dd1, dd2;

    //////////// cell sizes

    dy[0] = line.posf[1]-leftEdge;
    for(i=1; i<line.ngrd-1; i++)
        dy[i] = line.posf[i+1]-line.posf[i];
    dy[line.ngrd-1] = rightEdge-line.posf[line.ngrd-1];


    ///////////// Fill in integral quantities
    
    //---------- rhoK, rhoJ, UrhoK, UrhoJ

    for(i=0; i<line.ngrd; i++) {
        intRhoKi[i]  = K[i]*line.rho[i]*dy[i];
        intRhoJi[i]  = fabs(intRhoKi[i]);
        KK          += K[i]*K[i]*dy[i];         
        if(!odtP.Lspatial) {
            rhoK        += intRhoKi[i];
            rhoJ        += fabs(intRhoKi[i]);
            rhoKK       += K[i]*intRhoKi[i];
            rhoJK       += K[i]*intRhoJi[i];
            uRhoK[0]    += intRhoKi[i]*line.uvel[i];
            uRhoK[1]    += intRhoKi[i]*line.vvel[i];
            uRhoK[2]    += intRhoKi[i]*line.wvel[i];
            uRhoJ[0]    += intRhoJi[i]*line.uvel[i];
            uRhoJ[1]    += intRhoJi[i]*line.vvel[i];
            uRhoJ[2]    += intRhoJi[i]*line.wvel[i];
        }
        else {       // not mass, mass flux
            rhoK        += line.uvel[i] * intRhoKi[i];
            rhoJ        += line.uvel[i] * fabs(intRhoKi[i]);
            rhoKK       += line.uvel[i] * K[i]*intRhoKi[i];
            rhoJK       += line.uvel[i] * K[i]*intRhoJi[i];
            uRhoK[0]    += line.uvel[i] * intRhoKi[i]*line.uvel[i];
            uRhoK[1]    += line.uvel[i] * intRhoKi[i]*line.vvel[i];
            uRhoK[2]    += line.uvel[i] * intRhoKi[i]*line.wvel[i];
            uRhoJ[0]    += line.uvel[i] * intRhoJi[i]*line.uvel[i];
            uRhoJ[1]    += line.uvel[i] * intRhoJi[i]*line.vvel[i];
            uRhoJ[2]    += line.uvel[i] * intRhoJi[i]*line.wvel[i];
        }
    }

    KK /= eddySize*eddySize*eddySize; // Exact model formulation gives KK=4/27

    //////////// Compute Eddy Energy
        
    double AA;
    vector<double> BB(3,0.0);
    vector<double> CC(3,0.0);

        AA = 0.5*rhoKK + 0.5*rhoK*rhoK*rhoKK/rhoJ/rhoJ - 2*rhoK*rhoJK/rhoJ;
        for(i=0; i<3; i++) {
            BB[i] = partMomSrcInEddy[i]*rhoK*rhoKK/rhoJ/rhoJ - partMomSrcInEddy[i]*rhoJK/rhoJ - rhoK/rhoJ*uRhoJ[i] + uRhoK[i];    
            CC[i] = 0.5*rhoKK*partMomSrcInEddy[i]*partMomSrcInEddy[i]/rhoJ/rhoJ - partMomSrcInEddy[i]*uRhoJ[i]/rhoJ;
            Q[i] = -(CC[i] - BB[i]*BB[i]/4./AA);
        } 

        eKinEddy = KK*(Q[0]+Q[1]+Q[2]);            

    //////////// Compute Additional Surface Tension

    double eSurfTens = 0.0;

    if(odtP.eSurfTens != 0.0) {
        int iSurfTens = 0;

        for(i=0; i<line.ngrd-1; i++) {   
            if(line.phase[i] != line.phase[i+1]) 
                iSurfTens++; //USE POINTER TO line.phase
        }

        eSurfTens = odtP.eSurfTens * KK * (iSurfTens * 2 * 2. / 3.); 
        // KK for eKinEddy scaling, 2 for isotropic surface/vol, 2/3 for new interfaces
    }

    //////////// Compute Viscous Energy Penalty

    for(i=0; i<line.ngrd; i++) { 
        rhoEddy +=line.rho[i]*dy[i];
        viscEddy += dy[i]/line.molec[i];
    }
    rhoEddy  /= eddySize;
    viscEddy  = eddySize/viscEddy;

    double Ufavre;

    eViscPenalty = Z_value*0.5 * viscEddy*viscEddy/rhoEddy/eddySize; 

    if(odtP.Lspatial)  {
        Ufavre = eddyFavreAvgVelocity(line);
        eViscPenalty *= Ufavre;
    }

    //////////// Compute invTauEddy

    invTauEddy = 0.0;    
    dd = eKinEddy - eSurfTens - eViscPenalty;

    if(dd < 0.0) 
        return false;
    invTauEddy = sqrt(2.0*KK/rhoKK * dd);

    //////////// Compute Kernel Coefficients (THIS IS WASTEFUL)

    cCoef[0] = (-BB[0] + (BB[0]>0 ? 1.0 : -1.0)
            * sqrt(BB[0]*BB[0]-4*AA*(CC[0]+odtP.A_param*(Q[0]-0.5*Q[1]-0.5*Q[2])) - AA*eSurfTens/3.0))*0.5/AA;

    cCoef[1] = (-BB[1] + (BB[1]>0 ? 1.0 : -1.0)
            * sqrt(BB[1]*BB[1]-4*AA*(CC[1]+odtP.A_param*(Q[1]-0.5*Q[0]-0.5*Q[2])) - AA*eSurfTens/3.0))*0.5/AA;

    cCoef[2] = (-BB[2] + (BB[2]>0 ? 1.0 : -1.0)
            * sqrt(BB[2]*BB[2]-4*AA*(CC[2]+odtP.A_param*(Q[2]-0.5*Q[0]-0.5*Q[1])) - AA*eSurfTens/3.0))*0.5/AA;

    for(i=0; i<3; i++) {
        bCoef[i] = (-partMomSrcInEddy[i]-cCoef[i]*rhoK)/rhoJ;
    }
    return true;

}

///////////////////////////////////////////////////////////////////////////////

/** This is the constant property version of eddyTau, which is simpler and faster.
*  Does not have the surface tension implemented.
* 
*  @param line    \input odtline object.
*  @param odtP    \input parameters object.
*  @param Z_value \input large eddy suppression parameter.
*  @return false for implausible eddy, true for the usual, proper eddy.
*/
bool eddy::eddyTauCP(odtline &line, const odtParam &odtP, double Z_value) {


    if(Z_value <= 0) Z_value = odtP.Z_param; //default is suppression coefficient

    double         uK[3]={0,0,0};
#ifndef COMPSCI
    vector<double> dy(line.ngrd);
#else
    eddyTauDy.resize(line.ngrd);
    vector<double>& dy = eddyTauDy;
#endif

    double         KK = 0;            // equilvalent of the 4/27 factor

    int            i;

    //////////// cell sizes

    dy[0] = line.posf[1]-leftEdge;
    for(i=1; i<line.ngrd-1; i++)
        dy[i] = line.posf[i+1]-line.posf[i];
    dy[line.ngrd-1] = rightEdge-line.posf[line.ngrd-1];

    ///////////// Fill in integral quantities
    
    for(i=0; i<line.ngrd; i++) {                 //---------- K kernels: whole eddy
        uK[0] += K[i]*line.uvel[i]*dy[i]; 
        uK[1] += K[i]*line.vvel[i]*dy[i]; 
        uK[2] += K[i]*line.wvel[i]*dy[i]; 
        KK    += K[i]*K[i]*dy[i];         
    }

    //////////// Compute invTauEddy

    double uk = uK[0]/eddySize/eddySize;
    double vk = uK[1]/eddySize/eddySize;
    double wk = uK[2]/eddySize/eddySize;

    invTauEddy = 0.0;    

    double invTauSq;
    double eKinEddy     = uk*uk+vk*vk+wk*wk;
    double eViscPenalty = Z_value*odtP.visc_0*odtP.visc_0 /
                          (odtP.rho_0*odtP.rho_0*eddySize*eddySize);
    if(odtP.LheatedChannel){
        // FALKO debug: hier muss noch der Auftriebsterm hinzugefügt werden.
        //double rhoK         = 0.0;
        //for(i=0; i<line.ngrd; i++) {
        //    rhoK        += (line.rho[i]-odtP.rho_0) * (eddySize - 2*(line.pos[i] - rightEdge) )*dy[i];
        //rhoK = rhoK * 4.0 / 9.0 / eddySize;
        double ePotEddy     = 0.0; //8.0 / 27.0 * rhoK * odtP.Grav / odtP.rho_0;
        invTauSq = 1.0/(eddySize*eddySize) * (eKinEddy - ePotEddy - eViscPenalty);
    }
    else {
        invTauSq = 1.0/(eddySize*eddySize) * (eKinEddy - eViscPenalty);
    }

    if(invTauSq < 0.0) 
        return false; 
    invTauEddy = sqrt(invTauSq);
    
    if(odtP.Lspatial)
        invTauEddy = invTauEddy/eddyFavreAvgVelocity(line); // 1/s --> 1/m
    
    //////////// Compute Kernel Coefficients (THIS IS WASTEFUL)

    if(odtP.Lspatial)
        return true;               // vanilla odt for spatial formulation

    double cFac = eddySize*eddySize/KK; //AKeddy based on numerical kernel
    
    cCoef[0] = cFac * ( -uk + (uk>0 ? 1.0 : -1.0)
                              * sqrt( (1.0-odtP.A_param)*uk*uk
                                      + 0.5*odtP.A_param*(vk*vk + wk*wk) ) ); 

    cCoef[1] = cFac * ( -vk + (vk>0 ? 1.0 : -1.0)
                              * sqrt( (1.0-odtP.A_param)*vk*vk
                                      + 0.5*odtP.A_param*(uk*uk + wk*wk) ) ); 
    cCoef[2] = cFac * ( -wk + (wk>0 ? 1.0 : -1.0)
                              * sqrt( (1.0-odtP.A_param)*wk*wk
                                      + 0.5*odtP.A_param*(uk*uk + vk*vk) ) ); 

    return true; 

}

///////////////////////////////////////////////////////////////////////////////

/** This is the constant property version of eddyTau, which is simpler and faster.
*  Does not have the surface tension implemented.
* 
*  @param line    \input odtline object.
*  @param odtP    \input parameters object.
*  @param Z_value \input large eddy suppression parameter.
*  @return false for implausible eddy, true for the usual, proper eddy.
*/
bool eddy::eddyTauCPpartSrc(odtline &line, const odtParam &odtP, double Z_value) {


    if(Z_value <= 0) Z_value = odtP.Z_param; //default is suppression coefficient

    double         uJ[3]={0,0,0}, uK[3]={0,0,0};
    vector<double> dy(line.ngrd);    
    vector<double> intKi(line.ngrd);
    vector<double> intJi(line.ngrd); 

    double         sumintJ = 0, KK = 0, KJ = 0;            // equilvalent of the 4/27 factor

    int            i;
    
    double         Q[3]={0,0,0};

    for(i=0; i<3; i++) {
        bCoef[i] = 0;
        cCoef[i] = 0;
    }

    //////////// cell sizes

    dy[0] = line.posf[1]-leftEdge;
    for(i=1; i<line.ngrd-1; i++)
        dy[i] = line.posf[i+1]-line.posf[i];
    dy[line.ngrd-1] = rightEdge-line.posf[line.ngrd-1];

    ///////////// Fill in integral quantities


// cout << endl << fixed << setprecision(12) << "eddy posf ngrd " << line.posf[line.ngrd] << " rightEdge " << rightEdge << endl;
// cout << endl << fixed << setprecision(12) << "odtP.rho_0 " << odtP.rho_0 << endl;

    for(i=0; i<line.ngrd; i++) {
        intKi[i]  = K[i]*dy[i]/eddySize/eddySize;
        intJi[i]  = fabs(intKi[i]);
        sumintJ  += intJi[i];
        KK       += K[i]*K[i]*dy[i]/eddySize/eddySize/eddySize;         
        KJ       += K[i]*fabs(K[i])*dy[i]/eddySize/eddySize/eddySize;   
        uK[0]    += intKi[i]*line.uvel[i];
        uK[1]    += intKi[i]*line.vvel[i];
        uK[2]    += intKi[i]*line.wvel[i];
        uJ[0]    += intJi[i]*line.uvel[i];
        uJ[1]    += intJi[i]*line.vvel[i];
        uJ[2]    += intJi[i]*line.wvel[i];
    }

    //////////// particle coupling
    
    double AAhat;
    vector<double> BBhat(3,0.0);
    vector<double> CChat(3,0.0);
    vector<double> bCoefhat(3,0.0);
    vector<double> cCoefhat(3,0.0);
    vector<double> uvelhat(line.ngrd,0.0);
    vector<double> vvelhat(line.ngrd,0.0);
    vector<double> wvelhat(line.ngrd,0.0);

    AAhat = 0.5*odtP.rho_0*KK*eddySize*eddySize*eddySize; 

// cout << endl << "sumintJ " << sumintJ * eddySize *eddySize << endl;

    for(i=0; i<3; i++) {
        
        bCoefhat[i] = -partMomSrcInEddy[i]/odtP.rho_0/sumintJ/eddySize/eddySize; 
        BBhat[i] = (odtP.rho_0*bCoefhat[i]*KJ*eddySize + odtP.rho_0*uK[i])*eddySize*eddySize;
        CChat[i] = (0.5*odtP.rho_0*bCoefhat[i]*bCoefhat[i]*KK*eddySize + odtP.rho_0*bCoefhat[i]*uJ[i])*eddySize*eddySize + partErgSrcInEddy[i];
        cCoefhat[i] = (-BBhat[i] + (BBhat[i]>0 ? 1.0 : -1.0) * sqrt(BBhat[i]*BBhat[i]-4*AAhat*CChat[i]))*0.5/AAhat;
        
        if((BBhat[i]*BBhat[i]-4*AAhat*CChat[i]) < 0.) return false; 

    }

    for(i=0; i<line.ngrd; i++) {
        uvelhat[i] = line.uvel[i] + cCoefhat[0]*K[i] + bCoefhat[0]*fabs(K[i]);
        vvelhat[i] = line.vvel[i] + cCoefhat[1]*K[i] + bCoefhat[1]*fabs(K[i]);
        wvelhat[i] = line.wvel[i] + cCoefhat[2]*K[i] + bCoefhat[2]*fabs(K[i]);
    } 
    
    ///////////// (After adding particle source) Fill in integral quantities AGAIN

    for(i=0; i<3; i++) {
        uK[i] = 0;
    }

    for(i=0; i<line.ngrd; i++) {
        uK[0]    += intKi[i]*uvelhat[i];
        uK[1]    += intKi[i]*vvelhat[i];
        uK[2]    += intKi[i]*wvelhat[i];
    }
    
    //////////// Compute Kernel Coefficients (THIS IS WASTEFUL)
   
    double AA;
    vector<double> BB(3,0.0);
    
    AA = 0.5*odtP.rho_0*KK*eddySize*eddySize*eddySize; 
    
    for(i=0; i<3; i++) {
        
        BB[i] = odtP.rho_0*uK[i]*eddySize*eddySize;
        Q[i]  = BB[i]*BB[i]/4./AA;
    }

    cCoef[0] = (-BB[0] + (BB[0]>0 ? 1.0 : -1.0)
            * sqrt(BB[0]*BB[0]-4*AA*(odtP.A_param*(Q[0]-0.5*Q[1]-0.5*Q[2]))))*0.5/AA;

    cCoef[1] = (-BB[1] + (BB[1]>0 ? 1.0 : -1.0)
            * sqrt(BB[1]*BB[1]-4*AA*(odtP.A_param*(Q[1]-0.5*Q[0]-0.5*Q[2]))))*0.5/AA;

    cCoef[2] = (-BB[2] + (BB[2]>0 ? 1.0 : -1.0)
            * sqrt(BB[2]*BB[2]-4*AA*(odtP.A_param*(Q[2]-0.5*Q[0]-0.5*Q[1]))))*0.5/AA;

    
//     double cCoeftest0 = 1/KK/eddySize * ( -uK[0] + (uK[0]>0 ? 1.0 : -1.0)
//                               * sqrt( (1.0-odtP.A_param)*uK[0]*uK[0]
//                                       + 0.5*odtP.A_param*(uK[1]*uK[1] + uK[2]*uK[2]) ) ); 
// 
//     double cCoeftest1 = 1/KK/eddySize * ( -uK[1] + (uK[1]>0 ? 1.0 : -1.0)
//                               * sqrt( (1.0-odtP.A_param)*uK[1]*uK[1]
//                                       + 0.5*odtP.A_param*(uK[0]*uK[0] + uK[2]*uK[2]) ) ); 
//     double cCoeftest2 = 1/KK/eddySize * ( -uK[2] + (uK[2]>0 ? 1.0 : -1.0)
//                               * sqrt( (1.0-odtP.A_param)*uK[2]*uK[2]
//                                       + 0.5*odtP.A_param*(uK[0]*uK[0] + uK[1]*uK[1]) ) ); 
// 
// cout << endl << fixed << setprecision(6) << "cCoef[0] " << cCoef[0] << endl;
// cout << endl << fixed << setprecision(6) << "cCoef[1] " << cCoef[1] << endl;
// cout << endl << fixed << setprecision(6) <<"cCoef[2] " << cCoef[2] << endl;
// cout << endl << fixed << setprecision(6) <<"cCoeftest[0] " << cCoeftest0 << endl;
// cout << endl << fixed << setprecision(6) <<"cCoeftest[1] " << cCoeftest1 << endl;
// cout << endl << fixed << setprecision(6) <<"cCoeftest[2] " << cCoeftest2 << endl;

    bCoef[0] = bCoefhat[0];
    bCoef[1] = bCoefhat[1];
    bCoef[2] = bCoefhat[2];
    cCoef[0] = cCoef[0] + cCoefhat[0];
    cCoef[1] = cCoef[1] + cCoefhat[1];
    cCoef[2] = cCoef[2] + cCoefhat[2];

    //////////// Compute invTauEddy
    
    invTauEddy = 0.0;    

    double invTauSq;
    
    double eKinEddy = KK*(Q[0]+Q[1]+Q[2]);            

    double eViscPenalty = Z_value*odtP.visc_0*odtP.visc_0 /
                          (odtP.rho_0*eddySize)*0.5;
    invTauSq = 2.0/(odtP.rho_0*eddySize*eddySize*eddySize) * (eKinEddy - eViscPenalty);

    if(invTauSq < 0.0) 
        return false; 
    invTauEddy = sqrt(invTauSq);

    return true;

}

///////////////////////////////////////////////////////////////////////////////

/** Get the Favre average streamwise velocity for spatial formulations.
*  Used to convert from eddy timescale to the eddy spatial scale.
* 
*  @param line \input odtline object.
*  @return value of the favre avg velocity
*/
double eddy::eddyFavreAvgVelocity(odtline &line) {

    double ufavg = 0.0;
    double ravg  = 0.0;
    for(int i=0; i<line.ngrd; i++) {
        ufavg += line.rho[i]*line.uvel[i];
        ravg  += line.rho[i];
    }
    return ufavg/ravg;
}

///////////////////////////////////////////////////////////////////////////////

/**  Get eddy velocity v on the odtline through tripletMap of lagrangian tracer particles.
*  (Only moves tracer particles, if the tracer flag is on.  Also computes veddy, needed below).
*  \vc{
*  before TM   |   { a }   |   { b }   |   { c }   |   
*  after  TM   | a | b | c | c | b | a | a | b | c |
*  }
*  for example, if particle locates in the 1st third of b domain before TM, it will move to
*  the 1st b domain after TM. The fraction of position in the cell will NOT change.
*
*  Moves position of tracer particles.
*  For inertial particles, move particles as if tracers in order to compute vEddy for each particle.
*  At the end of the function (for inertial particles), reset particle positions to initial.
*  
*  @param odtP      \input odt Parameters object
*  @param line      \input odtline object
*  @param eddyLine  \input eddy line (also an odtline object)
*  @param part      \inout particles object
*  @param eddyStartTime \inout eddy Start Time for particle histories
*
*  @author Guangyuan Sun 08/2011
*/
void eddy::getEddyVvel_TMtracers(const odtParam &odtP, odtline &line, odtline &eddyLine, particles &part, double eddyStartTime) {
    
    if (!odtP.Iparticles) 
        return;

    vector<double> yPos0 = part.yPos;   
    vEddy.resize(part.nPart);
    double f13 = 1./3;
    double f23 = 2./3;

    for(int i=0; i<part.nPart; i++) {                

        if (!part.pActive[i]) 
            continue;
        
        if (LperiodicEddy) {
            if((part.yPos[i] >= line.posf[0]) && (part.yPos[i] <= (rightEdge-line.Ldomain))) {
                part.yPos[i] += line.Ldomain;
		        if (part.Ltracer) part.crossBound[i]--;
            }
        }
        
        double yPosInEddy = part.yPos[i];
        
        if ((part.yPos[i] >= leftEdge) && (part.yPos[i] <= rightEdge)) {            // only does particles in the eddy region
            
//            double iyPosInEddy = eddyLine.linePositionToIndex(part.yPos[i],true);
//            double fracC = (part.yPos[i] - eddyLine.posf[iyPosInEddy])  / 
//                (eddyLine.posf[iyPosInEddy+1] - eddyLine.posf[iyPosInEddy]);
            /********************************************************************************************/
//            part.set_iyPos();
//            part.iyPos[i] =  line.linePositionToIndex(part.yPos[i], true);
//            double fracC1 = (part.yPos[i] - line.posf[part.iyPos[i]])  / 
//                (line.posf[part.iyPos[i]+1] - line.posf[part.iyPos[i]]);
            /********************************************************************************************/
//            double fracL = (part.yPos[i] - leftEdge) / eddySize;
//            if ((fracC1 >= 0.0) && (fracC1 < f13)) 
//                part.yPos[i] = leftEdge + fracL * (f13*eddySize);    
//            else if ((fracC1 >= f13) && (fracC1 < f23))
//                part.yPos[i] = leftEdge + (2.-fracL)*(f13*eddySize);
//            else
//                part.yPos[i] = leftEdge + (2.+fracL)*(f13*eddySize);
            
            double fracL = (part.yPos[i] - leftEdge) / eddySize;
            
            double r = rnd->getRand();

            if (r >= 0.0 && r < f13) 
                part.yPos[i] = leftEdge + fracL * (f13*eddySize);    
            else if (r >= f13 && r < f23)
                part.yPos[i] = leftEdge + (2.-fracL)*(f13*eddySize);
            else
                part.yPos[i] = leftEdge + (2.+fracL)*(f13*eddySize);

            vEddy[i] = (part.yPos[i]-yPosInEddy)/eddyLife[i];

            if(LperiodicEddy)
                if(part.yPos[i] >= line.Ldomain){
                    part.yPos[i] -= line.Ldomain;
                    if (part.Ltracer) part.crossBound[i]++;
                }
            // JCH ------------------- adjust particle histories for Tracer particles    
            if ( part.Lhistories && part.Ltracer ) {      
                part.adjustEddyVelHistoriesTracers( i, eddyStartTime, eddyLife[i], vEddy[i] );
            }
        }

        part.yPosTracer_TM[i] = part.yPos[i];

        if (!part.Ltracer) {   ///inertia particle, new position and velocity will be calculated in eddy::getParticleUVWYafterEddy
            part.yPos[i] = yPos0[i];
            part.iyPos[i] = line.linePositionToIndex(part.yPos[i],true);
        }
    }
}


///////////////////////////////////////////////////////////////////////////////

/** Get eddy velocities u and w
*  by calculating the average of the velocity of each cell in the whole domain of eddy.
*  Note that eddy velocity v is fluid velocity in y direction seen by the particles during particle-eddy interaction
*  Therefore, eddy velocity v is different for different particles that interact with the same eddy
*  eddy velocity v is calculated in eddy::tripletMapParticle function
*
*  @author Guangyuan Sun 08/2011
*/

void eddy::getEddyUWvel(odtline &line, particles &part) {


    int iEdLeft  = line.linePositionToIndex(leftEdge,true);
    int iEdRight = line.linePositionToIndex(rightEdge,true);
    uEddyAvg = line.uvel[iEdLeft]*(line.posf[iEdLeft+1]-leftEdge) + line.uvel[iEdRight]*(rightEdge-line.posf[iEdRight]);
    vEddyAvg = line.vvel[iEdLeft]*(line.posf[iEdLeft+1]-leftEdge) + line.vvel[iEdRight]*(rightEdge-line.posf[iEdRight]);
    wEddyAvg = line.wvel[iEdLeft]*(line.posf[iEdLeft+1]-leftEdge) + line.wvel[iEdRight]*(rightEdge-line.posf[iEdRight]);
    for(int i=iEdLeft+1; i<=iEdRight-1; i++) {
        uEddyAvg = uEddyAvg + line.uvel[i]*(line.posf[i+1]-line.posf[i]);
        vEddyAvg = vEddyAvg + line.vvel[i]*(line.posf[i+1]-line.posf[i]);
        wEddyAvg = wEddyAvg + line.wvel[i]*(line.posf[i+1]-line.posf[i]);
    }
    uEddyAvg = uEddyAvg/eddySize;
    vEddyAvg = vEddyAvg/eddySize;
    wEddyAvg = wEddyAvg/eddySize;

    uEddy.resize(part.nPart,0.0);
    wEddy.resize(part.nPart,0.0);

//GYSun
//if(part.PeddyType == 2) {
//    if((part.yPos[0] >= leftEdge) && (part.yPos[0] <= rightEdge)) {
//        int iyPosInEddy = line.linePositionToIndex(part.yPos[0], true);
//        uEddyAvg = line.uvel[iyPosInEddy];
//        wEddyAvg = line.wvel[iyPosInEddy];
//    }
//}

    for(int i=0; i < part.nPart; i++) {
        if (!part.pActive[i]) {
            uEddy[i] = 0.0;
            wEddy[i] = 0.0;
            continue;                                     // skip inactive particles (e.g. wall collisions/outflow)
        }
        if ((part.yPos[i] >= leftEdge) && (part.yPos[i] <= rightEdge)) {            // only does particles in the eddy region
            int iyPosInEddy = line.linePositionToIndex(part.yPos[i],true);
            
            if(part.PeddyType == 1 || part.PeddyType == 3) { // ----- TypeI
                uEddy[i] = line.uvel[iyPosInEddy];
                wEddy[i] = line.wvel[iyPosInEddy];
            }
        }
    }

}

//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv Helper functions for eddy::getParticleUVWYafterEddy

double eddy::pPosRelToEdge(double Evel, 
                           double Gvel, 
                           double Pvel, 
                           double tauP, 
                           double f, 
                           double AG, 
                           double pPos0,
                           double edgePos0, 
                           double T) {

    // return (particle position) - (edge position)
    
    return pPos0 + (Gvel+AG*tauP/f)*T + tauP/f*(Pvel-Gvel-AG*tauP/f)*(1.0-exp(-T*f/tauP)) 
           - (edgePos0 + Evel*T);
}       

//-----------------------------------------------------------------------------------------

double eddy::ddt_pPosRelToEdge(double Evel, 
                               double Gvel, 
                               double Pvel, 
                               double tauP, 
                               double f, 
                               double AG, 
                               double T) {

    // return derivative with respect to time of [(particle position) - (edge position)]
    
    return (Gvel+AG*tauP/f) +(Pvel-Gvel-AG*tauP/f)*exp(-T*f/tauP) 
           - (Evel);
}       
//-----------------------------------------------------------------------------------------

double eddy::getEddyCrossingTime(double Evel, 
                                 double Gvel, 
                                 double Pvel, 
                                 double tauP, 
                                 double f, 
                                 double AG, 
                                 double pPos0,
                                 double edgePos0, 
                                 double T1, 
                                 double T2) {

    double Tcross;
    double f1 = pPosRelToEdge(Evel, Gvel, Pvel, tauP, f, AG, pPos0, edgePos0, T1);
    double f2 = pPosRelToEdge(Evel, Gvel, Pvel, tauP, f, AG, pPos0, edgePos0, T2);
    double Tn, fn, fp;
    if(f1*f2 > 0.0) {
        *proc.ostrm << endl << "# ERROR in eddy::getParticleUVWYafterEddy x left: no bracket" << endl;
        exit(0);
    }
    for(int it=1; it<=6; it++) {    
        Tn = (T1+T2)*0.5;
        fn = pPosRelToEdge(Evel, Gvel, Pvel, tauP, f, AG, pPos0, edgePos0, Tn);
        (f1*fn < 0.0) ? (T2 = Tn) : (T1=Tn);
    }   

    //----------- clean up solution with Newton

    int maxIt=1000;
    for(int it=1; it<=maxIt; it++) {
        fp = ddt_pPosRelToEdge(Evel, Gvel, Pvel, tauP, f, AG, Tn);
        Tn = Tn - fn/fp;
        fn = pPosRelToEdge(Evel, Gvel, Pvel, tauP, f, AG, pPos0, edgePos0, Tn);
        if( abs(fn/eddySize) < 1E-8 ) {
            Tcross = Tn;
            break;
        }
        if(it==maxIt) 
            *proc.ostrm << endl << "# warning in eddy::getParticleUVWYafterEddy, not converged" << endl;
    }

    return Tcross;

}

//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Helper functions for eddy::getParticleUVWYafterEddy

//-----------------------------------------------------------------------------------------

/* JCH -- For the X and Z positions, we might consider setting the particle at the upwind edge of the eddy instead of halfway.  This is because the transverse velocity correlation is half Lf for isotropic stationary flow.  Is this the right way of thinking about this?
 */

/// Helper function for eddy::getParticleUVWYafterEddy
double eddy::X1(double Euvel, double Guvel, double Puvel, double Tau, double f, double AGx, double T) {
    return -Euvel*T + 0.5*eddySize+(Guvel+AGx*Tau/f)*T+Tau/f*(Puvel-Guvel-AGx*Tau/f)*(1-exp(-T*f/Tau));
}
/// Helper function for eddy::getParticleUVWYafterEddy
double eddy::X2(double Euvel, double Guvel, double Puvel, double Tau, double f, double AGx, double T) {
    return -Euvel*T - 0.5*eddySize+(Guvel+AGx*Tau/f)*T+Tau/f*(Puvel-Guvel-AGx*Tau/f)*(1-exp(-T*f/Tau));
}
/// Helper function for eddy::getParticleUVWYafterEddy
double eddy::derX(double Euvel, double Guvel, double Puvel, double Tau, double f, double AGx, double T) {
    return -Euvel+(Guvel+AGx*Tau/f)+(Puvel-Guvel-AGx*Tau/f)*exp(-T*f/Tau);
}
/// Helper function for eddy::getParticleUVWYafterEddy  GSDB
double eddy::Y1(double PyPos, double Evvel, double Pvvel, double Tau, double f, double AGy, double T) {
    return -leftEdge+PyPos+(Evvel+AGy*Tau/f)*T+Tau/f*(Pvvel-Evvel-AGy*Tau/f)*(1-exp(-T*f/Tau));
}
// double eddy::Y1(double PyPos, double Pvvel, double velEd, double Tau, double f, double AGy, double T) {
//     return -leftEdge+PyPos+velEd*(T-Tau/f*(1-exp(-T*f/Tau)));
// }

/// Helper function for eddy::getParticleUVWYafterEddy GSDB
double eddy::Y2(double PyPos, double Evvel, double Pvvel, double Tau, double f, double AGy, double T) {
    return -rightEdge+PyPos+(Evvel+AGy*Tau/f)*T+Tau/f*(Pvvel-Evvel-AGy*Tau/f)*(1-exp(-T*f/Tau));
}
// double eddy::Y2(double PyPos, double Pvvel, double velEd, double Tau, double f, double AGy, double T) {
//     return -rightEdge+PyPos+velEd*(T-Tau/f*(1-exp(-T*f/Tau)));
// }

/// Helper function for eddy::getParticleUVWYafterEddy GSDB
double eddy::derY(double Evvel, double Pvvel, double Tau, double f, double AGy, double T) {
    return (Evvel+AGy*Tau/f)+(Pvvel-Evvel-AGy*Tau/f)*exp(-T*f/Tau);
}
// double eddy::derY(double Pvvel, double velEd, double Tau, double f, double AGy, double T) {
//     return velEd*(1-exp(-T*f/Tau));
// }

/// Helper function for eddy::getParticleUVWYafterEddy
double eddy::Z1(double Ewvel, double Gwvel, double Pwvel, double Tau, double f, double AGz, double T) {
    return -Ewvel*T + 0.5*eddySize+(Gwvel+AGz*Tau/f)*T+Tau/f*(Pwvel-Gwvel-AGz*Tau/f)*(1-exp(-T*f/Tau));
}
/// Helper function for eddy::getParticleUVWYafterEddy
double eddy::Z2(double Ewvel, double Gwvel, double Pwvel, double Tau, double f, double AGz, double T) {
    return -Ewvel*T - 0.5*eddySize+(Gwvel+AGz*Tau/f)*T+Tau/f*(Pwvel-Gwvel-AGz*Tau/f)*(1-exp(-T*f/Tau));
}
/// Helper function for eddy::getParticleUVWYafterEddy
double eddy::derZ(double Ewvel, double Gwvel, double Pwvel, double Tau, double f, double AGz, double T) {
    return -Ewvel + (Gwvel+AGz*Tau/f)+(Pwvel-Gwvel-AGz*Tau/f)*exp(-T*f/Tau);
}

///////////////////////////////////////////////////////////////////////////////

/** compute the index of particles that interact with eddies in type-I way when the eddy occurs
 *  iPartInEddyIC: vector stored for type-IC interaction, never update later
 *  iPartInEddy:   vector stored for type-I and -IC, update later
 *  xPartPos:      vector stored for particle position in x direction
 *  zPartPos:      vector stored for particle position in z direction
 *  @part: particle object
 *
 *  Guangyuan Sun 09/13
 */

void eddy::getIpartIxn(particles &part, double time) {

    iPartInEddyIC.clear();
    iPartInEddy.clear();
    xPartPos.clear();
    zPartPos.clear();

    for (int i=0; i<part.nPart; i++)
        if ((part.yPos[i] >= leftEdge) && (part.yPos[i] <= rightEdge)) {    // only does particles in the eddy region

            if( ((leftEdge == part.yPos[i])  && (part.vvel[i] < 0.0)) ||   // skip if part on edge and moving out of eddy box
                    ((rightEdge == part.yPos[i]) && (part.vvel[i] > 0.0))) 
                continue;
           
            iPartInEddyIC.push_back(i);
            iPartInEddy.push_back(i);
            xPartPos.push_back(0.);
            zPartPos.push_back(0.);
                       
//                        part.partIxn[i]++;
//                        if (time >= 0.2637 && time < 0.3637) part.partIxn1[i]++;
//                        if (time >= 0.3637 && time < 0.4637) part.partIxn2[i]++;
//                        if (time >= 0.4637 && time < 0.5637) part.partIxn3[i]++;
//                        if (time >= 0.5637 && time < 0.6637) part.partIxn4[i]++;
//                        if (time >= 0.6637 && time <= 0.7637) part.partIxn5[i]++;
//                        if (time >= 0.2637 && time <= 0.7637) part.partIxnFirstFiveSec[i]++;
//                        if (time >= 0 && time <= 0.7637) part.partIxnFirstFiveSecAll[i]++;
        }
}

///////////////////////////////////////////////////////////////////////////////

/** compute new veclocity and position of each particle in ODT line direction
 *  Newton's method
 *  Guangyuan Sun 09/11
 *
 *  Code computes the y-location (line direction) and u,v,w particle velocities for and
 *     instantaneous triplet map.
 *
 *  For particle interactions, eddies are assumed to have a cubical shape of size eddySize (L).
 *     Call this the eddy box.
 *  Line direction is y; x and z are streamwise and spanwise, resp.
 *  Particle initial position is L/2 for x and z, and the original y-location for y.
 *  Only particles in the box are considered.
 *  Particle initial velocity vector is as input to function.
 *
 *  Compute the interaction time of the particle in the notional eddy box.  
 *  Both the eddy and particle are assumed to move.  The eddy velocity is the average x,y,z fluid
 *     velocities in the eddy region.  Hence, we have particle motion in a moving eddy box.
 *     If the particle stays in the box, then the interaction time is the eddy time.
 *     If the particle leaves the box, then it is the exit time.
 *     Each direction is treated separately.
 *  The interaction is treated here in absolute coordinates (lab-coordinates).
 *                                                                              \cond
 *  dv_p/dt = (f/tau_p) * (v_p - v_eddy) + g      --> particle drag law
 *  dx/dt   = v_p                                 --> particle position
 *  
 *  The drag law is computed analytically to determine the interaction time.  Then, the 
 *     particle position (y) and velocity components are computed.
 *
 *  In the line direction, (y) we want only the effect of the instantaneous eddy triplet map,
 *     not the diffusional drag.  This is to avoid double counting with the subsequent 
 *     diffusion process (after eddies occur).  The y position and y-velocity are computed using the difference
 *     between the integrated drag law with and without the eddy velocity.
 *
 *  JCH added time of eddy occurrence as parameter to feed into particle histories.
 *                                                                              \endcond
 */

void eddy::tripletMapParticles(const odtParam &odtP, odtline &line, odtline &eddyLine, particles &part, double eddyStartTime, double time) {

//    for(int i=0; i<part.nPart; i++) {
//        if (!part.pActive[i]) continue;                                     // skip inactive particles (e.g. wall collisions/outflow)
//        part.vvel[i] = part.vEndLastEddy[i];
//    }

    eddyLife.resize(part.nPart, 0.0);    
    
    for(int i=0; i<part.nPart; i++) 
        eddyLife[i] = part.ParamEddylife[i]/invTauEddy;
    
    
    if(part.PeddyType == 1 || part.PeddyType == 3) // ----- TypeI or TypeIC
        getEddyVvel_TMtracers(odtP, line, eddyLine, part, eddyStartTime); // tripletMap tracer particles 
                                                       // and get eddy v velocity for inertial particles
    
    if(part.PeddyType == 2 || part.PeddyType == 3) // ----- TypeC or TypeIC
        getIpartIxn(part, time);

    if(!odtP.Iparticles || part.Ltracer) // do not do the rest if tracer particles
        return;

    getEddyUWvel(eddyLine, part);
   
    if(part.PeddyType == 1 || part.PeddyType == 3) { // ----- TypeI or TypeIC

        // ---------- ballistic particle

        if(part.Lballistic)
            tripletMapBallisticParticles(line, eddyLine, part, eddyStartTime, time);

        // ---------- inertial particle

        else 
            tripletMapInertiaParticles(line, eddyLine, part, eddyStartTime, time);
    }
} 

///////////////////////////////////////////////////////////////////////////////

void eddy::tripletMapBallisticParticles(odtline &line, odtline &eddyLine, particles &part, double eddyStartTime, double time) {
           
    double Txmin;         // interaction time in x, y, z, directions.
    double Tymin;
    double Tzmin;
    double Tmin;          // minimum interaction time = interaction time of particle with eddy.

    for(int i=0; i<part.nPart; i++) {

        if (!part.pActive[i]) continue;                                     // skip inactive particles (e.g. wall collisions/outflow)

        if (LperiodicEddy) {
            if((part.yPos[i] >= line.posf[0]) && (part.yPos[i] <= (rightEdge-line.Ldomain))) {
                part.yPos[i] += line.Ldomain;
		part.crossBound[i]--;
            }
        }

        if ((part.yPos[i] >= leftEdge) && (part.yPos[i] <= rightEdge)) {    // only does particles in the eddy region

            if ( ((leftEdge == part.yPos[i])  && (part.vvel[i] < 0.0)) ||   // skip if part on edge and moving out of eddy box
                    ((rightEdge == part.yPos[i]) && (part.vvel[i] > 0.0))) 
                continue;

            if((part.uvel[i]-uEddy[i]) < 0.0) 
                Txmin = 0.5*eddySize/(uEddy[i]-part.uvel[i]);
            else if((part.uvel[i]-uEddy[i]) > 0.0) 
                Txmin = 0.5*eddySize/(part.uvel[i]-uEddy[i]);
            else
                Txmin = eddyLife[i];

            if((part.vvel[i]-vEddy[i]) < 0.0) 
                Tymin = (part.yPos[i]-leftEdge)/(vEddy[i]-part.vvel[i]);
            else if((part.vvel[i]-vEddy[i]) > 0.0) 
                Tymin = (rightEdge-part.yPos[i])/(part.vvel[i]-vEddy[i]);
            else
                Tymin = eddyLife[i];

            if((part.wvel[i]-wEddy[i]) < 0.0) 
                Tzmin = 0.5*eddySize/(wEddy[i]-part.wvel[i]);
            else if((part.wvel[i]-wEddy[i]) > 0.0) 
                Tzmin = 0.5*eddySize/(part.wvel[i]-wEddy[i]);
            else
                Tzmin = eddyLife[i];

            double yPosOld = part.yPos[i];
            Tmin = min(min(Txmin,Tymin),min(Tymin,Tzmin));
// if(Tmin>(time-eddyStartTime)) Tmin = time-eddyStartTime;           
// if(Tmin>time) Tmin = time;           
            part.yPos[i] = part.yPos[i] + part.vvel[i]*Tmin;
            if(part.yPos[i] >= line.Ldomain) {
                part.yPos[i] -= line.Ldomain;
                part.crossBound[i]++;
            }

            //------------------- Store the effect of every eddy on each particle including 
            //                    (1)eddy end time (eddyEndTime); 
            //                    (2)relative velocity between particle and eddy (relativeVelEdPart); 

            double relativeVelEdPart = (part.yPos[i]-part.yPosTracer_TM[i])/Tmin;
            part.yPos[i] = yPosOld;
            double eddyEndTime = eddyStartTime + Tmin;

            part.updateEddyInfoArray(i, eddyEndTime, relativeVelEdPart);
        }    
    }
}

///////////////////////////////////////////////////////////////////////////////

void eddy::tripletMapInertiaParticles(odtline &line, odtline &eddyLine, particles &part, double eddyStartTime, double time) {

    double Txmin = 0.0;         // interaction time in x, y, z, directions.
    double Tymin = 0.0;
    double Tzmin = 0.0;
    double Tmin  = 0.0;          // minimum interaction time = interaction time of particle with eddy.

//     double partMomSrcInEddyReal[3] = {0,0,0}; 
    
    iPartInEddy.clear();

    for(int i=0; i<part.nPart; i++) {

        if (!part.pActive[i]) continue;                                     // skip inactive particles (e.g. wall collisions/outflow)
    
        part.set_f(i, uEddy[i], vEddy[i], wEddy[i]);
        part.set_TauP(i);

        if (LperiodicEddy) {
            if((part.yPos[i] >= line.posf[0]) && (part.yPos[i] <= (rightEdge-line.Ldomain))) {
                part.yPos[i] += line.Ldomain;
		        part.crossBound[i]--;
            }
        }

        if ((part.yPos[i] >= leftEdge) && (part.yPos[i] <= rightEdge)) {    // only does particles in the eddy region

            if( ((leftEdge == part.yPos[i])  && (part.vvel[i] < 0.0)) ||   // skip if part on edge and moving out of eddy box
                ((rightEdge == part.yPos[i]) && (part.vvel[i] > 0.0)) || 
                (vEddy[i] == 0.0)) 
                continue;

//            if(part.PeddyType == 3) { // TypeIC needs to know which particles in the eddy (diffuser::updateActiveEddy/checkActiveEddyeffect)
//                double size = iPartInEddy.size();
//                iPartInEddy.resize(size+1);
//                iPartInEddy[size] = i;
//            }

            double iyPosInEddy = eddyLine.linePositionToIndex(part.yPos[i],true);

            double Tcrit = -1;
            double dmb;
            double pPos_relEddyCenter;

            //=====================================================================================================
            //-------------------  Compute Txmin = minimum of (eddy life) and (time to cross x-eddy box boundary) (z is similar, not y)
            //=====================================================================================================

            Tcrit = -1;    // init to neg, then if a critical point exists, resets to positive; use +/- as test
            dmb   = -part.TauP[i]/part.f[i]*(part.uvel[i]-eddyLine.uvel[iyPosInEddy]-part.TauP[i]/part.f[i]*part.AGx)/
                     (eddyLine.uvel[iyPosInEddy] + part.TauP[i]/part.f[i]*part.AGx);
            if(dmb > 0.0)
                Tcrit = part.TauP[i]/part.f[i]*log(dmb);

            //-------------------

            if(Tcrit < 0.0 || Tcrit >= eddyLife[i]) {       // no critical point
                //if(true) { //doldb

                pPos_relEddyCenter = pPosRelToEdge(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, 0.0, 0.0, eddyLife[i]);

                if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) )     // in eddy box
                    Txmin = eddyLife[i];
                else if (pPos_relEddyCenter < -eddySize*0.5)   {                                         // outside on left
                    Txmin = getEddyCrossingTime(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, 0.0, -0.5*eddySize, 0.0, eddyLife[i]);
                }
                else {                                                                                  // outside on right
                    Txmin = getEddyCrossingTime(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, 0.0,  0.5*eddySize, 0.0, eddyLife[i]);
                }

            }
            //-------------------

            else {                  // has a critical point
                cout << endl << "made it x" << endl; //doldb

                pPos_relEddyCenter = pPosRelToEdge(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, 0.0, 0.0, Tcrit);

                if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) ) {   // in eddy box
                    pPos_relEddyCenter = pPosRelToEdge(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, 0.0, 0.0, eddyLife[i]);

                    if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) ) // in eddy box
                        Txmin = eddyLife[i];
                    else if (pPos_relEddyCenter < -eddySize*0.5) {                                      // outside on left
                        Txmin = getEddyCrossingTime(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, 0.0, -0.5*eddySize, 0.0, eddyLife[i]);
                    }
                    else {                                                                              // outside on right
                        Txmin = getEddyCrossingTime(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, 0.0,  0.5*eddySize, 0.0, eddyLife[i]);
                    }
                }
                else if (pPos_relEddyCenter < -eddySize*0.5)  {                                         // outside on left
                    Txmin = getEddyCrossingTime(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, 0.0, -0.5*eddySize, 0.0, Tcrit);
                }
                else {                                                                                 // outside on right
                    Txmin = getEddyCrossingTime(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, 0.0,  0.5*eddySize, 0.0, Tcrit);
                }
            

            }
            //Txmin = eddyLife; //doldb

            //=====================================================================================================
            //-------------------  Compute Tymin = minimum of (eddy life) and (time to cross y-eddy box boundary) 
            //=====================================================================================================

            Tcrit = -1;    // init to neg, then if a critical point exists, resets to positive; use +/- as test
            dmb   = -part.TauP[i]/part.f[i]*(part.vvel[i]-vEddy[i]-part.TauP[i]/part.f[i]*part.AGy)/
                     (vEddy[i] + part.TauP[i]/part.f[i]*part.AGy);
            if(dmb > 0.0)
                Tcrit = part.TauP[i]/part.f[i]*log(dmb);

            //-------------------

            pPos_relEddyCenter;

            if(Tcrit < 0.0 || Tcrit >= eddyLife[i]) {       // no critical point
                //if(true) { //doldb

                pPos_relEddyCenter = pPosRelToEdge(0.0, vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, part.yPos[i], 0.5*(leftEdge+rightEdge), eddyLife[i]);

                if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) )     // in eddy box
                    Tymin = eddyLife[i];
                else if (pPos_relEddyCenter < -eddySize*0.5) {                                           // outside on left
                    Tymin = getEddyCrossingTime(0.0, vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, part.yPos[i], leftEdge, 0.0, eddyLife[i]);
                }
                else  {                                                                                 // outside on right
                    Tymin = getEddyCrossingTime(0.0, vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, part.yPos[i],  rightEdge, 0.0, eddyLife[i]);
                }

            }
            //-------------------

            else {                  // has a critical point 
                cout << endl << "made it y" << endl; //doldb

                pPos_relEddyCenter = pPosRelToEdge(0.0, vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, part.yPos[i], 0.5*(leftEdge+rightEdge), Tcrit);

                if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) ) {   // in eddy box
                    pPos_relEddyCenter = pPosRelToEdge(0.0, vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, part.yPos[i], 0.5*(leftEdge+rightEdge), eddyLife[i]);

                    if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) ) // in eddy box
                        Tymin = eddyLife[i];
                    else if (pPos_relEddyCenter < -eddySize*0.5) {                                       // outside on left
                        Tymin = getEddyCrossingTime(0.0, vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, part.yPos[i], leftEdge,  0.0, eddyLife[i]);
                    }
                    else   {                                                                           // outside on right
                        Tymin = getEddyCrossingTime(0.0, vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, part.yPos[i], rightEdge, 0.0, eddyLife[i]);
                    }
                }
                else if (pPos_relEddyCenter < -eddySize*0.5) {                                          // outside on left
                    Tymin = getEddyCrossingTime(0.0, vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, part.yPos[i], leftEdge,  0.0, Tcrit);
                }
                else {                                                                                 // outside on right
                    Tymin = getEddyCrossingTime(0.0, vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, part.yPos[i], rightEdge, 0.0, Tcrit);
                }
            

            }
            //Tymin = eddyLife; //doldb

 //           double yPart = leftEdge + Y1(part.yPos[i], vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, eddyLife);/*{{{*//*{{{*/
 //           if ((yPart <= rightEdge) && (yPart >= leftEdge))
 //               Tymin = eddyLife;
 //           else if (yPart < leftEdge) {
 //               double Told = eddyLife;
 //               double Tnew = eddyLife;
 //               for(int j=1; j<=maxIt; j++) {
 //                   double yy1 = Y1(part.yPos[i], vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, Told);
 //                   double yy2 = derY(vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, Told);
 //                   Tnew = Told - yy1/yy2;
 //                   if (abs(Told-Tnew)<=1e-8 && Tnew > 0.) {
 //                       Tymin = Tnew;
 //                       break;
 //                   }    
 //                   Told = Tnew;
 //                   if(j==maxIt)
 //                       *proc.ostrm << endl << "# warning in eddy::getParticleUVWYafterEddy, not converged y left" << endl;
 //               }
 //           }
 //           else {
 //               double Told = eddyLife;
 //               double Tnew = eddyLife;
 //               for(int j=1; j<=maxIt; j++) {
 //                   double yy1 = Y2(part.yPos[i], vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, Told);
 //                   double yy2 = derY(vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, Told);
 //                   Tnew = Told - yy1/yy2;
 //                   if (abs(Told-Tnew)<=1e-8 && Tnew > 0.) {
 //                       Tymin = Tnew;
 //                       break;
 //                   }    
 //                   Told = Tnew;
 //                   if(j==maxIt)
 //                       *proc.ostrm << endl << "# warning in eddy::getParticleUVWYafterEddy, not converged y right" << endl;
 //               }
 //           }/*}}}*//*}}}*/
            
            //=====================================================================================================
            //-------------------  Compute Tzmin = minimum of (eddy life) and (time to cross z-eddy box boundary) (x is similar, not y)
            //=====================================================================================================

            Tcrit = -1;    // init to neg, then if a critical point exists, resets to positive; use +/- as test
            dmb   = -part.TauP[i]/part.f[i]*(part.wvel[i]-eddyLine.wvel[iyPosInEddy]-part.TauP[i]/part.f[i]*part.AGz)/
                     (eddyLine.wvel[iyPosInEddy] + part.TauP[i]/part.f[i]*part.AGz);
            if(dmb > 0.0)
                Tcrit = part.TauP[i]/part.f[i]*log(dmb);

            //-------------------

            pPos_relEddyCenter;

            if(Tcrit < 0.0 || Tcrit >= eddyLife[i]) {       // no critical point
                //if(true) { //doldb

                pPos_relEddyCenter = pPosRelToEdge(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, 0.0, 0.0, eddyLife[i]);

                if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) )     // in eddy box
                    Tzmin = eddyLife[i];
                else if (pPos_relEddyCenter < -eddySize*0.5) {                                          // outside on left
                    Tzmin = getEddyCrossingTime(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, 0.0, -0.5*eddySize, 0.0, eddyLife[i]);
                }
                else {                                                                                 // outside on right
                    Tzmin = getEddyCrossingTime(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, 0.0,  0.5*eddySize, 0.0, eddyLife[i]);
                }

            }
            //-------------------

            else {                  // has a critical point
                cout << endl << "made it z" << endl; //doldb

                pPos_relEddyCenter = pPosRelToEdge(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, 0.0, 0.0, Tcrit);

                if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) ) {   // in eddy box
                    pPos_relEddyCenter = pPosRelToEdge(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, 0.0, 0.0, eddyLife[i]);

                    if( (pPos_relEddyCenter >= -eddySize*0.5 ) && (pPos_relEddyCenter <= eddySize*0.5) ) // in eddy box
                        Tzmin = eddyLife[i];
                    else if (pPos_relEddyCenter < -eddySize*0.5)  {                                     // outside on left
                        Tzmin = getEddyCrossingTime(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, 0.0, -0.5*eddySize, 0.0, eddyLife[i]);
                    }
                    else {                                                                             // outside on right
                        Tzmin = getEddyCrossingTime(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, 0.0,  0.5*eddySize, 0.0, eddyLife[i]);
                    }
                }
                else if (pPos_relEddyCenter < -eddySize*0.5) {                                          // outside on left
                    Tzmin = getEddyCrossingTime(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, 0.0, -0.5*eddySize, 0.0, Tcrit);
                }
                else  {                                                                                // outside on right
                    Tzmin = getEddyCrossingTime(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, 0.0,  0.5*eddySize, 0.0, Tcrit);
                }
            

            }
            //Tzmin = eddyLife[i]; //doldb

            //double zEddy = wEddy[i] * eddyLife[i];/*{{{*/
            //double zPart = zEddy-0.5*eddySize + 
            //    Z1(wEddy[i], eddyLine.wvel[iyPosInEddy], 
            //            part.wvel[i], part.TauP[i], 
            //            part.f[i], part.AGz, eddyLife[i]);

            //if ((zPart <= (0.5*eddySize+zEddy)) && 
            //        (zPart >= (-0.5*eddySize+zEddy)))
            //    Tzmin = eddyLife;
            //else if (zPart < (-0.5*eddySize+zEddy)) {

            //    // ----- GSDB 06272013---- Find better staring point for Newton's method
            //                    
            //    double init1 = Z1(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, 0);
            //    double init2 = Z1(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, eddyLife/2);
            //    double init3 = Z1(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, eddyLife);
            //    
            //    double Told = 0.;
            //    double Tnew = 0.;
            //    if(init1*init2 < 0.) Told = eddyLife/4;
            //    else if(init2*init3 < 0.) Told = eddyLife*3/4;
            //    else Told = eddyLife;
            //    
            //    //                 double Told = eddyLife;
            //    //                 double Tnew = eddyLife;
            //    // ----- GSDB 06272013---- Find better staring point for Newton's method
            //    
            //    for(int j=1; j<=maxIt; j++) {
            //        double zz1 = Z1(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, Told);
            //        double zz2 = derZ(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, Told);
            //        Tnew = Told - zz1/zz2;
            //        if (abs(Told-Tnew)<=1e-8 && Tnew > 0.) {
            //            Tzmin = Tnew;
            //            break; 
            //        }
            //        Told = Tnew;
            //        if(j==maxIt)
            //            *proc.ostrm << endl << "# warning in eddy::getParticleUVWYafterEddy, not converged z left" << endl;
            //    }
            //}
            //else {

            //    // ----- GSDB 06272013---- Find better staring point for Newton's method
            //                    
            //    double init1 = Z2(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, 0);
            //    double init2 = Z2(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, eddyLife/2);
            //    double init3 = Z2(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, eddyLife);
            //    
            //    double Told = 0.;
            //    double Tnew = 0.;
            //    if(init1*init2 < 0.) Told = eddyLife/4;
            //    else if(init2*init3 < 0.) Told = eddyLife*3/4;
            //    else Told = eddyLife;
            //    
            //    //                 double Told = eddyLife;
            //    //                 double Tnew = eddyLife;
            //    // ----- GSDB 06272013---- Find better staring point for Newton's method
            //    
            //    for(int j=1; j<=maxIt; j++) {
            //        double zz1 = Z2(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, Told);
            //        double zz2 = derZ(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, Told);
            //        Tnew = Told - zz1/zz2;
            //        if (abs(Told-Tnew)<=1e-8 && Tnew > 0.) {
            //            Tzmin = Tnew;
            //            break;
            //        }
            //        Told = Tnew;
            //        if(j==maxIt)
            //            *proc.ostrm << endl << "# warning in eddy::getParticleUVWYafterEddy, not converged z right" << endl;
            //    }
            //}
            ////Tzmin = eddyLife; //doldb/*}}}*/

            //=====================================================================================================
            //=====================================================================================================

            //-------------------  

            Tmin = min(min(Txmin,Tymin),min(Tymin,Tzmin));
// if(Tmin>(time-eddyStartTime)) Tmin = time-eddyStartTime;           
// if(Tmin>time) Tmin = time;           
//                       part.partCross[i]++;
//                       if (Tmin != eddyLife[i]) {
//                       if (time >= 0.2637 && time < 0.3637) part.partCross1[i]++;
//                       if (time >= 0.3637 && time < 0.4637) part.partCross2[i]++;
//                       if (time >= 0.4637 && time < 0.5637) part.partCross3[i]++;
//                       if (time >= 0.5637 && time < 0.6637) part.partCross4[i]++;
//                       if (time >= 0.6637 && time <= 0.7637) part.partCross5[i]++;
//                       if (time >= 0.2637 && time <= 0.7637) part.partCrossFirstFiveSec[i]++;
//                       if (time >= 0 && time <= 0.7637) part.partCrossFirstFiveSecAll[i]++;
//                       }
//                       part.partIxn[i]++;
//                       if (time >= 0.2637 && time < 0.3637) part.partIxn1[i]++;
//                       if (time >= 0.3637 && time < 0.4637) part.partIxn2[i]++;
//                       if (time >= 0.4637 && time < 0.5637) part.partIxn3[i]++;
//                       if (time >= 0.5637 && time < 0.6637) part.partIxn4[i]++;
//                       if (time >= 0.6637 && time <= 0.7637) part.partIxn5[i]++;
//                       if (time >= 0.2637 && time <= 0.7637) part.partIxnFirstFiveSec[i]++;
//                       if (time >= 0 && time <= 0.7637) part.partIxnFirstFiveSecAll[i]++;

            double cons1 = part.TauP[i]/part.f[i];
            double cons2 = uEddy[i]+cons1*part.AGx;
            double cons3 = cons1*(1-exp(-Tmin/cons1));
            double cons4 = wEddy[i]+cons1*part.AGz;
            double cons5 = vEddy[i]+cons1*part.AGy;

            //------------------- adjust particle histories

            if ( part.Lhistories ) {
                part.adjustEddyVelHistories( i, eddyStartTime, Tmin, cons1, vEddy[i], part.vvel[i]);
            }

            //-------------------
            
            if(part.PeddyType == 1 || part.PeddyType == 3) { // ---- TypeI and Type IC

                // part.timeTotal[i] = part.timeTotal[i] + Tmin;  /// for test
                
                part.yPos[i] = part.yPos[i] + vEddy[i]*(Tmin-cons3);
                if(part.yPos[i] >= line.Ldomain) {    // dol todo: fix this
                    part.yPos[i] -= line.Ldomain;
                    part.crossBound[i]++;
                }
                part.vvel[i] = part.vvel[i]+vEddy[i]*(1-exp(-Tmin/cons1));

                //part.vEndLastEddy[i] = part.vvel[i]; 
                //part.vEndLastEddy[i] = vEddy[i] + part.TauP[i]/part.f[i]*part.AGy - 
                //                       (part.TauP[i]/part.f[i]*part.AGy + vEddy[i] - part.vEndLastEddy[i])*
                //                       exp(-Tmin/part.TauP[i]*part.f[i]);

                //------------------- Store the effect of every eddy on each particle including 
                //                    (1)eddy end time (eddyEndTime); 
                //                    (2)relative velocity between particle and eddy (relativeVelEdPart); 

                double relativeVelEdPart = (part.yPos[i]-part.yPosTracer_TM[i])/Tmin;
                double eddyEndTime = eddyStartTime + Tmin;
                part.updateEddyInfoArray(i, eddyEndTime, relativeVelEdPart);

            }
            
            else if(part.PeddyType == 2) {          // ----- TypeC
                
                double tPos = part.yPos[i];
                part.yPos[i] = part.yPos[i] + vEddy[i]*(Tmin-cons3);
                if(part.yPos[i] >= line.Ldomain) {    // dol todo: fix this
                    part.yPos[i] -= line.Ldomain;
                    part.crossBound[i]++;
                }

                //------------------- Store the effect of every eddy on each particle including 
                //                    (1)eddy end time (eddyEndTime); 
                //                    (2)relative velocity between particle and eddy (relativeVelEdPart); 

                double relativeVelEdPart = (part.yPos[i]-part.yPosTracer_TM[i])/Tmin;
                double eddyEndTime = eddyStartTime + Tmin;
                part.updateEddyInfoArray(i, eddyEndTime, relativeVelEdPart);

                part.yPos[i] = tPos; 
            }
        
        } 
        part.set_f(i, uEddy[i], vEddy[i], wEddy[i]);
        part.set_TauP(i);
    } // inertial particle loop end
}

///////////////////////////////////////////////////////////////////////////////

bool eddy::computePartSrcInEddy(const odtParam &odtP, odtline &line, odtline &eddyLine, particles &part, double eddyStartTime, double time, double Z_value) { 

    int maxIt=1000;
    
    getEddyUWvel(eddyLine, part);
    
    double Txmin;         // interaction time in x, y, z, directions.
    double Tymin;
    double Tzmin;
    double Tmin;          // minimum interaction time = interaction time of particle with eddy.
    
    for(int k=0; k<maxIt; k++) {  // calculate invTau and get converged
    
    eddyLife.resize(part.nPart, 0.0);    
    
    for(int i=0; i<part.nPart; i++) 
        eddyLife[i] = part.ParamEddylife[i]/invTauEddy;
    
    
    getEddyVvel_TMtracers(odtP, line, eddyLine, part, eddyStartTime); // tripletMap tracer particles 
                                                       // and get eddy v velocity for inertial particles
    partMomSrcInEddy[0] = 0.; 
    partMomSrcInEddy[1] = 0.;
    partMomSrcInEddy[2] = 0.;
    partErgSrcInEddy[0] = 0.; 
    partErgSrcInEddy[1] = 0.;
    partErgSrcInEddy[2] = 0.;
    
    Txmin = 0.0;         // interaction time in x, y, z, directions.
    Tymin = 0.0;
    Tzmin = 0.0;
    Tmin = 0.0;          // minimum interaction time = interaction time of particle with eddy.

    for(int i=0; i<part.nPart; i++) {
        if (!part.pActive[i]) continue;                                     // skip inactive particles (e.g. wall collisions/outflow)
        
        part.set_f(i, uEddy[i], vEddy[i], wEddy[i]);
        part.set_TauP(i);

        if (LperiodicEddy) {
            if((part.yPos[i] >= line.posf[0]) && (part.yPos[i] <= (rightEdge-line.Ldomain))) {
                part.yPos[i] += line.Ldomain;
            }
        }

        if ((part.yPos[i] >= leftEdge) && (part.yPos[i] <= rightEdge)) {    // only does particles in the eddy region

            if ( ((leftEdge == part.yPos[i])  && (part.vvel[i] < 0.0)) ||   // skip if part on edge and moving out of eddy box
                    ((rightEdge == part.yPos[i]) && (part.vvel[i] > 0.0)) || 
                    (vEddy[i] == 0.0)) {
                continue;
            }

            double iyPosInEddy = eddyLine.linePositionToIndex(part.yPos[i],true);

            //-------------------  Compute Txmin = minimum of (eddy life) and (time to cross x-eddy box boundary) (z is similar, not y)
            
            double xEddy = uEddy[i] * eddyLife[i];           // displacement of eddy box over time eddyLife
            double xPart = xEddy-0.5*eddySize +        // displacement of particle over time eddyLife
                X1(uEddy[i], eddyLine.uvel[iyPosInEddy], 
                        part.uvel[i], part.TauP[i], 
                        part.f[i], part.AGx, eddyLife[i]);

            if ((xPart <= (0.5*eddySize+xEddy)) &&     // check if particle is in the eddy box,
                    (xPart >= (-0.5*eddySize+xEddy)))      //       if it is, the interaction time is the eddy life
                Txmin = eddyLife[i];                   
            else if (xPart < (-0.5*eddySize+xEddy)) {  // else particle crosses on left edge
                double Told = eddyLife[i];
                double Tnew = eddyLife[i];
                for(int j=1; j<=maxIt; j++) {
                    double xx1 = X1(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, Told);
                    double xx2 = derX(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, Told);
                    Tnew = Told - xx1/xx2;
                    if (abs(Told-Tnew)<=1e-8 && Tnew > 0.) { 
                        Txmin = Tnew;
                        break; 
                    }
                    Told = Tnew;
                    if(j==maxIt) 
                        return false;
//                        *proc.ostrm << endl << "# warning in eddy::getParticleUVWYafterEddy, not converged x left in source function" << endl;
                }
            }
            else {                                    // else particle crosses on right edge
                double Told = eddyLife[i];
                double Tnew = eddyLife[i];
                for(int j=1; j<=maxIt; j++) {
                    double xx1 = X2(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, Told);
                    double xx2 = derX(uEddy[i], eddyLine.uvel[iyPosInEddy], part.uvel[i], part.TauP[i], part.f[i], part.AGx, Told);
                    Tnew = Told - xx1/xx2;
                    if (abs(Told-Tnew)<=1e-8 && Tnew > 0.) {
                        Txmin = Tnew;
                        break;
                    }                            
                    Told = Tnew;
                    if(j==maxIt) 
                        return false;
//                        *proc.ostrm << endl << "# warning in eddy::getParticleUVWYafterEddy, not converged x right in source function" << endl;
                }
            }

            //-------------------  Compute Tymin = minimum of (eddy life) and (time to cross y-eddy box boundary) 

            double yPart = leftEdge + Y1(part.yPos[i], vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, eddyLife[i]);
            if ((yPart <= rightEdge) && (yPart >= leftEdge))
                Tymin = eddyLife[i];
            else if (yPart < leftEdge) {
                double Told = eddyLife[i];
                double Tnew = eddyLife[i];
                for(int j=1; j<=maxIt; j++) {
                    double yy1 = Y1(part.yPos[i], vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, Told);
                    double yy2 = derY(vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, Told);
                    Tnew = Told - yy1/yy2;
                    if (abs(Told-Tnew)<=1e-8 && Tnew > 0.) {
                        Tymin = Tnew;
                        break;
                    }    
                    Told = Tnew;
                    if(j==maxIt) 
                        return false;
//                        *proc.ostrm << endl << "# warning in eddy::getParticleUVWYafterEddy, not converged y left in source function" << endl;
                }
            }
            else {
                double Told = eddyLife[i];
                double Tnew = eddyLife[i];
                for(int j=1; j<=maxIt; j++) {
                    double yy1 = Y2(part.yPos[i], vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, Told);
                    double yy2 = derY(vEddy[i], part.vvel[i], part.TauP[i], part.f[i], part.AGy, Told);
                    Tnew = Told - yy1/yy2;
                    if (abs(Told-Tnew)<=1e-8 && Tnew > 0.) {
                        Tymin = Tnew;
                        break;
                    }    
                    Told = Tnew;
                    if(j==maxIt) 
                        return false;
//                        *proc.ostrm << endl << "# warning in eddy::getParticleUVWYafterEddy, not converged y right source function" << endl;
                }
            }

            //-------------------  Compute Tzmin = minimum of (eddy life) and (time to cross z-eddy box boundary) (x is similar, not y)

            double zEddy = wEddy[i] * eddyLife[i];
            double zPart = zEddy-0.5*eddySize + Z1(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, eddyLife[i]);

            if ((zPart <= (0.5*eddySize+zEddy)) && (zPart >= (-0.5*eddySize+zEddy)))
                Tzmin = eddyLife[i];
            else if (zPart < (-0.5*eddySize+zEddy)) {
                double Told = eddyLife[i];
                double Tnew = eddyLife[i];
                for(int j=1; j<=maxIt; j++) {
                    double zz1 = Z1(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, Told);
                    double zz2 = derZ(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, Told);
                    Tnew = Told - zz1/zz2;
                    if (abs(Told-Tnew)<=1e-8 && Tnew > 0.) {
                        Tzmin = Tnew;
                        break; 
                    }
                    Told = Tnew;
                    if(j==maxIt) 
                        return false;
//                        *proc.ostrm << endl << "# warning in eddy::getParticleUVWYafterEddy, not converged z left in source function" << endl;
                }
            }
            else {
                double Told = eddyLife[i];
                double Tnew = eddyLife[i];
                for(int j=1; j<=maxIt; j++) {
                    double zz1 = Z2(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, Told);
                    double zz2 = derZ(wEddy[i], eddyLine.wvel[iyPosInEddy], part.wvel[i], part.TauP[i], part.f[i], part.AGz, Told);
                    Tnew = Told - zz1/zz2;
                    if (abs(Told-Tnew)<=1e-8 && Tnew > 0.) {
                        Tzmin = Tnew;
                        break;
                    }
                    Told = Tnew;
                    if(j==maxIt) 
                        return false;
//                        *proc.ostrm << endl << "# warning in eddy::getParticleUVWYafterEddy, not converged z rightin source function" << endl;
                }
            }

            //-------------------  

            Tmin = min(min(Txmin,Tymin),min(Tymin,Tzmin));
// if(Tmin>(time-eddyStartTime)) Tmin = time-eddyStartTime;           
// if(Tmin>time) Tmin = time;           
            
            double cons1 = part.TauP[i]/part.f[i];
            double cons2 = uEddy[i]+cons1*part.AGx;
            double cons3 = cons1*(1-exp(-Tmin/cons1));
            double cons4 = wEddy[i]+cons1*part.AGz;
            double cons5 = vEddy[i]+cons1*part.AGy;

            //-------------------

            double uvel = part.uvel[i];
            double vvel = part.vvel[i]+vEddy[i]*(1-exp(-Tmin/cons1));
            double wvel = part.wvel[i];
            double pMass = part.pDens0[i]*4./3.*3.14159*part.pRadi[i]*part.pRadi[i]*part.pRadi[i];
            partMomSrcInEddy[0] += pMass*part.nInPseudoPart[i]*(uvel-part.uvel[i]);
            partMomSrcInEddy[1] += pMass*part.nInPseudoPart[i]*(vvel-part.vvel[i]);
            partMomSrcInEddy[2] += pMass*part.nInPseudoPart[i]*(wvel-part.wvel[i]);
            partErgSrcInEddy[0] += 0.5*pMass*part.nInPseudoPart[i]*(uvel*uvel-part.uvel[i]*part.uvel[i]);
            partErgSrcInEddy[1] += 0.5*pMass*part.nInPseudoPart[i]*(vvel*vvel-part.vvel[i]*part.vvel[i]);
            partErgSrcInEddy[2] += 0.5*pMass*part.nInPseudoPart[i]*(wvel*wvel-part.wvel[i]*part.wvel[i]);
        } 
    }
    
    double invTauEddy0 = invTauEddy;
 
    if(odtP.LconstProp)  {
        if(!eddyTauCPpartSrc(eddyLine, odtP, Z_value))
            return false;
    }
    else {
        if(!eddyTauPartSrc(eddyLine, odtP, Z_value))
            return false;
    }

    if (abs(invTauEddy-invTauEddy0)/invTauEddy0<=1e-2 && invTauEddy > 0.) 
        {
        invTauEddy = invTauEddy0;
        break;
        }    
    
    if(k==maxIt) {
        *proc.ostrm << endl << "# warning in eddy::compute, not converged in two-way function" << endl;
        return false;   
    }
}
    return true;
}

/////////////////////////////DOXYGEN DOCUMENTATION//////////////////////////////////////////////////

/*! \fn void eddy::tripletMapParticles(const odtParam &odtP, odtline &line, odtline &eddyLine, particles &part, double eddyStartTime)
 *  \n
 *  \fun{\frac{dv_p}{dt} = \frac{f}{\tau_p} (v_p - v_{eddy}) + g}      --> particle drag law\n
 *  \fun{\frac{dx}{dt}   = v_p}                                        --> particle position\n\n
 *
 *  The drag law is computed analytically to determine the interaction time.  Then, the
 *     particle position (y) and velocity components are computed.
 *
 *  In the line direction, (y) we want only the effect of the instantaneous eddy triplet map,
 *     not the diffusional drag.  This is to avoid double counting with the subsequent
 *     diffusion process (after eddies occur).  The y position and y-velocity are computed using the difference
 *     between the integrated drag law with and without the eddy velocity.
 *
 *  JCH added time of eddy occurrence as parameter to feed into particle histories.
 *
 * @param odtP \input parameters object
 * @param line \input current odtline
 * @param eddyLine \input odtline object that contains info about the eddy line.
 * @param part \inout particles object
 * @param eddyStartTime \input eddy start time
 */

/////////////////////////////END DOXYGEN DOCUMENTATION////////////////////////////////////////////////
