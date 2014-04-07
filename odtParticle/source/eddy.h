/**
 * @file eddy.h
 * Header file for class eddy
 */

#ifndef EDDY_H
#define EDDY_H

#include <iostream>
#include <vector>
#include <fstream>
#include "anyline.h"
#include "odtline.h"
#include "odtParam.h"
#include "randomGenerator.h"
#include "particles.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////

/** An eddy class implementing all things eddy.  This includes triplet maps, 
 *  eddy timescales, large eddy supression,
 *  eddy acceptance, (or all things eddy).  Much of this happens on odt line segments
 *  that are separate objects from the main lines.
 *  
 *  @author David O. Lignell
 */ 

class eddy {

  ////////////////////// DATA MEMBERS /////////////////////

    public:

        double              leftEdge;           ///< left edge location of eddy
        double              rightEdge;          ///< right edge location of eddy
        double              invTauEddy;         ///< inverse eddy timescale 
        double              eddySize;           ///< size of eddy
        double              Pa;                 ///< eddy acceptance probability
        bool                LperiodicEddy;      ///< a wrap-around eddy
        bool                LmultiPhaseEddy;    ///< 
        vector<double>      cCoef;              ///< coefficient of K kernel
        vector<double>      bCoef;              ///< coefficient of J kernel
        vector<double>      K;                  ///< eddy kernel K
        vector<double>      partMomSrcInEddy;   ///< particle momentum source term (two-way coupling)
        vector<double>      partErgSrcInEddy;   ///< particle energy source term (two-way coupling)

        vector<double>      vEddy;              ///< v eddy velocity
        vector<double>      uEddy;              ///< u eddy velocity
        vector<double>      wEddy;              ///< w eddy velocity
        double              uEddyAvg;           ///< average u eddy velocity
        double              vEddyAvg;           ///< average v eddy velocity
        double              wEddyAvg;           ///< average w eddy velocity
        vector<double>      eddyLife;
        // double              eddyLife;
        
        vector<int>         iPartInEddyIC;      ///< index of particles that have interaction with eddy in type-C, never update
        vector<int>         iPartInEddy;        ///< index of particles that have interaction with eddy, update 
        vector<double>      xPartPos;           ///< x position of particles
        vector<double>      zPartPos;           ///< z position of particles

        randomGenerator     *rnd;            ///< random generator object

#ifdef COMPSCI
        vector<double>      eddyTauDy;          ///< used by eddyTauCP
#endif

        ////////////////////// MEMBER FUNCTIONS  /////////////////////
        
        void   tripMap(anyline &anyl);
        bool   eddyTau(odtline &line, const odtParam &odtP, double Z_value);
        bool   eddyTauPartSrc(odtline &line, const odtParam &odtP, double Z_value);
        bool   eddyTauCP(odtline &line, const odtParam &odtP, double Z_value);
        bool   eddyTauCPpartSrc(odtline &line, const odtParam &odtP, double Z_value);
        void   sampleEddySize(const odtParam &odtP);                 
        void   sampleEddyPosition(const odtParam &odtP, anyline &line);
        void   computeEddyAcceptanceProb(const odtParam &odtP, double &dtSample);
        void   applyVelocityKernels(odtline &line, const odtParam &odtP);
        void   fillKernel(anyline &line); 
        double eddyFavreAvgVelocity(odtline &line);

        //----------- for particles
        
        void   getIpartIxn(particles &part, double time);
        void   tripletMapParticles(const odtParam &odtP, odtline &line, odtline &eddyLine, particles &part, double eddyStartTime, double time);
        void   tripletMapBallisticParticles(odtline &line, odtline &eddyLine, particles &part, double eddyStartTime, double time);
        void   tripletMapInertiaParticles(odtline &line, odtline &eddyLine, particles &part, double eddyStartTime, double time);
        bool   computePartSrcInEddy(const odtParam &odtP, odtline &line, odtline &eddyLine, particles &part, double eddyStartTime, double time, double Z_value); 

    private : 

        void   getEddyVvel_TMtracers(const odtParam &odtP, odtline &line, odtline &eddyLine, particles &part, double eddyStartTime);
        void   getEddyUWvel(odtline &line, particles &part);
        //---- helper functions for getParticleUVWYafterEddy
        double X1(double Euvel, double Guvel, double Puvel, double Tau, double f, double AGx, double T);
        double X2(double Euvel, double Guvel, double Puvel, double Tau, double f, double AGx, double T); 
        double derX(double Euvel, double Guvel, double Puvel, double Tau, double f, double AGx, double T);
        double Y1(double PyPos, double Pvvel, double velEd, double Tau, double f, double AGy, double T);
        double Y2(double PyPos, double Pvvel, double velEd, double Tau, double f, double AGy, double T);
        double derY(double Pvvel, double velEd, double Tau, double f, double AGy, double T);
        double Z1(double Ewvel, double Gwvel, double Pwvel, double Tau, double f, double AGz, double T);
        double Z2(double Ewvel, double Gwvel, double Pwvel, double Tau, double f, double AGz, double T);
        double derZ(double Ewvel, double Gwvel, double Pwvel, double Tau, double f, double AGz, double T);
        double pPosRelToEdge(double Evel, double Gvel, double Pvel, double tauP, double f, double AG, double pPos0,double edgePos0, double T);
        double ddt_pPosRelToEdge(double Evel, double Gvel, double Pvel, double tauP, double f, double AG, double T);
        double getEddyCrossingTime(double Evel, double Gvel, double Pvel, double tauP, double f, double AG, double pPos0,double edgePos0, double T1, double T2);

        ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public :

        eddy();                    // constructor
        eddy(const eddy &ed);      // copy constructor subset
        ~eddy() {}                 // destructor

};

#endif
