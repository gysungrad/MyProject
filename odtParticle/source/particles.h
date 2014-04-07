/**
 * @file particles.h
 * Header file for class diffuser
 */

#ifndef PARTICLES_H
#define PARTICLES_H

#include <vector>
#include <string>
#include "odtParam.h"

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

//#include "eddy.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////

class odtline;
class diffuser;
/** Lagrangian particles class
 *
 */

class particles {

  ////////////////////// DATA MEMBERS /////////////////////

    public:

        odtParam                           *odtP;                ///< pointer to odtParam object
        odtline                            *line;                ///< pointer to line object
        diffuser                           *diff;                ///< back pointer to the diffuser

        int                                nPart;                ///< number of lagrangian particles (each has nInPseudoPart "real" particles)

        double                             dtStepCFL;            ///< basic step size for particles (for comparing to diffuser step size)
               
        vector<double>                     yPos;                 ///< particle position on the line
	    vector<int>                        crossBound;           ///< track the number of times each particle crosses a boundary
        vector<int>                        iyPos;                ///< index of particles in the line
        vector<double>                     fracC;                ///< particle location in the cell
        vector<double>                     uvel;                 ///< particle velocity in streamwise dir
        vector<double>                     vvel;                 ///< particle velocity in line dir
        vector<double>                     wvel;                 ///< particle velocity in spanwise dir

        vector<bool>                       pActive;              ///< flag of whether the particle is live

        int                                pShape;               ///< particle Shape (1 for cylinder, 2 for sphere, 0 for flat plate)
        vector<double>                     pRadi;                ///< particle radius
        double                             pLength;              ///< length of the particle(for cylindrical particles)
        double                             volOfFakePart;        ///< volume of the "fake" particle (all particles same, 4/3*pi*r^3 or pi*r^2*l) (a convenience var)
                  
        double                             AGx;                  ///< acceleration in x 
        double                             AGy;                  ///< acceleration in y 
        double                             AGz;                  ///< acceleration in z 
                  
        bool                               Ltracer;              ///< true if tracer particles (noninertial, follows the fluid)
        bool                               Lballistic;           ///< true if ballastic particles(huge inertial, never follow the fluid)
        bool                               Lhistories;           ///< true if collecting particle histories

        bool                               moveEddy;             ///< true if eddy "box" is moving in type-C and type-IC interaction

        int                                initPartLoc;          ///< mode of particle initialization: 0=rand distrib across domain;1=in order across domain
        double                             voidFrac;             ///< from input file: void fraction to be set in the line
        double                             initPartLocFrac;      ///< fraction of domain over which particles are initially placed
        double                             initPartLeftLocFrac;  ///< leftmost initial location of particles in domain
        int                                nPartInitLoc;         ///< mode 3; size of vector of initial particle location
        vector<double>                     partInitLoc;          ///< mode 3; vector of initial particle location

        int                                initPartVel;          ///< mode of particle initial velocity: 0=zero velocity; 1=line velocity; 2=hard coding in particle.inp
        double                             initPartUvel;         ///< mode 2 (hard coding) initial particle U velocity
        double                             initPartVvel;         ///< mode 2 (hard coding) initial particle V velocity
        double                             initPartWvel;         ///< mode 2 (hard coding) initial particle W velocity

        vector<vector<double> >            historiesTime;        ///< time for each history time step 
        vector<vector<double> >            historiesYpos;        ///< particle y position 
        vector<vector<double> >            historiesVpar;        ///< particle v velocity -- adjusted back in time 
        vector<vector<double> >            historiesTpar;        ///< particle temperature 
        vector<vector<double> >            historiesVslip;       ///< particle v slip velocity 
        vector<vector<double> >            historiesTgas;        ///< gas temp seen by particle 
        vector<vector<double> >            historiesZmix;        ///< mix frac seen by particle 
        vector<vector<double> >            historiesChi;         ///< scalar diss rate seen by particle 
        vector<vector<double> >            historiesDiff;        ///< mix frac diffusion seen by particle 
        vector<vector<double> >            historiesGradZ;       ///< mix frac gradient (ballistic mix frac rate )
        double                             timeNextHistoryPoint; ///< time when the next history time step is written 
	    double                             deltaTimeHistory;     ///< interval between particle histories--should be a few times less than tauP and Kolmogorov time 

        vector<double>                     f;                    ///< nonlinear correction factor
        vector<double>                     TauP;                 ///< response time of particles (Stokes number)
        int                                nParamEddylife;       ///< number of ParamEddylife 
        vector<double>                     ParamEddylife;        ///< Eddylife = ParamEddylife * invTauEddy (eddy::tripletMapParticle)
        class eddyInformation {                                  ///< for heat transfer coefficient (type C interaction)
            public:
            double endTime;
            double relativeVel;
            double eddyVvel;
        };
        vector<double>                     yPosTracer_TM;         ///< ?
        vector<vector<eddyInformation> >   eddyInfo;              ///< ?
        vector<double>                     vHT;                   ///< ?

        vector<double>                         PeddyUvel;       ///<gas velocity due to sum of mutiple eddy effect in diffustion(typeC two-way coupling for particles)
        vector<double>                         PeddyVvel;
        vector<double>                         PeddyWvel;
        int                                    PeddyType;       ///< 1 ---- type I;    2 ---- typeC;    3 ---- type IC 
        int                                    Ldeposition;     ///< 1 ---- deposition statistic; 0 ---- velocity statistic
        int                                    Ndpsn1;          ///< cumulative number part dep on left edge.
        int                                    Ndpsn2;          ///< cumulative number part dep on right edge.

        // -------------- internal particle heat equation solve

        int                                radial_ngrd;          ///< number of grid point in the radial direction of the particle
        vector< vector<double> >           cellTemp;             ///< particle temperature (sub-vector is the temperature down the radius)
        vector< vector<double> >           cellMass;             ///< particle cell mass
        vector< vector<double> >           cellEnth;             ///< particle cell enthalpy
        vector<vector<double> >            face_locs;            ///< face locations (first vector for each particle, each sub-vector contains the face locations)
        vector<vector<double> >            cell_locs;            ///< cell locations (first vector for each particle, each sub-vector contains the cell locations)
        double                             cell_size;            ///< cell size
        vector<vector<double> >            cellVolume;           ///< cell volumes
        vector<vector<double> >            cellArea;             ///< cell area
        vector<vector<double> >            cellFaceHeatFlux;     ///< heat flux for particles
        vector<vector<double> >            cellDens;             ///< particle density

        vector<vector<double> >            cellMassSource;       ///< particle mass source (negative if Lprxn flag is on)
        vector<vector<double> >            cellEnthSource;       ///< particle enthalpy source
        vector<vector<double> >            cellTotalEnthSource;  ///< particle total enthalpy source (d(m.p*h.p)/dt)
        
        // -------------- added for reacting particles
  
        int                               LuniformPart;           ///< flag to uniform particle; varing radius * density * paramEddylife
        int                                npRadiInitial;          ///< number of multiple initial particle radius
        vector<double>                     pRadi0;                 ///< initial particle radii
        int                                npDensInitial;          ///< number of multiple initial particle density
        vector<double>                     pDens0;                 ///< initial particle density
        double                             initPartMass_i;         ///< Initial individual particle mass

        vector<double>                     nInPseudoPart;          ///< # of given "particle" (each particle represents a collection at the given size)

        vector<double>                     particleMassSource;     ///< particle mass Source(needed to compute blowing factor)

        double                             emiss;                  ///< particle emissivities
        double                             Hf;                     ///< particle heat of formation (J/kg)
        double                             pCp;                    ///< particle heat capacity     (J/kg*K)

        int                                iPtMass;                ///< Pointer to Mass in gSource vector.
        int                                iPtEnth;                ///< Pointer to Enthalpy in gSource vector (when there is no Ieta)
        int                                iPtUvel;                ///< Pointer to U velocity in gSource vector.
        int                                iPtVvel;                ///< Pointer to V velocity in gSource vector.
        int                                iPtWvel;                ///< Pointer to W velocity in gSource vector.
        int                                iPtEta;                 ///< Pointer to Eta in gSource vector.
        int                                iPtMom;                 ///< Pointer to Mom in gSource vector.

        double                             totalVolatileMass;      ///< keep track of total volatile released to the gas
        double                             totalVolatileEnthalpy;  ///< keep track of total volatile enthalpy released to the gas
        double                             peakTemp;               ///< keep track of peak particle temperature
        double                             totalHeatDueToConvec;   ///< keep track of contribution of total heat to particle from convection
        double                             totalHeatDueToRad;      ///< keep track of contribution of total heat to particle from radiation

        bool                               LuseDiBlasiModel;       ///< flag to use Di Blasis particle combustion model
        bool                               LuseNunnsModel;         ///< flag to use Nunns particle combustion model
        
        bool                               LInitProfile;           ///< flag to use inital burnt profile for gas and particles
        string                             initProfileFile;        ///< variable that stores the name of the initial burnt profile file.
        std::vector<double>                initX;
        std::vector<double>                initU;

        //-------------------- added for particle moisture
        std::vector<vector<double> >       cellMoistureSource;    ///< particle moisture source
        std::vector<vector<double> >       moisture_tracker;      ///< particle moisture source tracker (when certain % of the mass is vaproized, cellMoistureSource=0;
        
        std::vector<vector<double> >       rho_M;                 ///< density of moisture

        double                             moist_content;         ///< particle moisture content
        double                             heatOfVaporization;    ///< heat of vaporization
     
        ////////////////////// MEMBER FUNCTIONS  /////////////////////
    
        void outputProperties(string fname);
        void readProperties(string fname);
        void readInitialProfile();
        void setParticles(string fname);
        void diffuseTracerParticle(vector<double> &Gvel); 
        void computeRHSFAndSetGasSource(vector<double> &Gvel,
                                        vector<vector<double> > &Prhs, 
                                        double dtStep,
                                        vector<double> & lambda_gas,
                                        vector<double>& dxML,
                                        vector<vector<double> >* gSource,
                                        double time);
        void set_f(int iPart, double Ug, double Vg, double Wg);
        void set_TauP(int iPart); 
        void checkBoundsSetInactive(vector<double> &Gvel, double time=-1);
        void setFracC();
        void updateEddyInfoArray(int iPart, double eddyEndTime, double relativeVelEdPart);
        void checkActiveEddy(double time);
        void getRelativeVel(double time);
    
        // ------------- added for reacting particles
    
        void set_iyPos();
        void setPartAndGasMassSource(vector<double> & gasMassSourceFromEachFakeParticle);
        
        void setPartEnthSource(vector<double> & lambda_gas,
                           vector<double> & Gvel);
        
        double getHeatCapFromTemp(double temp);
    
	    double limitTimeStepForHistories( double timeNow, double proposedTimeStep );
        void storeHistories( double time );
        void adjustEddyVelHistories( int iPart, 
                double timeOccurEddy, 
                double tInt, 
                double tPart, 
                double vEddy, 
                double vPrev );
        void adjustEddyVelHistoriesTracers( int iPart, 
                double timeOccurEddy, 
                double tInt, 
                double vEddy ); 
        void outputHistories( string fname );
        double getSurfTemp(int & i);
        double getPartSurrGasTemp(int & i);
        double getTotalParticleEnthalpy();
        double getTotalParticleMass();
        double getParticleEnthalpy(int &i);
        double getParticleMass(int & i);
        double getFakeParticleMass(int & i);
        double getFakeParticleEnth(int & i);
        double getFakeParticleIntensiveEnth(int & i);
        double setPartTemp(int i, int j);
        void   readUVelProfile(vector<double> &x, vector<double>&y, vector<double> &X, vector<double>&Y);
        void   moveUVelProfile();
        
    

    private:
        double computeHeatTranCoeff(int &    i,
                                    double & lambda_gas,
									double & gRelVel);
        double computeBlowingFactor (int & i, double& h_loc);
        double getThermalCond       (double & temp);
        vector<double> linspace     (double x1, 
                                     double x2, 
                                     int n);
        void setHeatFlux(vector<double> & lambda_gas,
                         vector<double>& Gvel);
        double getHeatRateDueToConvec(double & lambda_gas, 
                                              double & Gvel, 
                                              int & i);
        double getHeatRateDueToRad(int & i);
        double getHeatFluxDueToConvec(double & lambda_gas,
                                            double & Gvel,
                                            int & i);
        double getHeatFluxDueToRad(int & i);
        
        
    ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:

        particles(int pnp, 
                  odtParam *odtP_p,
                  odtline  *oline,
                  string fileName);
        particles(const particles &part);
        void operator=(const particles &part);
        ~particles(){}

};

#endif



