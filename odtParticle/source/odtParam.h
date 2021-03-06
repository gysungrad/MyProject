/**
 * @file odtParam.h
 * Header file for class odtParam
 */

#ifndef ODTPARAM_H
#define ODTPARAM_H

#include "inputFile.h"
#include <string>
#include <vector>
#include <map>

///////////////////////////////////////////////////////////////////////////////

/** A class containing all the basic parameters of the odt simulation.  Other
 * parameters are included in the streams class and odtline::setodtline.
 * Most parameters are user defined (with defaults), others are derived.
 *  
 *  @author David O. Lignell
 */ 

class odtParam {

   public:
    
    //////////////////CONSTANTS/////////////////////////////////

        static const int BCTYPE_DIRICHLET   = 0;    ///< bctype dirichlet
        static const int BCTYPE_PERIODIC    = 1;    ///< bctype periodic
        static const int BCTYPE_WALL        = 2;    ///< bctype wall
        static const int BCTYPE_OUTFLOW     = 3;    ///< bctype outflow
        static const int BCTYPE_INLET       = 4;    ///< bctype inlet
        static const int BCTYPE_WALLOUTFLOW = 5;    ///< bctype wall on left, outflow on right

        static const int IMOMTYPE_DEFAULT     = 1;  ///< bctype wall on left, outflow on right
        static const int IMOMTYPE_SOOT        = 2;  ///< bctype wall on left, outflow on right
        static const int IMOMTYPE_AQUEOUS     = 3;  ///< bctype wall 

        static const int IETATYPE_DEFAULT     = 1;  ///< bctype wall on left, outflow on right
        static const int IETATYPE_TABLELOOKUP = 2;  ///< bctype wall on left, outflow on right
        static const int IETATYPE_AQUEOUS     = 3;  ///< bctype Wall

        enum TableTypes { TT_NONE, TT_EQUILIBRIUM, TT_FLAMELET };


  ////////////////////// DATA MEMBERS /////////////////////

    public:

        int                 nOdtReals;       ///< number of ODT realizations to do
        int                 nStat;           ///< number of statistic intervals
        int                 nTseg;           ///< number of time segments
        int                 seed;            ///< random number seed
        double              tEnd;            ///< end simulation time  
        double              Pmax;            ///< maximum eddy acceptance probability
        double              Pav;             ///< Average eddy acceptance probability
        double              dtfac;           ///< maximum factor to increase dtSmean
        double              tdfac;           ///< factor between dtCUmax and dtCFL for temporal flows; DEFAULT = 1
        int                 nDtSmeanWait;    ///< number of eddy samples before increase dtSmean
        double              domainLength;    ///< length of the domain (see also anyline::Ldomain)
        int                 ngrd_0;          ///< initial number of grid points
        double              rho_0;           ///< default density
        double              visc_0;          ///< default dynamic viscosity (molecular viscosity [Pa * s])
        double              lambda_0;        ///< default thermal conductivity ([W / (m * K)])
        double              dPdx;            ///< pressure gradient in channel flow
        double              phase;           ///< the integer number of the starting phase (default 0.0)
        double              Z_third;         ///< suppress penalty coefficient for suppression test
        int                 LES_type;        ///< 0 off; 1 third; 2 elapsed time; 3 frac domain
        double              Z_param;         ///< viscouse penalty parameter
        double              A_param;         ///< Energy distribution parameter alpha
        double              C_param;         ///< eddy frequency parameter
        double              Lp;              ///< most prob eddy size (/Ldomain)
        double              Lmax;            ///< max eddy size (/Ldomain)
        double              Lmin;            ///< min eddy size (/Ldomain)
        bool                LconstProp;      ///< flag for constant property sim (rho, mu)
        double              diffCFL;         ///< factor (< 1) to multipy diffusive dt cfl
        double              gDens;           ///< grid Density for new mesher
        double              dxmin;           ///< smallest allowed grid dx / domain size
        double              dxmax;           ///< largest allowed grid dx  / domain size
        double              largeGradFrac;   ///< fraction of max delta phi (split cells greater)
        double              smallGradFrac;   ///< fraction of max delta phi (merge cells smaller)
        double              largeCurvFrac;   ///< fraction of max delta theta (split cells greater)
        double              smallCurvFrac;   ///< fraction of max delta theta (merge cells smaller)

        std::string         chemMechFile;    ///< chemical mechanism file name
        std::string         restartOdtL;     ///< restart OdtL file name
        std::string         restartPart;     ///< restart Particle file name
        std::string         caseInp;         ///< case input file name
        std::string         partInp;         ///< particle input file name
        std::string         tableInp;        ///< table input file name

        int                 Iparticles;      ///< flag to turn on/off lagrangian particles
        bool                Ltwoway;         ///< two-way coupling for momentum exchange between particles and gas
        bool                Lrxn;            ///< flag to turn on/off chemical reaction
        int                 Iradiation;      ///< radiation model: -1 off, 1 opthin, 2 two flux
        bool                LimplicitChem;   ///< flag to turn on/off implicit chemistry solve 
        bool                LsecondOrder;    ///< flag to turn on second order integration in time
        bool                Lstrang;         ///< flag turn on/off strang splitting of diffusion and chemistry in species equation with implicit diffusion and chemistry (should be used for premixed to run stable)
        bool                LlimitMassFrac;  ///< switch: 0 turn off, 1 turn on limiting of mass fraction, 2 turn on with error message
        bool                Lperiodic;       ///< e.g., | . | . | . |,
                                             ///< Ldomain = 3, posf[0]=0,
                                             ///< posf[3] = 3, first and last
                                             ///< faces overlap
        bool                LhasTemp;        ///< flag for temperature
        bool                Lfalko;          ///< temporary flag for Falko. Used in diffuser.cc. DEFAULT = false
        bool                LheatedChannel;  ///< temporary debug flag for Falko.
        bool                Ldebug;          ///< debugging option
        int                 odtPprobType;    ///< switch for probType; 0: default; 6: premixed combustion (-> same probType as in premixed.inp); 7: premixed combustion including calculation of mixture fraction

        bool                Lsubdomain;      ///< flag turn on/off subdomain decomposition (presently LEM only)
	double		    dtGatherSubdomains; ///< time step for gathering subdomains and write output and shift if  needed
	int                 numSubdomains;   ///< number of subdomains

        int                 eddyMinCells;    ///< ark
        double              DAtimeFac;       ///< ark
        int                 sLastDA;         ///< size of the lastDA vector, presently hard wired
        int                 nsgrd;           ///< number of points on the stats grid
        bool                Lrestart;        ///< 1 to restart, 0 otherwise
        int                 LperiRestart;    ///< periodic restart; 0 no restart; n restart after periode n. Only the data of periode n is saved, the calculation advances to the end. The next realisation starts from the saved periode.
        bool                LmultiPhase;     ///< flag: false = single phase   true = multi phase
        
        int                 bcType;          ///< 0 for dirichlet, 1 for periodic
        int                 shiftMethod;     ///< 0: no shifting; for premixed: 1: const speed; 2: flame speed calculated from species; 3: species mass target

        int                 modDump;         ///< number of accepted eddies before dump and output File
        int                 modDisp;         ///< number of accepted eddies before display to screen
        int                 modActv;         ///< number of accepted eddies before read the active file
        bool                Lbinary;         ///< binary output with netcdf   

        double              eSurfTens;       ///< surface tension (J/m2)

        bool                Llem;            ///< flag to turn on/off lem 
        double              Dt;              ///< turbulent diffusivity to calculate lambda for lem 

        bool                Llaminar;        ///< laminar 
        
        double              uBClo;           ///< lower bc u velocity
        double              uBChi;           ///< upper bc u velocity
        double              vBClo;           ///< lower bc v velocity
        double              vBChi;           ///< upper bc v velocity
        double              wBClo;           ///< lower bc w velocity
        double              wBChi;           ///< upper bc w velocity
        double              tempBClo;        ///< lower bc temperature
        double              tempBChi;        ///< upper bc temperature

        int                 nPjump;          ///< size of pJump vector (= anyline::nprops)
        std::vector<double> pJump;           ///< periodic jumps on variables


        bool                LalignwGrav;
        double              Grav;
        std::vector<int>                 sclBCType;     ///< scl[]
        std::vector<double>              sclPr;
        std::vector<double>              sclbeta;
        std::vector<double>              scldelta;
        std::vector<double>              sclref;
        std::vector<double>              sclBClo;
        std::vector<double>              sclBChi;
        std::vector<double>              flxSclBClo;
        std::vector<double>              flxSclBChi;


        int                 nSclBClo;       ///< number of scalar lower BCs
        int                 nSclBChi;       ///< number of scalar upper BCs
        int                 nFlxSclBClo;    ///< number of scalar flux lower BCs
        int                 nFlxSclBChi;    ///< number of scalar flux upper BCs


        int                 nInletYsp;       ///< inlet species mass fractions
        std::vector<double> inletYsp;        ///< inlet species mass fractions
        double              inletTemp;       ///< inlet gas temperature
        double              inletEquRatio;   ///< inlet equivalence ratio for premixed combustion
        int                 nInletFuel;      ///< flag turn on/off usage equivaence ratio and size of vector inletFuel
        int                 nInletOxidizer;  ///< flag turn on/off usage equivaence ratio and size of vector inletOxidizer
        std::vector<double> inletFuel;       ///< fuel mass fractions used for mixture calculated with equivalence ratio
        std::vector<double> inletOxidizer;   ///< oxidizer mass fractions used for mixture calculated with equivalence ratio

        int                 nFreeStreamYsp;  ///< inlet species mass fractions
        std::vector<double> freeStreamYsp;   ///< free stream species mass fractions
        double              freeStreamTemp;  ///< free stream gas temperature

        bool                Lbuoyant;        ///< 1 for bouyancy in diffuser 0 for not
        bool                Lspatial;        ///< 1 for spatial formulation, use bctype 3

        int                 Imom;            ///< number of particle moments
        int                 ImomType;        ///< 1 for default; 2 for soot
        int                 Ieta;            ///< number of transported parameters
        int                 IetaType;        ///< 1 for default; 2 for table lookup
        bool                Lprxn;           ///< flag for reacting particles

        int                 Iscl;            ///< number of additional scalars

        int                 ItableLookup;    ///< 1 table lookup (not an input flag: set by IetaType)

        double              shiftSpeed;      ///< in m/s if domain is shifted with constant speed (shiftMethod 1)
        int                 nflameSpeedSpecies;
        std::vector<std::string>    flameSpeedSpecies;  ///< Species used to calculate burning velocity (shiftMethod 2)
        std::string         specMassTarget;  ///< species mass target (shiftMethod 3)
        int                 timeInpSwitch;   ///< apply timeInput  0 = none, 1 = sinusoidal fluc. in phi, 2 = sinusoidal fluc. in temp; 3: time input file
        double              sinusFreq;       ///< freqency for sinusoidal time input
        double              sinusAmpl;       ///< amplitude for sinusoidal time input


//------------------ derived parameters
        /// @name derived parameters 
        //@{

        double              esdp1;           ///< -2*Lp*domainLength
        double              esdp2;           ///< Cmax - Cmin = exp(-2Lp/Lmax)-exp(-2Lp/Lmin)
        double              esdp3;           ///< Cmin
        double              esdp4;
        double              lemRateParam;    ///< rate parameter for LEM
        
        //@}

    private:

        std::map<std::string, int> input_map_;     ///< Maps variable names to line numbers for odtParam input file

        inputFile                  odtpInput_;     ///< Input file object that wraps odtParam input file
        inputFile                  bcInput_;       ///< Input file object that wraps boundary condition input file

        ////////////////////// MEMBER FUNCTIONS  /////////////////////

    public:

        void computeEddySizeDistParams();
        void unnormalizeGridAndEddyDistParams();
        void resetDomainLength(double domL);

    private:

        void readOdtParam(std::string fileName);
        void readBCinput(std::string fileName);
        void computeLemRateParam();
        ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public: 

        odtParam(std::string fileName, std::string bcFile);      // constructor


};

#endif
