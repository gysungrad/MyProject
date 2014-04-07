/**
 * @file diffuser.h
 * Header file for class diffuser
 */

#ifndef DIFFUSER_H
#define DIFFUSER_H

#include "odtline.h"
#include "odtParam.h"
#include "stats.h"
#include "batchRXR.h"
#include "dumper.h"
#include "radiation.h"
#include "particles.h"

using namespace std;
///////////////////////////////////////////////////////////////////////////////

/** Diffuser and reactor class advancing transport on the lines.  Transport 
 *  and reaction here occurs between eddy events.
 *  
 *  @author David O. Lignell
 */

class diffuser {

  ////////////////////// DATA MEMBERS /////////////////////

    public:

        vector<vector<double> >            k1rhs;           ///< rhs of all eqns, for multistage integ: [stage][var]
        vector<vector<double> >            rhsTrn;          ///< Transport part of rhs for all props
        vector<vector<double> >            rhsSrc;          ///< Source term part of rhs for all props
        vector<vector<double> >            rhsPartSrc;      ///< particle source term for two-way coupling
        vector<vector<double> >            flxProp;         ///< property fluxes [prop][grid]
        int                                neqns;           ///< number of equations solved (nprops)

        vector<vector<double> >            P1rhs;           ///< particle rhs function [y,u,v,w][ipart]
        
        int                                iptUvel;          ///< pointers to variables in krhs
        int                                iptVvel;
        int                                iptWvel;
        int                                iptTemp;
        int                                iptEnth;
        int                                iptYspc;
        int                                iptMom;
        int                                iptEta;
        double                             uBulk;            ///< the current bulk velocity in u direction calculated by odtStats
        double                             timeFS;

        vector<double>                     dxML;             ///< grid cell sizes of mixing line (odt line)
        
        vector<double>                     mmw;              ///< vector of mean molecular weights
        vector<double>                     invGamma;         ///< vector of 1/gamma
        vector<double>                     radSource_G;      ///< gas radiation source (W/m3)
        vector<double>                     radSource_P;      ///< particle radiation source (W/particle)
        vector<double>                     Gvel;             ///< gas velocity (for particles)

        vector<double>                     visc_f;           ///< face values of mu (dynamic viscosity)
        vector<double>                     lambda_f;         ///< face values of thermal conductivities
        vector<vector<double> >            DmixYs_f;         ///< face values of species diffusivities
        vector<vector<double> >            hsp;              ///< cell centered species enthalpies
        vector<vector<double> >            rrSpc;            ///< cell centered reaction rates, kg/m3*s
        
        double                             dpdt;             ///< dp/dt when needed

        double                             dtStepCFL;        ///< basic stepsize
        double                             Class_dtStep;     ///< dtStep using for diffuse single step

        int                                iN2;              ///< index of N2 in system

        odtParam                           *odtP;            ///< pointer to odtParam object
        odtline                            *odtl;            ///< pointer to odt line to diffuse
        particles                          *part;            ///< pointer to lagrangian particle object
        batchRXR                           brxr;             ///< integrates chemistry implicitly per cell
        dumper                             dmpr;             ///< object to output at specified times
        radiation                          rad;              ///< radiation module: optically thin and two-flux models
        vector<vector<double> >            *gSource;         ///< gas source term from particle reaction

        class eddyCombination {                                  ///< for heat transfer coefficient (type C interaction)
            public:
            double eddyLeftEdge;
            double eddyRightEdge;
            double eddyLength;
            double eddyLife;
            double eddyDeadTime;
            double eddyUvel;
            double eddyVvel;
            double eddyWvel;
            
            double eddyPvel;

            vector<int> eddyiPart;
            vector<int> eddyiPartIC;
            vector<double> eddyxPartPos;
            vector<double> eddyzPartPos;
        };
        
        vector<eddyCombination> ActiveEddy;

#ifdef IMPLICIT
        vector<vector<double> >            A;                ///< the lower diagonal of the tridiagonal system
        vector<vector<double> >            B;                ///< the middle diagonal of the tridiagonal system
        vector<vector<double> >            C;                ///< the upper diagonal of the tridiagonal system
        vector<vector<double> >            rhsImp;           ///< the right hand side of the tridiagonal system
#endif


//subdomain decomposition
	double dx_exp;
        double posNextSub_start;
        int indNextSub_start;
	double numberofsteps;             
	int subLeft_bc;				///< 0:Neumann; 1: Dirichlet
	int subRight_bc;			///< 0:Neumann; 1: Dirichlet
	
        ////////////////////// MEMBER FUNCTIONS  /////////////////////

        void diffuseOdtLine(double dtDiffuse, double time, stats &odtStats, const int &iStat);
        void diffuseOdtLine(double dtDiffuse, double time);
        void setGridSize();
        void updateGrids();

		double linearInterpToFace(int iface, const anyline & anyl,
                                  const vector<double> & vec);
        
        void updateActiveEddy(double leftEdge, double rightEdge, double eddySize, 
                            double eddyLife, double uEddy, double vEddy, double wEddy,
                            double eddyEndTime, double Pvel, vector<int> iPartInEddy,
                            vector<int> iPartInEddyIC, vector<double> xPartPos, vector<double> zPartPos
                            );
        void updatePartIxn(double time);
        void checkActiveEddyeffect(double time);
        void computeXZpartPos(double dtStep);
        void getEddyUWvel(double time);
        void getPartEddyVel(double time);
        
        void getMomTotal(odtline &line);

        void setRadHeatSources();
        void resetVarSizes();
        void setDiffusivities_HSP_RR();
        void setTimeStep();
        double setGasVelocity_return_dpdt();
	void set_diffuser(odtline *odtlp);

    private:

        void rhsf(vector<vector<double> > &krhs,
                  vector<vector<double> > &Prhs, double dtStep, double time);
        void setRhsTrn();
        void setRhsSrc(double dtStep);
        void setRhsPartSrc(vector<vector<double> > &Prhs);
        void getMomTotal();
        void outputAtSetTimes(double time, double timeD);
        void diffuseSingleStep(double &dtStep, double time); 
        
        void computeFluxes();

        void interpDiffusionCoeffToFaces(const anyline & anyl, vector<double> & diffc);
        void chopOutflowGrid();
        void gridInflow(double dtStep);
        void adaptGridsIfNeeded();

        void computeMatrix(vector<vector<double> > &A, vector<vector<double> > &B,
                vector<vector<double> > &C, vector<vector<double> > &rhs, double &dtStep);
        

        void updateRhoAndGrid (double & dtStep, vector<double> & uOld);         // enforcing continuity equation
	
	void   opposedJets(double timeStep);
	
	void strangSplitting(double dtStep, double time);

	void computeImplDiff(int k, double dtStep, std::vector<std::vector<double> > &flxProp );
    ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:

        diffuser(odtline *odtlp, odtParam &odtpp, particles *partp);  
        ~diffuser(){
            if(gSource) delete gSource;
        }

};

#endif
