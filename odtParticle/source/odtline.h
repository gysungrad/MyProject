/**
 * @file odtline.h
 * Header file for classes odtline 
 */

#ifndef ODTLINE_H
#define ODTLINE_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "anyline.h"
#include "streams.h"
#include "inputFile.h"
#include "adaptMesh.h"
#include "particles.h"
#include "ETA.h"
#include "MOM.h"

#ifdef DOCANTERA
#ifndef CANTERA21
#include "Cantera.h"
#endif
#include "IdealGasMix.h"
#include "transport.h"
#endif

#include "cantera_shell_functions.h"

using Cantera::Transport;
#ifdef CANTERA18
using Cantera_CXX::IdealGasMix;
#else
using Cantera::IdealGasMix;
#endif

using namespace std;
///////////////////////////////////////////////////////////////////////////////

/** A child class of anyline: use for combustion etc.  
 *  This line contains energy and composition (e.g.
 *  for combustion.
 *  It is separate from the odtline, which is used to hold velocities and
 *  perform eddy selection.  The grids may be different here than on the
 *  odtline and possibly other lines.
 *  
 *  @author David O. Lignell
 */

class odtline : public anyline {

  ////////////////////// DATA MEMBERS /////////////////////

    public:

        double                            momGas;
        double                            momPart;
        double                            momTotal;
        double                            ergGas;
        double                            ergPart;
        double                            ergTotal;

        vector<double>                    uvel;          ///< primary component
        vector<double>                    vvel;          ///< primary component
        vector<double>                    wvel;          ///< primary component
        
        vector<double>                    enth;          ///< J/kg, in props
        vector<vector<double> >           yspc;          ///< in props (yspc[k][i] is at pt i, sp k)
        double                            pres;          ///< Pressure (Pa)
        
        vector<double>                    temp;          ///< Auxiliary variable (Temperature)
        vector<double>                    mixf;          ///< Auxiliary variable (Mixture Fraction)
        vector<double>                    chi;           ///< Auxiliary variable (chi scalar dissipation rate)
        vector<double>                    gradZ;         ///< Auxiliary variable (mix frac gradient)
        vector<double>                    diffZ;         ///< Auxiliary variable (mix frac diffusion rate)

        vector<vector<double> >           mom;           ///< in props (mom[k][i] is at pt i, mom k)
        vector<vector<double> >           eta;           ///< in props (eta[k][i] is at pt i, param k)
        int                               nmom;          ///< number of moments
        int                               neta;          ///< number of paramters 
  
        int                               nscl;          ///< number of scalars 
        vector<vector<double> >           scl;           ///< in props (scl[k][i] is at pt i, scl k)
  
        IdealGasMix                       *gas;          ///< Cantera object
        Transport                         *tran;         ///< Cantera transport object
        streams                           *strm;         ///< pointer to stream object
  
        vector<string>                    spNames;       ///< species names
  
        bool                              LhasRxn;       ///< flag indicates h,Yi
        bool                              Lprxn;         ///< flag for particle reaction
        bool                              LhasVel;       ///< flag indicates u,v,w
        bool                              LhasTemp;      ///< flag indicated T (temperature)
        bool                              LhasMom;       ///< flag indicates particle moments
        bool                              LhasEta;       ///< flag indicates transported params (mixf heatloss chi etc)
        bool                              LhasScl;       ///< flag indicates additional scalars
        adaptMesh                         meshAdapter;   ///< mesh adapter object
        vector<double>                    voidFrac;      ///< cell void fraction (for particle cases)
        ETA*                              etaTools;      ///< Pointer to an eta object
        MOM*                              momTools;      ///< Pointer to a mom object

        int                               unifCompress;  ///< opposed jet: uniform compression or velocity based
        int                               nonPremixed;   ///< opposed jet: non-premixed=1, premixed=0
        double                            mixf_premixed; ///< opposed jet: mixture fraction for pre-mixed flow
        std::vector<std::vector<double> > Le_s;          ///< opposed jet: Lewis number of species. 	thermal diffusion / 	molecular diffusion
        std::vector<double>               Pr;            ///< opposed jet: Prandalt number of mixture	momentum transport/ 	heat transport
        std::vector<std::vector<double> > Sc_s;          ///< opposed jet: Schmidt number of species	momentum transport/ 	molecular diffusion
        std::vector<double>               heatLoss;      ///< heat loss
        // moved lambda to anyline (Falko)
        //std::vector<double>               lambda;        ///< thermal conductivity (specific for table lookup)

        int                               probType;     ///< ?

    private:

        inputFile                         caseInput_;  ///< inputFile object wrapping odtline input file

        ////////////////////// MEMBER FUNCTIONS  /////////////////////

    public:

        void                getYspVecAtPt(const int &i, vector<double> &yi);
        void                setTempVec();
        void                setRhoVec();
        void                setViscosity();
        void                setMixfVec();
        void                setChiVec();
        void                setDiffZVec();
        vector<double>      getMMW();
        void                normalizeSpecies(int ipos=-1);
        void                enforceYsBounds(std::vector<std::vector<double> > &Y_old);
        void                setVoidFrac(particles *part);

        virtual void        readProperties(  string fname);
        virtual void        outputProperties(string fname);
        virtual void        setAdptVar();
        double              getThermalDiff(int & i); 
        void                setOdtline(string fname);
        double              getTotalGasEtaEnthalpy();
        double              getTotalGasMass();

        void                setBCprops();
        void                expand(const double delz_in, const char LUC);

    private:
        void                setMesher(bool Ladpt);

        ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:
        
        odtline(int npts, double Ld, IdealGasMix *cantIG, 
                                     Transport *cantTran, 
                                     streams *strm_p, 
                                     odtParam *odtP_p,
                                     string caseSetupFile,
                                     bool LhasVel_p, bool LhasRxn_p);
        odtline(const odtline &odtl);                     // copy constructor
        odtline(const odtline *odtl, 
                int &i1, int &i2, bool Lwrap=false, bool Ladpt=true);      // copy constructor subset

        void initFromCopy(const odtline *odtl, 
                int &i1, int &i2, bool Lwrap=false, bool Ladpt=true);      // copy constructor subset


        void operator=(const odtline &odtl);
        
        ~odtline(){
            
        }

};


#endif

