/**
 * @file particles.cc
 * Source file for class particles
 */

#include "particles.h"
#include "processor.h"
#include "randomGenerator.h"
#include "inputFile.h"
#include "odtline.h"
#include "diffuser.h"
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std;

#ifdef CANTERA18
using Cantera_CXX::IdealGasMix;
#else
using Cantera::IdealGasMix;
#endif

extern processor proc;

///////////////////////////////////////////////////////////////////////////////
/** Constructor function
 * 
 *  @param pnp       \input number of particles
 *  @param odtP_p    \input odtParam object
 *  @param oline     \input #odtline object
 *  @param fileName  \input input file name for initializing particles
 */

particles::particles(int pnp, 
                     odtParam *odtP_p,
                     odtline  *oline,
                     string fileName) {

    nPart = pnp;
    
    pDens0.resize(nPart);
    
    odtP       = odtP_p;
    line       = oline;
    this->diff = diff;

    yPos.resize(nPart);
    crossBound.resize(nPart, 0); 
    uvel.resize(nPart);
    vvel.resize(nPart);
    wvel.resize(nPart);
    pRadi.resize(nPart);
    
    f.resize(nPart,1.0);
    TauP.resize(nPart,0.);

    pActive.resize(nPart, true);

    iyPos.resize(nPart, 0);
    fracC.resize(nPart, 0);
    
    PeddyUvel.resize(nPart, 0.);
    PeddyVvel.resize(nPart, 0.);
    PeddyWvel.resize(nPart, 0.);

    //----------- each particle gets a history
    
    historiesTime.resize(nPart);
    historiesYpos.resize(nPart);
    historiesVpar.resize(nPart);
    historiesTpar.resize(nPart);
    historiesVslip.resize(nPart);
    historiesTgas.resize(nPart);
    historiesZmix.resize(nPart);
    historiesChi.resize(nPart);
    historiesGradZ.resize(nPart);
    historiesDiff.resize(nPart);

    timeNextHistoryPoint = 0.0;
    deltaTimeHistory = 0.0;
    
    //------------

    AGx = -9.81;
    AGy = 0.0;
    AGz = 0.0;

    Ltracer    = false;                       // reset in setParticle
    Lballistic = false;                       
    Lhistories = false;                       // reset in setParticle
    initPartLoc= 0;                       
    initPartVel = 0;

    // ParamEddylife = 1.0;                      // reset in setParticle
    ParamEddylife.resize(nPart);                      // reset in setParticle

    yPosTracer_TM.resize(nPart);
    eddyInfo.resize(nPart);
    vHT.resize(nPart, 0.0);

    nInPseudoPart.resize(nPart);
    pRadi0.resize(nPart);
    particleMassSource.resize(nPart, 0.0);

    emiss = 1.0;
    Hf    = 0;
    // Hf    = -5344000;
    pCp   = 2500;

    Ndpsn1 = 0;
    Ndpsn2 = 0;

    //---------- set pointers for the gas source

    int iPtDumb = 3;

    iPtUvel = 0;
    iPtVvel = 1;
    iPtWvel = 2;

    if (odtP -> Ieta){
        iPtEta  = iPtDumb;
        iPtDumb = iPtEta + line -> neta;
    }

    if (odtP -> Imom){
        iPtMom  = iPtDumb;
        iPtDumb = iPtMom + line -> nmom;
    }

    if(!odtP->Ieta){
        iPtEnth=iPtDumb;
        iPtDumb++;
    }
    iPtMass = iPtDumb;
    
    //----------------- set particle position, velocity, ...

    initPartMass_i = 0;

    setParticles(fileName);

    if(!odtP->Lprxn)
        *proc.ostrm << "\nWARNING: particles without reaction, computing lambda with T=300 K, P=1 atm" << endl;
    
    //--------------------
    
    totalVolatileMass = 0.0;
    totalVolatileEnthalpy = 0.0;
    peakTemp = 0;
    LuseDiBlasiModel = false;
    LuseNunnsModel = true;
    totalHeatDueToConvec = 0.0; //--These varaibles are used to calculate the contribution to total heat.
    totalHeatDueToRad = 0.0;
    
    //Compute the basic time step for the particle
    //alpha = k/(rho*Cp)
    double alpha = 0.15/pDens0[0]/1300;
    dtStepCFL = pow(cell_size,2.0) / alpha;
    *proc.ostrm << endl << "Particle step size--> " << dtStepCFL;    
}

///////////////////////////////////////////////////////////////////////////////
/** Constructor function (copy constructor)
 * 
 *  @param part       \input particles object to copy
*/

particles::particles(const particles &part) {
    
    radial_ngrd           = part.radial_ngrd;
    face_locs             = part.face_locs;
    cell_locs             = part.cell_locs;
    cell_size             = part.cell_size;
    cellVolume            = part.cellVolume;
    cellArea              = part.cellArea;
    cellFaceHeatFlux      = part.cellFaceHeatFlux;
    pDens0                = part.pDens0;
    dtStepCFL             = part.dtStepCFL;
        
    odtP                  = part.odtP;   
    line                  = part.line;
    
    nPart                 = part.nPart;   
   
    yPos                  = part.yPos;
    crossBound            = part.crossBound;
    uvel                  = part.uvel;
    vvel                  = part.vvel;
    wvel                  = part.wvel;
    cellTemp              = part.cellTemp;
    cellMass              = part.cellMass;
    cellEnth              = part.cellEnth;
    pRadi                 = part.pRadi;
    cellDens              = part.cellDens;
    f                     = part.f;
    TauP                  = part.TauP;
   
    pActive               = part.pActive;
    
    iyPos                 = part.iyPos;
    fracC                 = part.fracC;
   
    AGx                   = part.AGx;
    AGy                   = part.AGy;
    AGz                   = part.AGz;
   
    Ltracer               = part.Ltracer;
    Lballistic            = part.Lballistic;
    Lhistories            = part.Lhistories;
    initPartLoc           = part.initPartLoc;
    initPartVel           = part.initPartVel;

    ParamEddylife         = part.ParamEddylife;
   
    yPosTracer_TM         = part.yPosTracer_TM;
    eddyInfo              = part.eddyInfo;
    vHT                   = part.vHT;
   
    PeddyUvel             = part.PeddyUvel;
    PeddyVvel             = part.PeddyVvel;
    PeddyWvel             = part.PeddyWvel;
   
    nInPseudoPart         = part.nInPseudoPart;
    pRadi0                = part.pRadi0;
    particleMassSource    = part.particleMassSource;
    emiss                 = part.emiss;
    Hf                    = part.Hf;
    pCp                   = part.pCp;
    iPtMass               = part.iPtMass;
    iPtUvel               = part.iPtUvel;
    iPtVvel               = part.iPtVvel;
    iPtWvel               = part.iPtWvel;
    iPtEta                = part.iPtEta;
    iPtMom                = part.iPtMom;
   
    historiesTime         = part.historiesTime;
    historiesYpos         = part.historiesYpos;
    historiesVpar         = part.historiesVpar;
    historiesTpar         = part.historiesTpar;
    historiesVslip         = part.historiesVslip;
    historiesTgas         = part.historiesTgas;
    historiesZmix         = part.historiesZmix;
    historiesChi          = part.historiesChi;
    historiesGradZ         = part.historiesGradZ;
    historiesDiff         = part.historiesDiff;
    timeNextHistoryPoint  = part.timeNextHistoryPoint;
    deltaTimeHistory      = part.deltaTimeHistory;
    
    totalVolatileMass     = part.totalVolatileMass;
    totalVolatileEnthalpy = part.totalVolatileEnthalpy;
    volOfFakePart         = part.volOfFakePart;
    peakTemp              = part.peakTemp;
    LuseDiBlasiModel       = part.LuseDiBlasiModel;
    LuseNunnsModel         = part.LuseNunnsModel;
       
    pLength               = part.pLength;
    pShape                = part.pShape;
    voidFrac              = part.voidFrac;
    cellMassSource        = part.cellMassSource;
    cellEnthSource        = part.cellEnthSource;
    initPartMass_i        = part.initPartMass_i;
    totalHeatDueToConvec  = part.totalHeatDueToConvec;
    totalHeatDueToRad     = part.totalHeatDueToRad;  
}

///////////////////////////////////////////////////////////////////////////////

/** Overloaded assignment operator.
 *
 *  @param part \input particles object to use for assignment.
 */

void particles::operator=(const particles &part) {

    radial_ngrd           = part.radial_ngrd;
        
    face_locs             = part.face_locs;
    cell_locs             = part.cell_locs;
    cell_size             = part.cell_size;
    cellVolume            = part.cellVolume;
    cellArea              = part.cellArea;
    cellFaceHeatFlux      = part.cellFaceHeatFlux;
    pDens0                = part.pDens0;
    dtStepCFL             = part.dtStepCFL;
       
    odtP                  = part.odtP;   
    line                  = part.line;
   
    nPart                 = part.nPart;   
   
    yPos                  = part.yPos;
    crossBound            = part.crossBound;
    uvel                  = part.uvel;
    vvel                  = part.vvel;
    wvel                  = part.wvel;
    cellTemp              = part.cellTemp;
    cellMass              = part.cellMass;
    cellEnth              = part.cellEnth;
    pRadi                 = part.pRadi;
    cellDens              = part.cellDens;
    f                     = part.f;
    TauP                  = part.TauP;
   
    pActive               = part.pActive;
  
    iyPos                 = part.iyPos;
    fracC                 = part.fracC;
   
    AGx                   = part.AGx;
    AGy                   = part.AGy;
    AGz                   = part.AGz;
   
    Ltracer               = part.Ltracer;
    Lballistic            = part.Lballistic;
    Lhistories            = part.Lhistories;
    initPartLoc           = part.initPartLoc;
    initPartVel           = part.initPartVel;

    ParamEddylife         = part.ParamEddylife;
    
    historiesTime         = part.historiesTime;
    historiesYpos         = part.historiesYpos;
    historiesVpar         = part.historiesVpar;
    historiesTpar         = part.historiesTpar;
    historiesVslip         = part.historiesVslip;
    historiesTgas         = part.historiesTgas;
    historiesZmix         = part.historiesZmix;
    historiesChi          = part.historiesChi;
    historiesGradZ         = part.historiesGradZ;
    historiesDiff         = part.historiesDiff;
    timeNextHistoryPoint  = part.timeNextHistoryPoint;
    deltaTimeHistory      = part.deltaTimeHistory;
       
    yPosTracer_TM         = part.yPosTracer_TM;
    eddyInfo              = part.eddyInfo;
    vHT                   = part.vHT;
   
    nInPseudoPart         = part.nInPseudoPart;
    pRadi0                = part.pRadi0;
    particleMassSource    = part.particleMassSource;
    emiss                 = part.emiss;
    Hf                    = part.Hf;
    pCp                   = part.pCp;
    iPtMass               = part.iPtMass;
    iPtUvel               = part.iPtUvel;
    iPtVvel               = part.iPtVvel;
    iPtWvel               = part.iPtWvel;
    iPtEta                = part.iPtEta;
    iPtMom                = part.iPtMom;
    totalVolatileMass     = part.totalVolatileMass;
    totalVolatileEnthalpy = part.totalVolatileEnthalpy;
    volOfFakePart         = part.volOfFakePart;
    peakTemp              = part.peakTemp;
    LuseDiBlasiModel       = part.LuseDiBlasiModel;
    LuseNunnsModel         = part.LuseNunnsModel;
       
    pLength               = part.pLength;
    pShape                = part.pShape;
    voidFrac              = part.voidFrac;
    cellMassSource        = part.cellMassSource;
    cellEnthSource        = part.cellEnthSource;
    initPartMass_i        = part.initPartMass_i;
    totalHeatDueToConvec  = part.totalHeatDueToConvec;
    totalHeatDueToRad     = part.totalHeatDueToRad;
    
}


///////////////////////////////////////////////////////////////////////////////

/** Output the particle properties.  
 *
 * @param fname \input output file name.
 */

void particles::outputProperties(string fname) {

//   setParticleTemperature();

   ofstream ofile(fname.c_str()); 
   if(!ofile) 
       *proc.ostrm << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;

   ofile << "# particles = " << nPart;
   ofile << endl;
   ofile << setw(19) << "# Particle#"
         << setw(19) << "2_pos"
         << setw(19) << "3_uvel"
         << setw(19) << "4_vvel"
         << setw(19) << "5_wvel"
         << setw(19) << "6_temp"
         << setw(19) << "7_fractionLeft"
         << setw(19) << "VHT[i]"
         << setw(19) << "gVel[iyPos[i]]";
         

   ofile << scientific;
   ofile << setprecision(10);
   
   for(int i=0; i<nPart; i++) {
       ofile << endl;
//       if(!pActive[i])
//           continue;
       double T_surf = getSurfTemp(i);
       double fractionLeft = getFakeParticleMass(i) / initPartMass_i;
       ofile << setw(19) << i 
             << setw(19) << yPos[i] 
             << setw(19) << uvel[i]
             << setw(19) << vvel[i]
             << setw(19) << wvel[i]
             << setw(19) << T_surf
             << setw(19) << fractionLeft
             << setw(19) << vHT[i]
             << setw(19) << diff->Gvel[iyPos[i]];
   }

   ofile.close();
}

///////////////////////////////////////////////////////////////////////////////

/** Reads properties from an input file.  Used for restarts, etc.
 *  Note, this does not change the propNames, or the property bounds, which 
 *  are retained from the constructor.
 *
 *  @param fname \input name of file to read.
 */

void particles::readProperties(string fname) {

   *proc.ostrm << endl << "# Reading particles property file " << fname << endl;

   ifstream ifile(fname.c_str()); 
   if(!ifile) {
       *proc.ostrm << "\n\n# ***************** ERROR OPENING FILE " << fname << endl << endl;
       exit(0);
   }
   
   string       s1;
   stringstream ss1;

   getline(ifile, s1);                        // read line "# grid point = 100"
   ss1.str(s1);
   ss1 >> s1 >> s1 >> s1 >> nPart;
   ss1.clear();
   getline(ifile, s1);                        // read header line # 1_pos ...

   yPos.resize(nPart);   
   crossBound.resize(nPart,0);
   uvel.resize(nPart);          
   vvel.resize(nPart);         
   wvel.resize(nPart);        
   cellTemp.resize(nPart);

   int idumb;                  
   for(int i=0; i<nPart; i++) {
       ifile >> idumb 
           >> yPos[i]
           >> uvel[i]
           >> vvel[i]
           >> wvel[i] 
           >> cellTemp[i][radial_ngrd-1];

	crossBound[i] = int(yPos[i] / line->Ldomain);
	yPos[i] = yPos[i] - crossBound[i] * line->Ldomain;
   }

   ifile.close();

}

////////////////////////////////////////////////////////////////////////////
/**
 * reads the initial burning solution. That is, a mixture fraction profile is set
 * and the corresponding particle profile is solved assuming the mixf profile comes
 * from the particle off-gassing. 
 */
void particles::readInitialProfile() {
    cout << endl << "Reading inital particle profile. Please wait." << endl;
    //string pathToFile = "../input/particles_homoTurb/getRestartFile/eta0_eta1_mFakePart_hFakePart_u.dat";
    string pathToFile = "../input/particles_homoTurb/getRestartFile/" + initProfileFile;
    ifstream ifile(pathToFile.c_str());
    if(!ifile){
        *proc.ostrm << endl << "\n\n ********ERROR OPENING FILE " << pathToFile << endl;
        exit(EXIT_FAILURE);
    }
    
    // ---------------------------Read in the necessary variables in appropriate vectors
    vector<double> mFakePart = vector<double> (line->ngrd,0.0);
    vector<double> hFakePart = vector<double> (line->ngrd,0.0);
    
    initX = vector<double> (line->ngrd,0.0);
    initU = vector<double> (line->ngrd,0.0);
    
    for(int i=0;i<line->ngrd;i++){
        
        ifile >> line->eta[0][i] >> line->eta[1][i] >> mFakePart[i] >> hFakePart[i] >> line->uvel[i];
        line->uvel[i]=0;
        initX[i] = line->pos[i];
        initU[i] = line->uvel[i];
    }
    //tools::plot(line->pos,line->uvel,"uvel");
    // ---------------------------Update the line properties
    line->etaTools->updateOdtLineVecs(true);
    
    // ---------------------------Update the particle properties
    for (int i = 0; i < line->ngrd; i++) {
        double totalMass = 0;
        for (int j = 0; j < radial_ngrd; j++) {
            cellDens[i][j] = mFakePart[i]/volOfFakePart;
            cellMass[i][j] = cellDens[i][j] * cellVolume[i][j];
            totalMass +=cellMass[i][j];
            cellEnth[i][j] = hFakePart[i];
            cellTemp[i][j] = setPartTemp(i,j);
            if(cellTemp[i][j]>(95+273))
                moisture_tracker[i][j] = moist_content/(100*radial_ngrd) * initPartMass_i/radial_ngrd;
            else
                moisture_tracker[i][j] = 0;
            
        }
    }
    diff->resetVarSizes();
    diff->setGridSize();
    diff->setDiffusivities_HSP_RR();
    diff->setTimeStep();
    diff->setGasVelocity_return_dpdt();
    
    setHeatFlux(diff->lambda_f,diff->Gvel);
    
}

void particles::readUVelProfile(vector<double> &x, vector<double>&y, vector<double> &X, vector<double>&Y){
//void particles::readUVelProfile(){
    
//    string pathToDir = "../input/particles_homoTurb/getRestartFile/";
//    double gridPoints = line->ngrd;
//    string pathToFile = pathToDir + "u_"+tools::convertNumberToString(gridPoints)+".dat";
//    ifstream ifile(pathToFile.c_str());
//    if(!ifile){
//        *proc.ostrm << endl << "\n\n ********ERROR OPENING FILE " << pathToFile << endl;
//        exit(EXIT_FAILURE);
//    }
//    for (int i=0; i<gridPoints;i++){
//        ifile >> line->uvel[i];
//    }
    
//    Interpolation
    int N2 = X.size();
    int n_old = x.size();

    int j = 0;
    for (int i = 0; i < N2; i++)
        if (X[i] <= x[j + 1] || j + 2 == n_old)
            Y[i] = y[j] + (y[j + 1] - y[j]) / (x[j + 1] - x[j]) *(X[i] - x[j]);
        else {
            j = j + 1; // makes it move on to the next course grid point
            i = i - 1; // Makes it redo the current step
        }
}

void particles::moveUVelProfile(){
    double stored = 0;
    double toCopy = 0;
    int    myIndex     = 0;
    for (int i=0; i<(int)initU.size(); i++){
        if (i+1 == (int)initU.size()){
            myIndex = 0;
        }
        else
            myIndex = i+1;
        if (i==0){
            toCopy = initU[i];
        }
        stored = initU[myIndex];
        initU[myIndex] = toCopy;
        toCopy = stored;
    }
}
    
///////////////////////////////////////////////////////////////////////////////

/** 
* Read input file and set particle velocity, position, etc.
*
* We now have an input option to chose among three particle initializations                                                      \n
*   initPartLoc = 0 : randomly distribute particles across the domain (DEFAULT)                                                  \n
*   initPartLoc = 1 : place particles in order across the domain (see particleLocationFraction and leftParticleLocation)         \n
*   initPartLoc = 5 : place particles one per grid point                                                                         \n
*   initPartLeftLocFrac : initialize particles with leftmost particle at this fraction of the domain                             \n
*   initPartLocFrac : initialize particles across this fraction of the domain                                                    \n
*
* @param fname \input particle file name
*/

void particles::setParticles(string fname) {

    // double pDensInitial;           // local uniform property read from file to populate the vector pDens0;
    vector<double> pDensInitial;           // local uniform property read from file to populate the vector pDens0;
    // double pRadiInitial;
    vector<double> pRadiInitial;
    double nInPseudoPartInitial;
    double pTempInitial;
    vector<double> ParamEddylifeInitial;

    inputFile partInput_;          // inputFile object wrapping particle input file
    
    partInput_.setFile(fname);
    string myString = "";

    partInput_.getParameter("Ltracer",             &Ltracer,              false);
    partInput_.getParameter("initPartLoc",         &initPartLoc,          0);
    partInput_.getParameter("initPartLeftLocFrac", &initPartLeftLocFrac,  0.0);
    partInput_.getParameter("initPartLocFrac",     &initPartLocFrac,      1.0);
    
    partInput_.getParameter("nPartInitLoc",        &nPartInitLoc,         0);
    if (nPartInitLoc > 0 && nPartInitLoc == nPart) partInput_.getVector("nPartInitLoc",     &partInitLoc);
    
    partInput_.getParameter("initPartVel",         &initPartVel,          0);
    partInput_.getParameter("initPartUvel",        &initPartUvel,         0.0);
    partInput_.getParameter("initPartVvel",        &initPartVvel,         0.0);
    partInput_.getParameter("initPartWvel",        &initPartWvel,         0.0);
    partInput_.getParameter("Lballistic",          &Lballistic,           false);
    partInput_.getParameter("part_Length",         &pLength,              1.);      //m
    partInput_.getParameter("part_Shape",          &pShape,               2);
    partInput_.getParameter("void_Frac",           &voidFrac,             0.5);     //void fraction on the line
    
    partInput_.getParameter("LuniformPart",        &LuniformPart,       1);    
    partInput_.getParameter("nPartRadi",           &npRadiInitial,      1);
    partInput_.getVector("nPartRadi",              &pRadiInitial);          //m
    partInput_.getParameter("nPartDens",           &npDensInitial,      1);
    partInput_.getVector("nPartDens",              &pDensInitial);          //kg/m^3
    partInput_.getParameter("nParamEddylife",     &nParamEddylife,      1);
    partInput_.getVector("nParamEddylife",          &ParamEddylifeInitial);

    partInput_.getParameter("Ninitial",            &nInPseudoPartInitial, 1.);
    partInput_.getParameter("pTempInitial",        &pTempInitial,         298.15);  // K
    partInput_.getParameter("Lhistories",          &Lhistories,           false);
    partInput_.getParameter("deltaTimeHistory",    &deltaTimeHistory,     0.0);
    partInput_.getParameter("ngrd",                &radial_ngrd,          1);
    partInput_.getParameter("PeddyType",           &PeddyType,            3); 
    partInput_.getParameter("AGx",                 &AGx,                  0.0); 
    partInput_.getParameter("AGy",                 &AGy,                  0.0); 
    partInput_.getParameter("AGz",                 &AGz,                  0.0); 
    partInput_.getParameter("moist_content",       &moist_content,        10.0);      // percentage
    partInput_.getParameter("LInitProfile",        &LInitProfile,         false);
    partInput_.getParameter("initProfileFile",     &initProfileFile,      myString);
   
    face_locs = vector<vector<double> > (nPart, vector<double>(radial_ngrd + 1, 0.0));
    cell_locs = vector<vector<double> > (nPart, vector<double>(radial_ngrd, 0.0));
    cellVolume = vector<vector<double> > (nPart, vector<double>(radial_ngrd, 0.0));
    cellArea = vector<vector<double> > (nPart, vector<double>(radial_ngrd + 1, 0.0));
    cellFaceHeatFlux = vector<vector<double> > (nPart, vector<double>(radial_ngrd + 1, 0.0));
    
    moisture_tracker = vector<vector<double> > (nPart, vector<double>(radial_ngrd, 0.0));
    rho_M            = vector<vector<double> > (nPart, vector<double>(radial_ngrd, 0.0));
    heatOfVaporization = -2.44E6; //J/kg
    
    cellTemp = vector<vector<double> > (nPart, vector<double>(radial_ngrd, 0.0));
    cellMass = vector<vector<double> > (nPart, vector<double>(radial_ngrd, 0.0));
    cellEnth = vector<vector<double> > (nPart, vector<double>(radial_ngrd, 0.0));
    cellDens = vector<vector<double> > (nPart, vector<double>(radial_ngrd, 0));
    
    //----------------- check...
    
    if (Ltracer && Lballistic) {
        *proc.ostrm << "\nERROR: Ltracer and Lballistic are true, Terminating the program....." << endl;
        exit(0);
    }

    if (PeddyType != 1 && ParamEddylife.size() > 1) {
        *proc.ostrm << "\nERROR: multiple ParamEddylife ONLY apply for type-I interaction, Terminating the program....." << endl;
        exit(0);
    }

    if (ParamEddylifeInitial.size() < 1) {
        *proc.ostrm << "\nERROR: ParamEddylifeInitial size < 1. Go to particle.inp and set up ParamEddylife vector, Terminating the program....." << endl;
        exit(0);
    }
    
    if (pRadiInitial.size() < 1) {
        *proc.ostrm << "\nERROR: pRadiInitial size < 1. Go to particle.inp and set up pRadiInitial vector, Terminating the program....." << endl;
        exit(0);
    }

    if (pDensInitial.size() < 1) {
        *proc.ostrm << "\nERROR: pDensInitial size < 1. Go to particle.inp and set up pDensInitial vector, Terminating the program....." << endl;
        exit(0);
    }

    if (LuniformPart == 0) {
        if (nPart != npRadiInitial*npDensInitial*nParamEddylife) {
            *proc.ostrm << "\nERROR: when LuniformPart is false, nPart should be equal to npartRadi * npartDens * nParamEddylife " << endl;
            *proc.ostrm << "Check odtParam.inp and particle.inp, Terminating the program....." << endl;
            exit(0);
        }
    } 
    else if (LuniformPart == 1) {
        if (npDensInitial != 1 || npRadiInitial != 1 || nParamEddylife != 1) {
            *proc.ostrm << "\nERROR: when LuniformPart is 1, npartRadi = 1 && npartDens = 1 && nParamEddylife = 1 " << endl;
            *proc.ostrm << "Check particle.inp, Terminating the program....." << endl;
            exit(0);
        }   
    }
    else if (LuniformPart == 2) {
        if (nPart != npRadiInitial || nPart != npDensInitial || npRadiInitial != npDensInitial || nParamEddylife != 1) {
            *proc.ostrm << "\nERROR: when LuniformPart is true, nPart = npartRadi && nPart = npartDens && npartRadi = npartDens && nParamEddylife = 1 " << endl;
            *proc.ostrm << "Check odtParam.inp and particle.inp, Terminating the program....." << endl;
            exit(0);
        }
    }
    else {
        *proc.ostrm << "\nERROR: LuniformPart is 0 or 1 or 2 " << endl;
        *proc.ostrm << "Check odtParam.inp and particle.inp, Terminating the program....." << endl;
        exit(0);
    }

    //--------------------------

    //--------- Set the inital mass
    
    if (pShape == 2) //if particle shape is sphere
        volOfFakePart = 4.0 / 3.0 * M_PI * pow(pRadiInitial[0], 3.0); //TODO: mutiple particle properties (pDens[i]) do not apply here.
        // volOfFakePart = 4.0 / 3.0 * M_PI * pow(pRadiInitial, 3.0);
    else if (pShape == 1)//if particle shape is cylinder
        volOfFakePart = M_PI * pow(pRadiInitial[0], 2.0) * pLength; //TODO: mutiple particle properties (pDens[i]) do not apply here.
        // volOfFakePart = M_PI * pow(pRadiInitial, 2.0) * pLength;
    else{
        *proc.ostrm << endl << "Unrecognized particle shape--> " << pShape << " \n....Exiting out of the program. If you think there are"
                "insufficient flags, you can change me in file " << __FILE__ << " at line # " << __LINE__ << endl;
        exit(0);
    }
    // initPartMass_i = volOfFakePart * pDensInitial; 
    initPartMass_i = volOfFakePart * pDensInitial[0]; //TODO: mutiple particle properties (pDens[i]) do not apply here.
    
    //-----------------------------
    
    if ( initPartLeftLocFrac + initPartLocFrac > 1.0 ) {
        *proc.ostrm << "\nError in particles::setParticles. initPartLeftLocFrac + initPartLocFrac > 1.0 "
                       "and particles we be placed beyond domain boundaries" << endl;
        exit(0);
    }

    //------------ Initialize radius, density, and paramEddylife of ununiform particles 
        
    if (LuniformPart == 0) { // nonuniform particles
        // ununiform particles 
        // [#radius , #density , #paramEddylife]
        // for example, if 2 different radius, 3 different density, 4 different paramEddylife
        // then there are 2 * 3 * 4 = 24 different types of particles
        // [1, 2, 3] represents particle type has radius[1], density[2] and paramEddylife[3]
        // particle index is 1 * 3 * 4 + 2 * 4 + 3 = 23 

        for(int i = 0; i < npRadiInitial; i++)
            for (int j = 0; j < npDensInitial; j++)
                for (int k = 0; k < nParamEddylife; k++) {
                    int iPart = i * npDensInitial * nParamEddylife + j * nParamEddylife + k;
                    pRadi0[iPart] = pRadiInitial[i];
                    pRadi[iPart]  = pRadiInitial[i];
                    pDens0[iPart] = pDensInitial[j];
                    ParamEddylife[iPart] = ParamEddylifeInitial[k];
                }
    }

    //------------ Initialize particle positions to be randomly distributed

    randomGenerator      rr(-1);          // random generator object

    for(int i=0; i<nPart; i++)  {

        //------------ Initialize radius, density, and paramEddylife of uniform particles (1)
        if (LuniformPart == 1) { 
            pRadi0[i] = pRadiInitial[0];
            pRadi[i]  = pRadiInitial[0];
            pDens0[i] = pDensInitial[0];
            ParamEddylife[i] = ParamEddylifeInitial[0];
        }
        //------------ Initialize radius, density, and paramEddylife of uniform particles (2)
        if (LuniformPart == 2) { 
            pRadi0[i] = pRadiInitial[i];
            pRadi[i]  = pRadiInitial[i];
            pDens0[i] = pDensInitial[i];
            ParamEddylife[i] = ParamEddylifeInitial[0];
        }

        //------------ Initialize particles location 
        if ( initPartLoc == 0 ) {
            ///< randomly initialize the position of particles
            yPos[i] = odtP->domainLength*( initPartLeftLocFrac + initPartLocFrac * rr.getRand() );
            nInPseudoPart[i] = nInPseudoPartInitial;
        } 
        else if ( initPartLoc == 1 ) { ///< order the initial position of particles 
            yPos[i] = odtP->domainLength 
                * ( initPartLeftLocFrac + ( (double)(i)/ (double)(nPart-1) ) * initPartLocFrac );
            nInPseudoPart[i] = nInPseudoPartInitial;
        }
        else if ( initPartLoc == 2) { ///< initialize the position of particles in the middle of domain
            yPos[i] = 0.5 * odtP->domainLength;
            nInPseudoPart[i] = nInPseudoPartInitial;
        }
        else if ( initPartLoc == 3) {

            //doldbpart yPos[i] = ...;
            yPos[i] = partInitLoc[i];
            nInPseudoPart[i] = nInPseudoPartInitial;

            if (nPart != nPartInitLoc) {
                *proc.ostrm << "\nERROR: " << endl;
                *proc.ostrm << "\nparticle.inp: if iniPartLoc = 3, nPartInitLoc should be equal to Iparticles (odtParam.inp) and size of partInitLoc " << endl;
                exit(0);
            }
        }
        else if(initPartLoc == 5){
            /*Evenly distributed particles over the whole domain (with user selected void fraction)
              Number of particles must be equal to the initial number of grid points.*/
            yPos[i] = line->pos[i];
            line->voidFrac[i] = voidFrac;

            nInPseudoPart[i] = (1.0-line->voidFrac[i]) *(line->posf[i+1]-line->posf[i]) / 
                                       volOfFakePart;
            line->voidFrac[i] = 1.0;
            
        }
        else {
            *proc.ostrm << "\nERROR: Unknown initPartLoc" << endl;
             exit(0);
        }
        
        if (initPartLoc != 5)  
            line->setVoidFrac(this); //set and check the void fraction 
        
        iyPos[i] = line->linePositionToIndex(yPos[i], true);

    //---- initialize particle velocity -------------------
        
        if (initPartVel == 0) { // zero velocities
            uvel[i] = 0.0; 
            vvel[i] = 0.0;
            wvel[i] = 0.0;
        }
        else if (initPartVel == 1) { // line velocities
            uvel[i] = line->uvel[iyPos[i]]; 
            vvel[i] = line->vvel[iyPos[i]]; 
            wvel[i] = line->wvel[iyPos[i]]; 
        }
        else if (initPartVel == 2) { // read velocity values from particle.inp
            uvel[i] = initPartUvel; 
            vvel[i] = initPartVvel; 
            wvel[i] = initPartWvel; 
        }
        else {
            *proc.ostrm << "\nERROR: Unknown initPartVel" << endl;
            *proc.ostrm << "\nIn particle.inp, choose initPartVel as 0 for zero velocity; 1 for line velocity; 2 for hard coding or specify particle velocity (initUvel, initVvel, initWvel) in particle.inp" << endl;
             exit(0);
        }

        //---------------- set vars for particle internal heat equation solution
        
        face_locs[i] = linspace(0, pRadi[i], radial_ngrd + 1);
        if (pShape == 2) {
            cellArea[i][radial_ngrd] = 4.0 * M_PI * pow(face_locs[i][radial_ngrd], 2.0); //The last value
        } else if (pShape == 1) {
            cellArea[i][radial_ngrd] = 2.0 * M_PI * face_locs[i][radial_ngrd] * (face_locs[i][radial_ngrd] + pLength);
        }

        for (int j = 0; j < radial_ngrd; j++) {
            cellTemp[i][j] = pTempInitial;
            // cellDens[i][j] = pDensInitial;
            cellDens[i][j] = pDens0[i];
            cell_locs[i][j] = (face_locs[i][j] + face_locs[i][j + 1])*0.5;
            if (pShape == 2) {
                cellVolume[i][j] = 4. / 3. * M_PI * (pow(face_locs[i][j + 1], 3.0) - pow(face_locs[i][j], 3.0));
                cellArea[i][j] = 4.0 * M_PI * pow(face_locs[i][j], 2.0);
            } else if (pShape == 1) {
                cellVolume[i][j] = M_PI * (pow(face_locs[i][j + 1], 2.0) - pow(face_locs[i][j], 2.0)) * pLength;
                cellArea[i][j] = 2.0 * M_PI * face_locs[i][j] * (face_locs[i][j] + pLength);
            }
            cellMass[i][j] = cellVolume[i][j] * cellDens[i][j];
            cellEnth[i][j] = Hf + getHeatCapFromTemp(pTempInitial) * (pTempInitial - 298.15);

        }
        
    }
    cell_size = face_locs[0][1] - face_locs[0][0];

} 

///////////////////////////////////////////////////////////////////////////////////////// 
/** Diffuse tracer particles 
 *  Sets tracer particle locations and properties to be those of the grid.
 *  Move tracers by cell and fraction of original location in cell.
 *  Need to call particles::setFracC() first.  (Done in diffuser at the 
 *  beginning of the step, then advance the grid, then call this function to
 *  get the new location by simple stretching:
 *  \vc{
 *  | *    |
 *  |   *       |
 *  }
 *            Particle is the same fractional location before and after.
 *
 * Line velocity is the gas velocity not the ODT line velocity.  If you want to track
 * the ODT line velocity, then change here.
 *
 * @param Gvel \input vector of gas velocities (for particles)
 */
    
void particles::diffuseTracerParticle(vector<double> &Gvel) {    

    if(!Ltracer)
        return;
        
        for(int i=0; i<nPart; i++) {
            
            if(!pActive[i]) continue;

            if(!odtP->Lrxn) {
                iyPos[i] = line->linePositionToIndex(yPos[i],true);
                uvel[i] = line->uvel[iyPos[i]];
                vvel[i] = 0;
                wvel[i] = line->wvel[iyPos[i]];
            }   

        else {
            yPos[i] = line->posf[iyPos[i]] + fracC[i]*(line->posf[iyPos[i]+1]-line->posf[iyPos[i]]);
            iyPos[i] = line->linePositionToIndex(yPos[i],true);
            uvel[i] = line->uvel[iyPos[i]];
            vvel[i] = Gvel[iyPos[i]];
            wvel[i] = line->wvel[iyPos[i]];
        }   
    }
}

///////////////////////////////////////////////////////////////////////////////////////// 

/**
 *  Right hand side of langrangian particles.
 *  Update tracer and ballistic particles: make sure they work with ballistic particles and grid motion.
 *  Prhs comes from the diffuser.
 *
 * @param Gvel       \input  Gas velocity from the diffuser
 * @param Prhs       \output vector that holds the right hand side of the transport equations.
 * @param dtStep     \input  time step
 * @param lambda_gas \input  vector of thermal conductivity of the gas at each face
 * @param dxML       \input  Grid Size
 * @param gSource    \output vector that holds the source terms for gas from particles.
 *
 * Solving for the following equations
 *
 * Gas momentum source term
 * \f[
 *      S_{gas,\,mom} = \frac{uS_\phi}{\rho}
 * \f]
 * Gas enthalpy source term =
 * \f[
 *      S_{gas,\,enth} = -\frac{hS_\phi}{\rho} - \frac{1}{\rho}\sum\limits_{i=1} 
        \tilde{h}_i\theta_iA_{p,\,i}(T_g-T_{p,\,i})n_i 
 * \f]
 * Gas mixture fraction source term =
 * \f[
 *      S_{gas,\,mixf} = (1-{\xi})\frac{S_\phi}{\rho}
 * \f]
 */
void particles::computeRHSFAndSetGasSource(vector<double> &Gvel, 
                                           vector<vector<double> > &Prhs,
                                           double dtStep,
                                           vector<double> & lambda_gas, 
                                           vector<double> & dxML, 
                                           vector<vector<double> >* gSource,
                                           double time) {

    if (Ltracer)           // tracer particles done elsewhere
        return;

     if (!Lballistic && !Ltracer) 
         diff->getPartEddyVel(time);
    //--------------- loop over all particles
    
    for (int i = 0; i < nPart; i++){

        if (!pActive[i]) 
            continue;

        //-------------------  ballistic particles

        if (Lballistic){
            Prhs[0][i] = vvel[i];
            Prhs[1][i] = 0;
            Prhs[2][i] = 0;
            Prhs[3][i] = 0;
        }

    //-------------------  inertial particles

        else if(!Lballistic && !Ltracer) {

            set_f(i, PeddyUvel[i], PeddyVvel[i], PeddyWvel[i]); 
            set_TauP(i);

            Prhs[0][i] = (PeddyVvel[i] + AGy*TauP[i]/f[i])                 // y
                + TauP[i]/f[i]*(vvel[i]-PeddyVvel[i]-AGy*TauP[i]/f[i]) 
                * (1-exp(-dtStep*f[i]/TauP[i])) / dtStep;
            Prhs[1][i] = -(uvel[i] - PeddyUvel[i])*f[i]/TauP[i] + AGx;     // v
            Prhs[2][i] = -(vvel[i] - PeddyVvel[i])*f[i]/TauP[i] + AGy;     // u
            Prhs[3][i] = -(wvel[i] - PeddyWvel[i])*f[i]/TauP[i] + AGz;     // w
        }           // end if inertial particles calculation
    }               // end loop over particles

    //------------- add particle reaction stuff

    if (odtP->Lrxn) {
        cellMassSource = vector<vector<double> >(nPart,vector<double>(radial_ngrd,0.0)); // kg/s per particle
        cellEnthSource = vector<vector<double> > (nPart,vector<double>(radial_ngrd,0.0)); // J/kg*s per particle
        cellTotalEnthSource = vector<vector<double> > (nPart,vector<double>(radial_ngrd,0.0)); // J/s per particle
        vector<double> gasMassSourceFromEachFakeParticle = vector<double>(nPart,0.0);

        cellMoistureSource = vector<vector<double> >(nPart,vector<double>(radial_ngrd,0.0));
        rho_M        = vector<vector<double> >(nPart,vector<double>(radial_ngrd,0.0));

        if(odtP->Lprxn)
            setPartAndGasMassSource(gasMassSourceFromEachFakeParticle);

        //-------- Compute particle enthalpy source (Particle enthalpy source term should be computed even if there is no particle reaction)
        setPartEnthSource(lambda_gas,Gvel);

        totalHeatDueToConvec = 0.0;
        totalHeatDueToRad = 0.0;
        for (int i = 0; i < nPart; i++){

            if(!pActive[i])
                continue;

            //----------gas Mass source term

            double gasMassSourceFromEachParticle = 0.0;
            if(odtP->Lprxn){
                totalVolatileMass += dtStep * gasMassSourceFromEachFakeParticle[i] * nInPseudoPart[i];
                gasMassSourceFromEachParticle = nInPseudoPart[i] * gasMassSourceFromEachFakeParticle[i]
                    / (dxML[iyPos[i]] * line->voidFrac[iyPos[i]]);
                (*gSource)[iPtMass][iyPos[i]] += gasMassSourceFromEachParticle;
            }

            //----------- gas enthalpy source term (-S_phi/rho *(h_p + h_gas)+heatLossFromConvec + heatLossFromRad

            int enthalpyPointer=-1; 
            if(odtP->Ieta)
                enthalpyPointer=iPtEta+1;
            else if(!odtP->Ieta)
                enthalpyPointer=iPtEnth;

            (*gSource)[enthalpyPointer][iyPos[i]] += (getFakeParticleIntensiveEnth(i) - line->enth[iyPos[i]]) * gasMassSourceFromEachParticle
                / line->rho[iyPos[i]] - getHeatRateDueToConvec(lambda_gas[iyPos[i]], Gvel[iyPos[i]], i) 
                * nInPseudoPart[i] / (line->rho)[iyPos[i]]/(dxML[iyPos[i]]*line->voidFrac[iyPos[i]]); //(J/s kg)

            //----------gas Momentum source terms

            (*gSource)[iPtUvel][iyPos[i]] += -(uvel[i] + line->uvel[iyPos[i]]) * gasMassSourceFromEachParticle 
                / line->rho[iyPos[i]];
            (*gSource)[iPtVvel][iyPos[i]] += -(vvel[i] + line->vvel[iyPos[i]]) * gasMassSourceFromEachParticle
                / line->rho[iyPos[i]];
            (*gSource)[iPtWvel][iyPos[i]] += -(wvel[i] + line->wvel[iyPos[i]]) * gasMassSourceFromEachParticle 
                / line->rho[iyPos[i]];

            //----------gas moment source terms
            for (int k=0; k<line->nmom; k++){
                (*gSource)[iPtMom + k][iyPos[i]] += -(line->mom)[k][iyPos[i]] * gasMassSourceFromEachParticle
                    / line->rho[iyPos[i]];
            }
            //---------gas mixture fraction source term

            if(odtP->Ieta){ //only if Ieta flag is on
                if(line->eta[0][iyPos[i]] > 1.){
                    (*proc.ostrm).precision(10);
                    *proc.ostrm << endl << "line mixf is greater than 1. Correcting it from " 
                        << line->eta[0][iyPos[i]] << " to 1." << endl;
                    line->eta[0][iyPos[i]] = 1.;
                }

                (*gSource)[iPtEta][iyPos[i]] += (1.0 - (line->eta)[0][iyPos[i]]) * gasMassSourceFromEachParticle 
                    / (line->rho)[iyPos[i]];
            }
        } // end loop particles
    } // Lrxn end
}

///////////////////////////////////////////////////////////////////////////////

/** Calculate position index iyPos
 */
void particles::set_iyPos(){
    for (int i = 0; i < nPart; i++){
        if (!pActive[i]){
            continue;
        }

        iyPos[i] = line -> linePositionToIndex(yPos[i], true);
    }
}

///////////////////////////////////////////////////////////////////////////////

/** Calculate position index iyPos and nonlinear correaction factor f 
 */

void particles::set_f(int iPart, double Ug, double Vg, double Wg) { 

    if(odtP->ItableLookup)
        line->etaTools->updateOdtLineVecs();
    
    double diffu;
    double diffv;
    double diffw;
    double rootDiffVel;
    double ReSlip;

    iyPos[iPart] = line->linePositionToIndex(yPos[iPart], true);

    diffu = uvel[iPart] - Ug;
    diffv = vvel[iPart] - Vg;
    diffw = wvel[iPart] - Wg;
    rootDiffVel = sqrt(diffu*diffu + diffv*diffv +diffw*diffw);
    ReSlip = line->rho[iyPos[iPart]]*2*pRadi[iPart]*rootDiffVel/line->molec[iyPos[iPart]];

    if(ReSlip <= 1000.)
         f[iPart] = 1.+ 0.15*pow(ReSlip,0.687);
//          f[iPart] = 1.; //GYSun
    else 
        *proc.ostrm << "\nWARNING: LARGE REYNOLDS NUMBER > 1000 " << endl; 
}

///////////////////////////////////////////////////////////////////////////////

/** Calculate position index iyPos and Stokes # TauPi 
 */

void particles::set_TauP(int iPart) {
    
    double CunFac;
    double gasMeanFreePath;
    
    iyPos[iPart] = line->linePositionToIndex(yPos[iPart], true);
    
    if(odtP->Lrxn || odtP->Lprxn) 
        gasMeanFreePath = 1.380648813E-23*line->temp[iyPos[iPart]]/(sqrt(2.0)*cellArea[iPart][radial_ngrd]*line->pres);
    else
        gasMeanFreePath = 1.380648813E-23*300.0               /(sqrt(2.0)*cellArea[iPart][radial_ngrd]*101325.);
    
    CunFac = 1 + gasMeanFreePath / pRadi[iPart]* 1; ///Cunningham slip factor
    
    //Calculate average cell density
    double totalDens = 0;
    for (int j = 0; j < radial_ngrd; j++) {
        totalDens += cellDens[iPart][j];
    }
    double averageDens = totalDens/radial_ngrd;
    // TauP[iPart] = averageDens*4 * pRadi[iPart] * pRadi[iPart] * CunFac / 18 / line->molec[iyPos[iPart]];
    TauP[iPart] = pDens0[iPart]*4 * pRadi[iPart] * pRadi[iPart] * CunFac / 18 / line->molec[iyPos[iPart]];

}

///////////////////////////////////////////////////////////////////////////////

/** Check if particles are out of bounds and inactivate if they are or periodic case
 *
 * @param Gvel \input vector of gas velocities (for particles)
 *
 */

void particles::checkBoundsSetInactive(vector<double> &Gvel, double time) {

    for(int i=0; i<nPart; i++){

        if(!pActive[i]) continue;

        if( yPos[i] <= line->posf[0] && odtP->Lperiodic ) {
            yPos[i] += line->Ldomain; 
            crossBound[i]--;
            iyPos[i] = line->linePositionToIndex(yPos[i], true);
            if (Ltracer) {                          // periodic boudnary condition
                uvel[i] = line->uvel[iyPos[i]];
                wvel[i] = line->wvel[iyPos[i]];
                if(odtP->Lrxn) 
                    vvel[i] = line->vvel[iyPos[i]];
                else
                    vvel[i] = Gvel[iyPos[i]];
            }
        }
        else if( yPos[i] <= (line->posf[0]+ pRadi[i]) && odtP->bcType == 2 ) {

            if( Ldeposition ) { // wall-bounded and deposition boundary
                pActive[i] = false;
                yPos[i] = 0.0;     
                uvel[i] = 0.0;
                vvel[i] = 0.0;
                wvel[i] = 0.0;

                if(time >=0.0) {
                    Ndpsn1++;
                    *proc.ostrm << endl << "# Ndep left, tot, myid, time = " 
                        << Ndpsn1 << "  " << Ndpsn1+Ndpsn2 << "  " 
                        << proc.myid << "  " 
                        << setprecision(10)
                        << time; 
                }

            }
            else if( !Ldeposition) { // wall-bounded and elastic collision
                yPos[i] = 0.0;
                vvel[i] = -vvel[i]; 
                wvel[i] = -wvel[i]; 
            }
        }
        else if( yPos[i] <= line->posf[0] && !odtP->Lperiodic ) {
            // out flow 
            pActive[i] = false;
            yPos[i] = 0.0;     
        }
        else if( ( yPos[i] >= line->posf[0] + line->posf[line->ngrd] )
                && odtP->Lperiodic ) {
            yPos[i] -= line->Ldomain; 
            crossBound[i]++;
            iyPos[i] = line->linePositionToIndex(yPos[i], true);
            if (Ltracer) {
                uvel[i] = line->uvel[iyPos[i]];
                wvel[i] = line->wvel[iyPos[i]];
                if(odtP->Lrxn) 
                    vvel[i] = line->vvel[iyPos[i]];
                else
                    vvel[i] = Gvel[iyPos[i]];
            }
        }
        else if( ( yPos[i] >= ( line->posf[0] + line->posf[line->ngrd] - pRadi[i] ) )
                && odtP->bcType == 2 ) {
            if( Ldeposition ) { // wall-bounded and deposition boundary
                pActive[i] = false;
                yPos[i] = odtP->domainLength;
                uvel[i] = 0.0;
                vvel[i] = 0.0;
                wvel[i] = 0.0;

                if(time >=0.0) {
                    Ndpsn2++;
                    *proc.ostrm << endl << "# Ndep right, tot, myid, time = " 
                        << Ndpsn2 << "  " << Ndpsn1+Ndpsn2 << "  " 
                        << proc.myid << "  " 
                        << setprecision(10)
                        << time; 
                }
            }
            else if( !Ldeposition ) { // wall-bounded and elastic collision
                yPos[i] = odtP->domainLength;
                vvel[i] = -vvel[i]; 
                wvel[i] = -wvel[i]; 
            }
        }
        else if( ( yPos[i] >= line->posf[0] + line->posf[line->ngrd] )
                && !odtP->Lperiodic ) {
            // out flow 
            pActive[i] = false;
            yPos[i] = odtP->domainLength;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

/** 
 * Sets the fracC vector for later use. 
*/

void particles::setFracC() {

    if (odtP->Lrxn && Ltracer) {
        for (int i = 0; i < nPart; i++) {
            if (!pActive[i]) continue;
            fracC[i] = (yPos[i] - line->posf[iyPos[i]]) /
                    (line->posf[iyPos[i] + 1] - line->posf[iyPos[i]]);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////// 

/** 
 * Adjust time step to store particle histories at specified times
 */

double particles::limitTimeStepForHistories( double timeNow, double proposedTimeStep ) {

  if ( deltaTimeHistory < 1e-30 )  return proposedTimeStep; //don't want to set time step to zero
  
  //adjust time step if we are going too far with the proposed step.
  if ( timeNow + proposedTimeStep > timeNextHistoryPoint ) {
    return timeNextHistoryPoint - timeNow;
  } else {
    return proposedTimeStep;
  }
}

///////////////////////////////////////////////////////////////////////////////////////// 

/** 
 * Store particle histories
 */

void particles::storeHistories( double time ) {    

    if ( time < timeNextHistoryPoint ) return;
    timeNextHistoryPoint += deltaTimeHistory;

    for(int i=0; i<nPart; i++) {

        //if(!pActive[i]) continue;

        historiesTime[i].push_back( time );

        if(pActive[i]) {
            iyPos[i] = line->linePositionToIndex(yPos[i],true);
            historiesYpos[i].push_back( yPos[i] + crossBound[i] * line->Ldomain );
            historiesVpar[i].push_back( vvel[i] );
            if(odtP->Lrxn) {
                historiesVslip[i].push_back( vvel[i] );
                //historiesVgas[i].push_back( line->vvel[iyPos[i]] );
                historiesTgas[i].push_back( line->temp[iyPos[i]] );
                historiesZmix[i].push_back( line->mixf[iyPos[i]] );
                historiesChi[i].push_back( line->chi[iyPos[i]] );
                historiesGradZ[i].push_back( line->gradZ[iyPos[i]] );
                historiesDiff[i].push_back( line->diffZ[iyPos[i]] );
                //historiesTpar[i].push_back( cellTemp[i][0] );
                // cout << endl << "storehistory getSurfTemp(i) " << getSurfTemp(i) << endl;
                historiesTpar[i].push_back( getSurfTemp(i) );
            }  
        } else {
            iyPos[i] = line->linePositionToIndex(yPos[i],true);
            historiesYpos[i].push_back( yPos[i] + crossBound[i] * line->Ldomain );
            historiesVpar[i].push_back( 0.0 );

            if(!odtP->Lrxn) {
                historiesVslip[i].push_back( 0.0 );
                historiesTgas[i].push_back( -1.0 );
                historiesZmix[i].push_back( -1.0 );
                historiesChi[i].push_back( 0.0 );
                historiesGradZ[i].push_back(0.0 );
                historiesDiff[i].push_back( 0.0 );
                // historiesTpar[i].push_back( -1.0 );
            }
        }
    }
}

/////////////////////////////////////////////////////////////////////////////////////////

/** 
 * The eddy acts instantaneously in its displacement, 
 * but we will add the velocity that caused that 
 * displacement to the existing particle history.
 *  @param iPart input: particle number acted on by eddy
 *  @param time  input: eddy occurence time
 *  @param tInt  input: particle-eddy interaction time
 *  @param tPart input: particle time constant
 *  @param vEddy input: eddy velocity
 *  @param vPrev input: particle velocity prior to eddy
 *
 *  NOTE: As of 5/31/12 the timeOccurEddy is not aligned with the historiesTime vector: 
 *        there is a big gap waiting to be caught up in the diffusion step.  
 *        How should we fix this?  Perhaps set timeOccurEddy = historiesTime[iPart][size()-1] 
 */

void particles::adjustEddyVelHistories( int iPart, 
        double timeOccurEddy, 
        double tInt, 
        double tPart, 
        double vEddy, 
        double vPrev ) {

    //return;

    int j = historiesTime[iPart].size();
    if ( j == 0 ) return;

    double tInEddy;
    tInEddy = tInt - ( timeOccurEddy - historiesTime[iPart][--j] );
    
    while ( tInEddy > 0.0 ) {
        if ( j < 0 ) {
            cout << endl << "Warning: j < 0 in adjustEddyVelHistories" << endl;
	    break;
        } 
        else {
            historiesVpar[iPart][j] += vEddy * ( 1 - exp( - tInEddy / tPart ) );
	    //right now only need slip velocity for reacting flows
	    if ( odtP->Lrxn )
	      historiesVslip[iPart][j] += - vEddy * exp( - tInEddy / tPart ) ;
        }
        j--;
        tInEddy = tInt - ( timeOccurEddy - historiesTime[iPart][j] );
    };

}

/////////////////////////////////////////////////////////////////////////////////////////

/** 
 * The eddy acts instantaneously in its displacement, 
 * but we will add the velocity that caused that 
 * displacement to the existing particle history.
 *  @param iPart input: particle number acted on by eddy
 *  @param time  input: eddy occurence time
 *  @param tInt  input: particle-eddy interaction time
 *  @param tPart input: particle time constant
 *  @param vEddy input: eddy velocity
 *  @param vPrev input: particle velocity prior to eddy
 *
 *  NOTE: As of 5/31/12 the timeOccurEddy is not aligned with the historiesTime vector: 
 *        there is a big gap waiting to be caught up in the diffusion step.  
 *        How should we fix this?  Perhaps set timeOccurEddy = historiesTime[iPart][size()-1] 
 */

void particles::adjustEddyVelHistoriesTracers( int iPart, 
        double timeOccurEddy, 
        double tInt, 
        double vEddy ) {

    //return;

    int j = historiesTime[iPart].size();
    if ( j == 0 ) return;

    double tInEddy;
    tInEddy = tInt - ( timeOccurEddy - historiesTime[iPart][--j] );

    while ( tInEddy > 0.0 ) {
        if ( j < 0 ) {
            cout << endl << "Warning: j < 0 in adjustEddyVelHistories" << endl;
	    break;
        } 
        else {
	  historiesVpar[iPart][j] += vEddy;
        }
        j--;
        tInEddy = tInt - ( timeOccurEddy - historiesTime[iPart][j] );
    };

}

////////////////////////////////////////////////////////////////////////////////
// helper function 
// convert int to string in particles::outputHistories
string IntToString (int a)
{
    ostringstream temp;
    temp << a;
    return temp.str();
}

////////////////////////////////////////////////////////////////////////////////

/** Output the particle properties.  
 *
 * @param fname input: output file name.
 * by John Hewson
 * modified by Guangyuan Sun 03/2014
 */

void particles::outputHistories(std::string fname) {

    ofstream histfile(fname.c_str());
    if (!histfile.is_open())
        *proc.ostrm << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;

    histfile << "#" 
        << setw(18) << "time";
    for(int n = 0; n < nPart; n++) {
        string s1 = "Part_" + IntToString(n) + " 1_Ypos";
        string s2 = "Part_" + IntToString(n) + " 2_Vpar";
        histfile << setw(19) << s1 
            << setw(19) << s2; 
        if(odtP->Lrxn) {
            string s3 = "Part_" + IntToString(n) + " 3_Tpar";
            string s4 = "Part_" + IntToString(n) + " 4_Vslip";
            string s5 = "Part_" + IntToString(n) + " 5_Zmix";
            string s6 = "Part_" + IntToString(n) + " 6_Chi";
            string s7 = "Part_" + IntToString(n) + " 7_GradZ";
            string s8 = "Part_" + IntToString(n) + " 8_Diff";
            string s9 = "Part_" + IntToString(n) + " 9_Tgas";
            histfile << setw(19) << s3 
                << setw(19) << s4 
                << setw(19) << s5 
                << setw(19) << s6 
                << setw(19) << s7 
                << setw(19) << s8 
                << setw(19) << s9; 
        }
    }
    histfile << endl;

    int nHistTime = historiesTime[0].size();
    for(int i = 0; i < nHistTime; i++) { 
            histfile << setw(19) << historiesTime[0][i];         // history time
        for(int n = 0; n < nPart; n++) {
            histfile << setw(19) << historiesYpos[n][i]         // history Ypos 
                << setw(19) << historiesVpar[n][i];        // history Vpar
            if(odtP->Lrxn) {
                histfile << setw(19) << historiesTpar[n][i]     // history Tpar  
                    << setw(19) << historiesVslip[n][i]    // history Vslip
                    << setw(19) << historiesZmix[n][i]     // history Zmix
                    << setw(19) << historiesChi[n][i]      // history Chi
                    << setw(19) << historiesGradZ[n][i]    // history GradZ
                    << setw(19) << historiesDiff[n][i]     // history Diff
                    << setw(19) << historiesTgas[n][i];    // history Tgas
            }
        }
        histfile << endl;
    }
    histfile.close();

//     for( int n = 0; n < nPart; n++ ) {/*{{{*/
//         
//         stringstream ss1;
//         ss1 << n+1 << "_time.dat";
//         string    fnameTime = fname + ss1.str();
// 
//         ofstream ofileTime(fnameTime.c_str()); 
//         if(!ofileTime) 
//             *proc.ostrm << "\n\n***************** ERROR OPENING FILE " << fnameTime << endl << endl;
// 
//         //   ofileTime << "# particles = " << nPart;
// 
//         //-------------- iterate over times
// 
//         for( int i = 0; i < historiesTime[n].size() ; i++ ) { 
//             ofileTime << historiesTime[n][i] << endl ;
//         }
//         ofileTime.close();
// 
// 
//         stringstream ssY;
//         ssY << n+1 << "_Y.dat";
//         string fnameY = fname + ssY.str();
// 
//         ofstream  ofileY(fnameY.c_str()); 
//         if(!ofileY) 
//             *proc.ostrm << "\n\n***************** ERROR OPENING FILE " << fnameY << endl << endl;
// 
//         //------------- iterate over history output
// 
//         for ( int i = 0; i < historiesYpos[n].size() ; i++ ) { 
//             ofileY << historiesYpos[n][i] << endl ;
//         }
// 
//         ofileY.close();
// 
// 
//         stringstream ssV;
//         ssV << n+1 << "_V.dat";
//         string fnameV = fname + ssV.str();
// 
//         ofstream  ofileV(fnameV.c_str()); 
//         if(!ofileV) 
//             *proc.ostrm << "\n\n***************** ERROR OPENING FILE " << fnameV << endl << endl;
// 
//         //-------------- iterate over history output
// 
//         for( int i = 0; i < historiesVpar[n].size() ; i++ ) { 
//             ofileV << historiesVpar[n][i] << endl ;
//         }
// 
//         ofileV.close();
// 
//         //--------------reacting particles have additional properties
// 
//         if(odtP->Lrxn) {
// 
//             stringstream ssTp;
//             ssTp << n+1 << "_Tpar.dat";
//             string fnameTp = fname + ssTp.str();
// 
//             ofstream  ofileTp(fnameTp.c_str()); 
//             if(!ofileTp) 
//                 *proc.ostrm << "\n\n***************** ERROR OPENING FILE " << fnameTp << endl << endl;
// 
//             //-------------- iterate over history output
// 
//             for( int i = 0; i < historiesTpar[n].size() ; i++ ) { 
//                 ofileTp << historiesTpar[n][i] << endl ;
//             }
// 
//             ofileTp.close();
// 
//             stringstream ssSlip;
//             ssSlip << n+1 << "_Slip.dat";
//             string fnameSlip = fname + ssSlip.str();
// 
//             ofstream  ofileSlip(fnameSlip.c_str()); 
//             if(!ofileSlip) 
//                 *proc.ostrm << "\n\n***************** ERROR OPENING FILE " << fnameSlip << endl << endl;
// 
//             //-------------- iterate over history output
// 
//             for( int i = 0; i < historiesVslip[n].size() ; i++ ) { 
//                 ofileSlip << historiesVslip[n][i] << endl ;
//             }
// 
//             ofileSlip.close();
// 
//             stringstream ssZ;
//             //     ss1 << "Z_" << n+1 << ".dat";
//             ssZ << n+1 << "_Z.dat";
//             string fnameZ = fname + ssZ.str();
// 
//             ofstream  ofileZ(fnameZ.c_str()); 
//             if(!ofileZ) 
//                 *proc.ostrm << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
// 
//             for( int i = 0; i < historiesZmix[n].size() ; i++ ) { 
//                 ofileZ << historiesZmix[n][i] << endl ;
//             }
// 
//             //iterate over history output
//             ofileZ.close();
// 
//             //---------------------
//             //scalar dissipation histories
//             stringstream ssX;
//             ssX << n+1 << "_Chi.dat";
//             string fnameX = fname + ssX.str();
// 
//             ofstream  ofileX(fnameX.c_str()); 
//             if(!ofileX) 
//                 *proc.ostrm << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
// 
//             for( int i = 0; i < historiesChi[n].size() ; i++ ) { 
//                 ofileX << historiesChi[n][i] << endl ;
//             }
// 
//             ofileX.close();
// 
// 	    //---------------------
//  	    //mixture fraction gradient times particle velocity (large Stokes mix frac rate)
//             stringstream ssG;
//             ssG << n+1 << "_GradZ.dat";
//             string fnameG = fname + ssG.str();
// 
//             ofstream  ofileG(fnameG.c_str()); 
//             if(!ofileG) 
//                 *proc.ostrm << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
// 
//             for( int i = 0; i < historiesGradZ[n].size() ; i++ ) { 
//                 ofileG << historiesGradZ[n][i] << endl ;
//             }
// 
//             ofileG.close();
// 
// 
// 	    //---------------------
//  	    //mixture fraction diffusion histories
//             stringstream ssD;
//             ssD << n+1 << "_DiffZ.dat";
//             string fnameD = fname + ssD.str();
// 
//             ofstream  ofileD(fnameD.c_str()); 
//             if(!ofileD) 
//                 *proc.ostrm << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
// 
//             for( int i = 0; i < historiesDiff[n].size() ; i++ ) { 
//                 ofileD << historiesDiff[n][i] << endl ;
//             }
// 
//             ofileD.close();
// 
// 
// 	    //---------------------
//  	    //temperature histories (temperature of gas environment seen by particles)
//             stringstream ssT;
//             ssT << n+1 << "_Tgas.dat";
//             string fnameT = fname + ssT.str();
// 
//             ofstream  ofileT(fnameT.c_str()); 
//             if(!ofileT) 
//                 *proc.ostrm << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
// 
//             for( int i = 0; i < historiesTgas[n].size() ; i++ ) { 
//                 ofileT << historiesTgas[n][i] << endl ;
//             }
// 
//             ofileT.close();
// 	    
//         }
// 
//     }
/*}}}*/
}

///////////////////////////////////////////////////////////////////////////////

/** Set the gas phase mass source and particle mass source
 *  gasMassSource = S_phi (=) kg/m3*s  = ydot * (4/3 *pi*r3) * rho_part * nInPseudoPart.
 * \f[
 *      S_{\phi}(=)\frac{kg}{m^3/s} = \dot{y}_p \rho_p (1-\gamma)
 *      
 * \f]
 * where \f[
 *              \gamma = 1- \frac{m_p}{\rho_p}
 *      \f]
 *
 *  @param cellMassSource   \inout particle mass Source
 */

void particles::setPartAndGasMassSource(vector<double> & gasMassSourceFromEachFakeParticle) {

    /**********************************************************************************************
     * particle model (1) is adapted from Colomba Di Blasi and Carmen Branca,
       Kinetics of Primary Product Formation from Wood Pyrolysis,Ind.Eng.Chem.Res.2001,40,5547-5556
       Particle model (2) is adapted from NUnn et al. [1994]
     **********************************************************************************************/
    if(!odtP->Lprxn)
        return;
    particleMassSource = vector<double>(nPart,0.0);
    set_iyPos();
    
    double Rg = 8.314; //universal gas constant (=J/mol*K)
    if(LuseDiBlasiModel) {
        double k_i = 0.0; // overall rate constant at each cell
        double k_LG_i = 0.0; // gas and liquid rate constant at each cell
        double Y_V_i = 0.0; // volatile yield at each cell
        double m_V_rate_i = 0.0; // mass rate of production of volatiles at each cell
        double Y_i_rate; // rate of change in solid mass fraction at each cell
        double Y_i = 0; // solid mass fraction at each cell
        for (int i = 0; i < nPart; i++) {
            if (!pActive[i])
                continue;
            for (int j = 0; j < radial_ngrd; j++) {
                k_i = 4.4E9 * exp(-141000.0 / (Rg * cellTemp[i][j]));
                k_LG_i = 1.5E10 * exp(-149000.0 / (Rg * cellTemp[i][j]));
                Y_V_i = (k_LG_i) / k_i;

                Y_i = cellDens[i][j] / pDens0[i]; //Assuming constant volume
                Y_i_rate = -k_i * Y_i;
                //dm/dt = -ypdot*rhoP*volP
                cellMassSource[i][j] = Y_i_rate * cellVolume[i][j] * cellDens[i][j]; //(=)kg/s
                m_V_rate_i = -cellMassSource[i][j] * Y_V_i; //(=)kg/s
                particleMassSource[i] += cellMassSource[i][j];
                gasMassSourceFromEachFakeParticle[i] += m_V_rate_i;
            }
        }
    }
    else if(LuseNunnsModel){
        double Y_V_inf = 0.9297; //ultimate Yield
        double Y_V_i = 0; // volatile Yield at each cell
        double Y_V_i_rate = 0; // rate of volatile yield at each cell
        for(int i=0; i<nPart; i++){
            if (!pActive[i])
                continue;
            for (int j = 0; j < radial_ngrd; j++) {
                double densityRatio = cellDens[i][j]/pDens0[i];
                /*If the ratio is less than (1-Y_V_inf), then the Y_V_i is negative. If it is negative
                 make it equal to the ultimate yield so that the Y_V_i_rate equals 0*/
                if(densityRatio > (1-Y_V_inf)){  
                    Y_V_i = 1.0 - densityRatio;
                }
                else{
                    Y_V_i = Y_V_inf;
                }
                Y_V_i_rate = 33884.0 * exp (-69000.0 / (Rg * cellTemp[i][j])) * (Y_V_inf - Y_V_i);
                cellMassSource[i][j] = -Y_V_i_rate * cellVolume[i][j] * cellDens[i][j];
                particleMassSource[i] += cellMassSource[i][j];
                gasMassSourceFromEachFakeParticle[i] += -cellMassSource[i][j];

                // -- Particle moisture.
                double originalMoistureMass = moist_content/(100*radial_ngrd) * (initPartMass_i/radial_ngrd);
                double mass_moisture = originalMoistureMass - moisture_tracker[i][j];
                rho_M[i][j] = mass_moisture/cellVolume[i][j];
                double A_moist = 5.13E10;
                double Ea_moist = 88E3;
                cellMoistureSource[i][j]=0.0; 

                double rateConstant_moist = A_moist * exp(-Ea_moist/(Rg*cellTemp[i][j]));
                if (rho_M[i][j] > 0 && cellTemp[i][j] > (95+273))
                    cellMoistureSource[i][j] = rateConstant_moist * rho_M[i][j] * cellVolume[i][j];
                else
                    cellMoistureSource[i][j] = 0.0;
                

            }
            
        }
    }
    else{
        *proc.ostrm << "Kinetic model for particle combustion is not specified. Program terminated for the reason.";
        exit(1);
    }
}

////////////////////////////////////////////////////////////////////////////
/**
 * sets particle enthalpy
 * @param lambda_gas thermal conductivity of the gas
 * @param Gvel gas velocity because of cell contraction and expansion
 */
void particles::setPartEnthSource(vector<double> & lambda_gas, vector<double> & Gvel){

    //----------------- Compute heat flux
    setHeatFlux(lambda_gas,Gvel);
    for (int i = 0; i < nPart; i++) {
        if (!pActive[i])
            continue;
        double cellVolatileEnth;
        for (int j = 0; j < radial_ngrd; j++){
            cellVolatileEnth = cellEnth[i][j]; //------ Neglecting heat of volatilization 
            cellEnthSource[i][j] = (cellFaceHeatFlux[i][j] * cellArea[i][j] - cellFaceHeatFlux[i][j + 1] * cellArea[i][j + 1] +
                    (-cellMassSource[i][j]) * (cellEnth[i][j] - cellVolatileEnth)); // cellMass[i][j];

            cellEnthSource[i][j] = cellEnthSource[i][j]+heatOfVaporization*cellMoistureSource[i][j];
        }
    }
    
}

///////////////////////////////////////////////////////////////////////////////

/**
 * Computes the heat transfer coefficient for the ith particle.
 *
 * \cond
 *     (h = Nu * lambda_gas / 2 * pRadi[i]) 
 * \endcond
 *
 * \f[
 *     h = \nu \frac{\lambda_{gas}}{2} r_{p,i}
 * \f]
 *
 * @param i          \input Index corresponding the cell
 * @param lambda_gas \input thermal conductivity of the gas at the face value (information is contained
 * in the diffuser)
 * @param gRelVel    \input g Relative Velocity
 * @return heat transfer coefficient corresponding to the ith cell in the domain
 */


double particles::computeHeatTranCoeff(int& i, double &lambda_gas, double &gRelVel) {
   
    double cp_gas = 0;
    double mu_gas = 0;
	double thermalCond_gas = lambda_gas;    

    
    if (odtP->ItableLookup) {
        cp_gas = line->etaTools->luTable->getValAtGridPoint(iyPos[i], line->etaTools->luTable->index.cp_index);
        mu_gas = line->etaTools->luTable->getValAtGridPoint(iyPos[i], line->etaTools->luTable->index.visc_index);
	}
    else {
        line->gas->setState_HP(line->enth[iyPos[i]], line->pres, 1E-10);
        cp_gas = line->gas->cp_mass();
        mu_gas = line->molec[iyPos[i]];
    }

//    double u_gas_square = pow(line->uvel[iyPos[i]],2.0);
//    double v_gas_square = pow(line->vvel[iyPos[i]], 2.0);
//    double w_gas_square = pow(line->wvel[iyPos[i]],2.0);
//    double gas_vel = pow(u_gas_square + v_gas_square + w_gas_square, 1.0 / 2.0);
    
    double u_part_square = pow(uvel[i],2.0);
    double v_part_square = pow(vvel[i], 2.0);
    double w_part_square = pow(wvel[i], 2.0);
    double part_vel = pow(u_part_square+v_part_square+w_part_square,1.0/2.0);
    //added by apaudel to force vHT
    //vHT[i] = 232;
    double relVel  = vHT[i] + gRelVel - part_vel;
        
//    double relVel = gas_vel + gRelVel - part_vel;
    double rho_gas = line->rho[iyPos[i]];
    double Re = rho_gas *relVel * 2 * pRadi[i] / mu_gas;
    double Pr = cp_gas * mu_gas / thermalCond_gas;
    double Nu = 0;
    if(pShape == 2){ //sphere
        /*Ranz and Marshall*/
        Nu = 2 + 0.6 * pow(abs(Re), 0.5) * pow(Pr, 1.0/3.0);
    }
    else if(pShape ==1){
        /*Churchill and Bernstein.*/
        Nu = 
        0.3 + (0.62*pow(abs(Re),0.5)*pow(Pr,1.0/3.0))/(pow(1+pow(0.4/Pr,2.0/3.0),1.0/4.0)) *pow(1+pow(abs(Re)/282000.0,5.0/8.0),4.0/5.0);
    }

    return Nu * thermalCond_gas / (2 * pRadi[i]);
}

///////////////////////////////////////////////////////////////////////////////

/**
 * Computes the blowing factor of particle 'i'
 * @param i \input index of the particle of which blowing factor is to be computed.
 * @return the blowing factor corresponding to the ith cell in the domain
 */
double particles::computeBlowingFactor(int &i, double& h_loc) {
    double bf = 1.0;
        if (line->voidFrac[iyPos[i]] < 1.0) {
            if (abs(particleMassSource[i]) > 0) {
                double phi = (-particleMassSource[i] / cellArea[i][radial_ngrd]) * pCp / h_loc;
                bf = phi / (exp(phi) - 1.0);
            }
        }
    return bf;
}

/////////////////////////////////////////////////////////////////////
/**
 * Computes the thermal conductivity at a given temperature
 * 
 * @param temp
 * @return
 * adapted from Hong Lu
 * \f[
 *      \lambda = =0.1272*(1+0.00205*(temp-273))*0.986     (W/mK)
 * \f]
 * 
 */
double particles::getThermalCond(double& temp){
    return 0.08831;
    // return 0.1272*(1+0.00205*(temp-273.15))*0.986;
}


////////////////////////////////////////////////////////////////////////////

/**
 * Sets the heat flux vector (For particle temperature distribution)
 */
void particles::setHeatFlux(vector<double> & lambda_gas,vector<double>& Gvel){
    
    set_iyPos();
    for (int i = 0; i < nPart; i++) {
        if(!pActive[i])
            continue;
        //boundary Condition dT/dR|r=0 = 0;
        cellFaceHeatFlux[i][0] = 0.0;
        
        //boundary condition dT/dr|r=R = -heatChangeDueToConvection - heatChangeDueToRadiation
        
        cellFaceHeatFlux[i][radial_ngrd] =  -1 * (getHeatFluxDueToConvec(lambda_gas[iyPos[i]],Gvel[iyPos[i]],i) + getHeatFluxDueToRad(i));
        
        //compute the heat flux at other positions q = - k * dT/dr
        for (int j = 1; j < radial_ngrd; j++) {
            cellFaceHeatFlux[i][j] = -getThermalCond(cellTemp[i][j])/cell_size * (cellTemp[i][j] - cellTemp[i][j-1]);
        }
    }
}

////////////////////////////////////////////////////////////////////////////

/**
 * Computes a linearly spaced vector (same functionality as matlab's linspace)
 * @param x1 beginning value
 * @param x2 ending value
 * @param n number of points between x1 and x2
 * @return linearly spaced vector between x1 and x2 with n points
 */

vector<double> particles::linspace(double x1, double x2, int n){
    vector<double> linSpacedVector(floor(n),1);
    if(floor(n)<=0){
        return linSpacedVector;
    }
    if(n==1){
        linSpacedVector[0] = x2;
    }
    else{
        linSpacedVector[0]=x1;
        for(int i=1;i<floor(n);i++){
            double dx  = (x2-x1)/(double(floor(n))-1);
            linSpacedVector[i] = x1+dx*i;
        }
    }
    return linSpacedVector;

}

///////////////////////////////////////////////////////////////////////

/**Computes the heat capacity at a given temperature
 * @author apaudel
 * @param temp \input temperature at which the heat capacity is to be computed
 * @return the heat capacity at the temperature temp
 * 
 * 
 * \f[
 * 
 *      c_{p,\,i} = 7.6658E(-7) * T^(3) - 3.6335E(-3) * T^(2) + 5.9279 * T
 * \f]
 */

double particles::getHeatCapFromTemp(double temp){
    
    if(temp > 2000){ //Above 2000K, the heat capacity has a constant profile
        temp = 2000;
    }
    return 817.69;
    // return 7.6658E-7 * pow(temp,3) - 3.6335E-7 * pow(temp,2) + 5.9279 * temp
    //         -326.96;
    
}

////////////////////////////////////////////////////////////////////////////
/**
 * Returns the surface temperature of the particle based on the interpolation
 * @param i particle number
 * @return the surface temperature of the particle.
 */

double particles::getSurfTemp(int& i){
    if(radial_ngrd > 1)
        return cellTemp[i][radial_ngrd-1] + (cellTemp[i][radial_ngrd-1] - cellTemp[i][radial_ngrd-2])/ 2;
    else
        return cellTemp[i][radial_ngrd-1];
}

////////////////////////////////////////////////////////////////////////////
/**
 * returns the surrounding gas temperature of the particle
 * @param i the particle number of which the surrounding gas temperature is to be found.
 * @return surrounding gas temperature based on lagrange interpolating polynomial
 * \f[
 *      P(x) = \sum_{j=1}^{n}P_{j}(x), where P_{j}(x) = y_{j} \prod_{k=1}^{n}\frac{x-x_{k}}{x_{j}-x_{k}}
 * \f]
 */
double particles::getPartSurrGasTemp(int & i){
    /*Goal: Find y.*/
    double x1 = 0; double x2 = 0; double x3 = 0;
    double y1 = 0; double y2 = 0; double y3 = 0;
    double x  = 0; //double y = 0; //  !!!!!  unused variable
    //Find x
    x = yPos[i];
    
    if(iyPos[i]==0){ //location of the particle is at the beginning
        x1 = line->pos[0];
        y1 = line->temp[0];
        x2 = line->pos[1];
        y2 = line->temp[1];
        x3 = line->pos[2];
        y3 = line->temp[2];
    }
    else if(iyPos[i]==line->ngrd-1){ //location of the particle is at the end
        x1 = line->pos[line->ngrd-1];
        y1 = line->temp[line->ngrd-1];
        x2 = line->pos[line->ngrd-2];
        y2 = line->temp[line->ngrd-2];
        x3 = line->pos[line->ngrd-3];
        y3 = line->temp[line->ngrd-3];
    }
    else { //location of the particle is elsewhere
        if(!pActive[i])
            return line->temp[0];
        x1 = line->pos[iyPos[i]];
        y1 = line->temp[iyPos[i]];
        x2 = line->pos[iyPos[i]-1];
        y2 = line->temp[iyPos[i]-1];
        x3 = line->pos[iyPos[i]+1];
        y3 = line->temp[iyPos[i]+1];
    }
    
    return ( (x-x2) * (x-x3) )/( (x1-x2) * (x1-x3) ) * y1
            + ( (x-x1) * (x-x3) )/( (x2-x1) * (x2-x3) ) * y2
            + ( (x-x1) * (x-x2) )/( (x3-x1) * (x3-x2) ) * y3;
}


////////////////////////////////////////////////////////////////////////////
/**
 * Method to calculate the total enthalpy the particle
 * @return the instantaneous total enthalpy of the particle (J)
 */
double particles::getTotalParticleEnthalpy() {
    double totalPartEnth = 0;
    for (int i = 0; i < nPart; i++) {
        totalPartEnth += getParticleEnthalpy(i);
    }
    return totalPartEnth;
    
}

////////////////////////////////////////////////////////////////////////////

/**
 * Method to calculate the total mass of the particles
 * @return the instantaneous mass of the particles (kg)
 */

double particles::getTotalParticleMass() {
    double totalPartMass = 0;
    for (int i = 0; i < nPart; i++) {
        totalPartMass += getParticleMass(i);
    }
    return totalPartMass;
}

////////////////////////////////////////////////////////////////////////////

/**
 * Method to calculate the total enthalpy the particle i
 * @return the instantaneous total enthalpy of the particle (J)
 */

double particles::getParticleEnthalpy(int &i) {
    return getFakeParticleEnth(i) * nInPseudoPart[i];
}

////////////////////////////////////////////////////////////////////////////

/**
 * Method to calculate the total mass of the particle i
 * @return the instantaneous mass of the particles (kg)
 */

double particles::getParticleMass(int & i) {
    return getFakeParticleMass(i)*nInPseudoPart[i];
}

////////////////////////////////////////////////////////////////////////////

/**
 * calculation of the mass of a fake particle
 * @param i Particle number of whose the mass is to be calculated
 * @return 
 */

double particles::getFakeParticleMass(int& i){
    double fakePartMass = 0;
    for (int j = 0; j < radial_ngrd; j++) {
        fakePartMass += cellMass[i][j];
    }
    return fakePartMass;

}

////////////////////////////////////////////////////////////////////////////

/**
 * Calculation of the enthalpy of a fake particle
 * @param i Particle number of whose the enthalpy is to be calculated
 * @return fake particle's enthalpy (J)
 */

double particles::getFakeParticleEnth(int & i){
    double fakePartEnth = 0;
    for (int j = 0; j < radial_ngrd; j++){
        fakePartEnth += cellMass[i][j] * cellEnth[i][j];
    }
    return fakePartEnth;
}

////////////////////////////////////////////////////////////////////////////

/**
 * Calculation of the intensive enthalpy of a fake particle
 * @param i Particle number of whose the enthalpy is to be calculated
 * @return fake particle's enthalpy J/kg
 */

double particles::getFakeParticleIntensiveEnth(int& i) {
    return getFakeParticleEnth(i) / getFakeParticleMass(i);
}

////////////////////////////////////////////////////////////////////////////

/*
 * Helper function to compute temperature using newton's method
 * Method setPartTemp uses this function
 */

double tempSolver(double Hf, double cellEnth, double temp) {
    double a = 7.6658E-7;
    double b = -3.6335E-3;
    double c = 5.9279;
    double d = -326.96;
    double refTemp = 298.15;
    
    return ( Hf - cellEnth)+(a / 4 * pow(temp, 4.0) + b / 3 * pow(temp, 3.0) + c / 2 * pow(temp, 2.0) + d * temp) -
            (a / 4 * pow(refTemp, 4.0) + b / 3 * pow(refTemp, 3.0) + c / 2 * pow(refTemp, 2.0) + d * refTemp);

}

////////////////////////////////////////////////////////////////////////////

/**
 * Solves particle temperature using newtons method. The integration of temperature dependent Cp
 * is in function tempSolver.
 * @param i particle number
 * @param j grid point number
 * @param guessValue guess value for the temperature
 * @return temperature of ith particle at jth radial point.
 */

double particles::setPartTemp(int i, int j){
     double originalGuessValue = (cellEnth[i][j] - Hf) / getHeatCapFromTemp(cellTemp[i][j]) + 298.15;
     double gVal = originalGuessValue; //guess Value
     double tol = 0.001; //tolerance
     int maxIter = 100; //maximum iteration
     int itrNum = 1; //iteration Number
     double dSize = 0.0005; //step size
     double gNew = gVal; //new guess value
     while (itrNum < maxIter) {
         
         gNew = gVal - tempSolver(Hf, cellEnth[i][j], gVal) /
                 ((tempSolver(Hf, cellEnth[i][j], gVal + dSize) - tempSolver(Hf, cellEnth[i][j], gVal)) / dSize);
         if (abs(gNew - gVal) < tol) {
             if (gNew > peakTemp){
                 peakTemp = gNew;
             }
             return gNew;
         }
         itrNum++;
         gVal = gNew;
     }
     cout << endl << "WARNING:: Maximum iteration reached while setting the particle temperature. " << endl;
     cout << "Guess Value--> " << originalGuessValue << ". End value--> " << gNew << endl;

     return gNew;   

//     return cellEnth[i][j]/getHeatCapFromTemp(cellTemp[i][j]) + 298.15;
}


////////////////////////////////////////////////////////////////////////////

/**
 * calculation of total heat change due to convection
 * @param lambda_gas thermal conductivity of the gas
 * @param Gvel gas velocity due to cell expansion and contraction
 * @param i Particle number of whose the heat change is to be calculated
 * @return heat change due to convection (positive if particle Temp less than cell Temp), negative otherwise
 */

double particles::getHeatRateDueToConvec(double & lambda_gas, 
                                                 double & Gvel, 
                                                 int& i) {
    return getHeatFluxDueToConvec(lambda_gas,Gvel,i)*cellArea[i][radial_ngrd];
}

////////////////////////////////////////////////////////////////////////////

/**
 * calculation of total heat change due to radiation
 * @param i Particle number of whose the heat change is to be calculated
 * @return heat change due to radiation.
 */

double particles::getHeatRateDueToRad(int& i){
    return getHeatFluxDueToRad(i)*cellArea[i][radial_ngrd];
}

////////////////////////////////////////////////////////////////////////////

/**
 * to calculate the heat flux due to convection at the surface of the particle (positive
 * flux if cell temperature > surface temperature of the particle, negative otherwise)
 * @param lambda_gas thermal conductivity of the gas
 * @param Gvel gas velocity due to cell compression and cell expansion
 * @param i particle index
 * @return the heat flux at the surface of the particle due to covection
 */

double particles::getHeatFluxDueToConvec(double& lambda_gas, double& Gvel, int& i){
    //(1) Convection
    
    double gThermalCond = odtP->ItableLookup ?
            line->etaTools->luTable->getValAtGridPoint(iyPos[i], line->etaTools->luTable->index.lambda_index)
            : lambda_gas;
    double hp = computeHeatTranCoeff(i, gThermalCond, Gvel);
//cout << endl << "hp = " << hp << endl;    
    double bf = computeBlowingFactor(i, hp);
// cout << endl << "line->temp = " << line->temp[iyPos[i]] << endl;     
// cout << endl << "Surftemp = " << getSurfTemp(i) << endl;     
// cout << endl << "blowingfactor = " << bf << endl;     
    double heatFluxDueToConvection = hp * bf * (line->temp[iyPos[i]] - getSurfTemp(i));
    return heatFluxDueToConvection;
}

////////////////////////////////////////////////////////////////////////////
/**
 * Calculates the heat flux due to radiation at the surface of the particle
 * @param i particle index
 * @return heat flux due to radiation of particle i
 */

double particles::getHeatFluxDueToRad(int & i){
    if (isnan(abs(diff->radSource_P[i])))
        cout << endl << "Rad source is nan."<< endl;
    return diff->radSource_P[i]/cellArea[i][radial_ngrd];
}

///////////////////////////////////////////////////////////////////////

/** helper function for sort function in particles::updateEddyInfoArray
 *
 *  @Guangyuan Sun 07/2012
 */

bool compare(particles::eddyInformation a, particles::eddyInformation b) {
    return a.endTime < b.endTime;
}

///////////////////////////////////////////////////////////////////////

/**Store the effect of active eddy on ONE particle. This includes                  \n
 * (1)eddy end time (eddyEndTime)                                                  \n
 * (2)relative velocity between particle and eddy (relativeVelEdPart)
 *
 * @author Guangyuan Sun 07/2012
 */

void particles::updateEddyInfoArray(int iPart, double eddyEndTime, double relativeVelEdPart) {

    int size = eddyInfo[iPart].size();
    eddyInfo[iPart].resize(size+1);
    eddyInfo[iPart][size].endTime = eddyEndTime; 
    eddyInfo[iPart][size].relativeVel = relativeVelEdPart;
    sort (eddyInfo[iPart].begin(),eddyInfo[iPart].end(),compare);
}

///////////////////////////////////////////////////////////////////////

/**Check active eddies for each particle at current time
 * remove unactive eddies from the vector
 * @param time \input current input time
 *
 * @Guangyuan Sun 06/2012
 *
 */

void particles::checkActiveEddy(double time) {
    
    for(int iPart=0; iPart<nPart; iPart++) { 
        if(!pActive[iPart]) continue;
        if(eddyInfo[iPart].empty()) 
            continue;    

        // --------- find out how many eddies are dead at current time
        std::vector<eddyInformation>::iterator it = eddyInfo[iPart].begin();
        std::vector<eddyInformation>::iterator end = eddyInfo[iPart].begin();
        for ( it=eddyInfo[iPart].begin() ; it < eddyInfo[iPart].end(); it++) {
            if (it->endTime <= time)
                continue;
            else {
                end = it;
                break;
            }
        }
        eddyInfo[iPart].erase(eddyInfo[iPart].begin(),end);
    }
}

///////////////////////////////////////////////////////////////////////

/**Calculate the sum of relative velocities for each particle.
 * \b Note that if high accurancy (there might be several eddies "dead" during dtDiffuse in diffuseSingelstep), further work is needed
 * @param time \input current input time
 *
 * @author Guangyuan Sun 
 * @date   06/2012
 *
 */

void particles::getRelativeVel(double time) {
    checkActiveEddy(time);

    for(int iPart=0; iPart<nPart; iPart++) {
        if(!pActive[iPart]) continue;

        vHT[iPart] = 0; 
        for(int i=0; i<eddyInfo[iPart].size(); i++) {
            vHT[iPart] += eddyInfo[iPart][i].relativeVel;
        }
    }
}

