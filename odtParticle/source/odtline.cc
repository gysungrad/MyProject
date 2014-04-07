/**
 * @file odtline.cc
 * Source file for class odtline
 */

#include <iomanip>
#include <cmath>          // fabs
#include <sstream>
#include <cstdlib>
#include "processor.h"
#include "odtline.h"

using namespace std;

#ifdef DOCANTERA
#include "equilibrium.h"
#endif

extern processor proc;

///////////////////////////////////////////////////////////////////////////////

/** Overloaded assignment operator.
 *
 *  @param odtl \input odtline to use for assignment.
 */

void odtline::operator=(const odtline &odtl) {

    LhasRxn = odtl.LhasRxn;
    Lprxn   = odtl.Lprxn;
    LhasVel = odtl.LhasVel;
    LhasTemp= odtl.LhasTemp;
    LhasMom = odtl.LhasMom;
    LhasEta = odtl.LhasEta;

    LhasScl = odtl.LhasScl;

    uvel = odtl.uvel;
    vvel = odtl.vvel;
    wvel = odtl.wvel;

    enth    = odtl.enth;
    yspc    = odtl.yspc;
    temp    = odtl.temp;
    mixf    = odtl.mixf;
    chi     = odtl.chi;
    gradZ   = odtl.gradZ;
    diffZ   = odtl.diffZ;
    pres    = odtl.pres;
    gas     = odtl.gas;
    odtP    = odtl.odtP;
    tran    = odtl.tran;
    strm    = odtl.strm;
    spNames = odtl.spNames;
    adptVar = odtl.adptVar;

    mom     = odtl.mom;  
    eta     = odtl.eta;  
    nmom    = odtl.nmom;
    neta    = odtl.neta;

    scl     = odtl.scl;  
    nscl    = odtl.nscl;


    ngrd      = odtl.ngrd;
    ngrdf     = odtl.ngrdf;
    Ldomain   = odtl.Ldomain;
    pos       = odtl.pos;
    posf      = odtl.posf;
    rho       = odtl.rho;
    molec     = odtl.molec;
    lambda    = odtl.lambda;
    phase     = odtl.phase;
    nprops    = odtl.nprops;
    bcprops   = odtl.bcprops;
    propNames = odtl.propNames;
    upBrho    = odtl.upBrho;
    loBrho    = odtl.loBrho;
    upBmolec  = odtl.upBmolec;
    loBmolec  = odtl.loBmolec;
    upBlambda = odtl.upBlambda;
    loBlambda = odtl.loBlambda;
    upBprops  = odtl.upBprops;
    loBprops  = odtl.loBprops;
    upBuptr   = odtl.upBuptr;
    loBuptr   = odtl.loBuptr;

    nspc        = odtl.nspc;
    LpropsHasYi = odtl.LpropsHasYi;
    iPropsYi    = odtl.iPropsYi;
    etaTools    = odtl.etaTools;
    momTools    = odtl.momTools;
    if(odtl.odtP->Iparticles)
        voidFrac = odtl.voidFrac;

    for(int k=0; k<nprops; k++) 
        *props[k] = *(odtl.props[k]);

    uptr        = &uvel;
    heatLoss    = odtl.heatLoss;
    lambda      = odtl.lambda;

    setMesher(true);

}

///////////////////////////////////////////////////////////////////////////////

/**Constructor function.
 * \b Note: this function sets properties for itself and also for its
 * base class #anyline \n
 * \b Note: mixf is not consistent with species and enthalpy/temp here. \n
 * Just creating/sizing objects, not setting them per se.  
 * Use setOdtline for that.
 *
 * @param npts          \input initial number of points on the line.
 * @param Ld            \input initial line length.
 * @param cantIG        \input pointer to an #IdealGasMix object to assign member pointer.
 * @param cantTran      \input pointer to a #Transport object to assign member pointer.
 * @param strm_p        \input pointer to a #streams object to assign member pointer.
 * @param odtP_p        \input pointer to a #odtParam object to assign member pointer.
 * @param caseSetupFile \input name of case input file to read.
 * @param LhasVel_p     \input flag indicates u,v,w velocities are used
 * @param LhasRxn_p     \input flag indicates h,Yi are used
 */

odtline::odtline(int npts, double Ld, IdealGasMix *cantIG, 
                                      Transport   *cantTran, 
                                      streams     *strm_p,
                                      odtParam    *odtP_p,
                                      string caseSetupFile,
                                      bool LhasVel_p, bool LhasRxn_p) : anyline(npts,Ld) {


    pres=-1;

    odtP    = odtP_p;

    LhasVel = LhasVel_p;
    LhasRxn = LhasRxn_p;
    LhasTemp= odtP->LhasTemp;
    LhasMom = odtP->Imom>0 ? true : false; 
    LhasEta = odtP->Ieta>0 ? true : false; 
    Lprxn   = odtP->Lprxn;
    LhasScl = odtP->Iscl>0 ? true : false; 	

    if(!LhasVel && !LhasRxn) {
        cerr << endl << "Error, LhasVel and LhasRxn both false" << endl;
        //exit(0); ///Alan commented this on 10-13-13: why not nonreacting LEM?
    }

    ngrd    = npts;
    ngrdf   = ngrd+1;
    Ldomain = Ld;
    pos     = vector<double>(ngrd,  0.0);
    posf    = vector<double>(ngrdf, 0.0);

    gas     = cantIG;
    tran    = cantTran;
    strm    = strm_p;

    nspc    = (LhasRxn || odtP->ItableLookup) ? gas->nSpecies() : 0;
    spNames = (LhasRxn || odtP->ItableLookup) ? gas->speciesNames() : vector<string>(0);
    pres    = strm->pres;

    nmom    = LhasMom ? odtP->Imom : 0;
    neta    = LhasEta ? odtP->Ieta : 0;

    nscl    = LhasScl ? odtP->Iscl : 0;

    iPropsYi = -1;
    if(LhasRxn) { iPropsYi = LhasVel ? 4 : 1; }

    //----------- set enthalpy, density, temperature, mixf

    rho     = vector<double>(ngrd, odtP->rho_0);
    molec   = vector<double>(ngrd, odtP->visc_0); 
    phase   = vector<double>(ngrd, odtP->phase);
    if(LhasTemp) {
        lambda = vector<double>(ngrd, odtP->lambda_0);
    }
    else
        lambda = vector<double>(ngrd, 0.0);

    if(LhasVel) {
        uvel = vector<double>(ngrd,0.0);
        vvel = vector<double>(ngrd,0.0);
        wvel = vector<double>(ngrd,0.0);
        
        uptr = &uvel;
    }
    
    if(LhasTemp) {
        temp    = vector<double>(ngrd, 0.0);
    }
    
    //----------- set default species mass fractions
    if(LhasRxn || Lprxn || odtP->ItableLookup) {
        temp    = vector<double>(ngrd, 300.);
        enth    = vector<double>(ngrd);
        mixf    = vector<double>(ngrd,  0.0);
        chi     = vector<double>(ngrd,  0.0);
        gradZ   = vector<double>(ngrd,  0.0);
        diffZ   = vector<double>(ngrd,  0.0);

        yspc.resize(nspc);
        for(int k=0; k<nspc; k++)
            yspc[k] = vector<double>(ngrd, 0.0);
    }

    //----------- set particle moments

    if(LhasMom) {
        mom.resize(nmom);
        for(int k=0; k<nmom; k++)
            mom[k] = vector<double>(ngrd, 0.0);
    }

    //----------- set eta parameters

    if(LhasEta) {
        eta.resize(neta);
        for(int k=0; k<neta; k++)
            eta[k] = vector<double>(ngrd, 0.0);
    }

    //----------- set scl parameters

    if(LhasScl) {
        scl.resize(nscl);
        for(int k=0; k<nscl; k++)
            scl[k] = vector<double>(ngrd, 0.0);
    }

    //---------- set props array pointers
    
    nprops = 0;
    if(LhasVel) nprops += 3;
    if(LhasTemp) nprops++;
    if(LhasRxn) nprops += 1+nspc;
    if(LhasEta) nprops += neta;
    if(LhasMom) nprops += nmom;

    if(LhasScl) nprops += nscl;

    props.resize(nprops);

    int idmb = 0;
    if(LhasVel) {
        props[0] = &uvel;
        props[1] = &vvel; 
        props[2] = &wvel;
        idmb=3;
    }
    if(LhasTemp) {
        props[0+idmb] = &temp;
        idmb++;
    }
    if(LhasRxn) {
        props[0+idmb] = &enth;
        idmb++;
        for(int k=0; k<nspc; k++)
            props[k+idmb] = &yspc[k];
        idmb += nspc;
    }
    if(LhasEta) {
        for(int k=0; k<neta; k++)
            props[k+idmb] = &eta[k];
        idmb += neta;
    }
    if(LhasMom) {
        for(int k=0; k<nmom; k++)
            props[k+idmb] = &mom[k];
        idmb += nmom;
    }

    if(LhasScl) {
        for(int k=0; k<nscl; k++)
            props[k+idmb] = &scl[k];
        idmb += nscl;
    }

    //---------- set props names

    propNames = vector<string>(nprops, "");
    idmb=0;
    if(LhasVel) {
        propNames[0] = "uvel";
        propNames[1] = "vvel";
        propNames[2] = "wvel";
        idmb=3;
    }
    if(LhasTemp) {
        propNames[idmb]   = "temp";
        idmb++;
    }
    if(LhasRxn) {
        propNames[0+idmb] = "enth";
        idmb++;
        for(int k=0; k<nspc; k++)
            propNames[k+idmb] = "Y_" + spNames[k];
        idmb += nspc;
    }
    string       s1;
    stringstream ss1;
    if(LhasEta) {
        for(int k=0; k<neta; k++) {
            ss1.clear(); ss1 << k; ss1 >> s1;
            propNames[k+idmb] = "Eta_" + s1;
        }
        idmb += neta;
    }
    if(LhasMom) {
        for(int k=0; k<nmom; k++) {
            ss1.clear(); ss1 << k; ss1 >> s1;
            propNames[k+idmb] = "Mom_" + s1;
        }
        idmb += nmom;
    }

    if(LhasScl) {
        for(int k=0; k<nscl; k++) {
            ss1.clear(); ss1 << k; ss1 >> s1;
            propNames[k+idmb] = "Scl_" + s1;
        }
        idmb += nscl;
    }

    //---------- set prop bounds

    upBrho    = 1.0E20;
    loBrho    = 0.0;
    upBmolec  = 1.0E20;
    loBmolec  = 0.0;
    upBlambda = 1.0E20;
    loBlambda = 0.0;
    upBuptr   = 1.0E20;
    loBuptr   = -1.0E20;
    upBprops  = vector<double>(nprops, 1.0E20);
    loBprops  = vector<double>(nprops, -1.0E20);
    if(LhasRxn) {
        for(int k=iPropsYi; k<nprops; k++) {  // Yi
            upBprops[k] = 1.0E20;      // used in anyline::splitCell for rho*Yi (props is Yi)
            loBprops[k] = 0.0;
        }
    }
    
    //---------- set props boundary conditions
    
    this->setBCprops();
    
    //---------- set cell positions
    
    double dp = Ldomain/ngrd;
    for(int i=1; i<ngrdf; i++)
        posf[i] = posf[i-1] + dp;
    pos[0] = dp/2;
    for(int i=1; i<ngrd; i++)
        pos[i] = pos[i-1] + dp;
    if(odtP->IetaType == odtP->IETATYPE_TABLELOOKUP){
        heatLoss = vector<double>(ngrd,0.0);
        lambda   = vector<double>(ngrd,0.0);
    }

    //---------- setup the mesh adapter variables

    setMesher(true);

    if (odtP->Iparticles)
        voidFrac.resize(ngrd);
    
    if(caseSetupFile=="")
        return;
    
    if(caseSetupFile=="")
      return;
   
    if(LhasRxn || odtP->IetaType==odtP->IETATYPE_TABLELOOKUP || odtP->IetaType==odtP->IETATYPE_DEFAULT || odtP->IetaType==odtP->IETATYPE_AQUEOUS || LhasVel)
        setOdtline(caseSetupFile);    //DH

}

///////////////////////////////////////////////////////////////////////////////

/**Copy constructor.  \b Note: this function sets properties for itself and its
 * objects inherited from base class anyline.
 * 
 * @param odtl \input odtline object to copy.
 */

odtline::odtline(const odtline &odtl) {

    ngrd        = odtl.ngrd;   
    ngrdf       = odtl.ngrdf; 
    Ldomain     = odtl.Ldomain;
    pos         = odtl.pos;   
    posf        = odtl.posf;  
    rho         = odtl.rho;   
    molec       = odtl.molec;
    lambda      = odtl.lambda;
    phase       = odtl.phase;
    nprops      = odtl.nprops;
    bcprops     = odtl.bcprops;
    propNames   = odtl.propNames;
    nspc        = odtl.nspc;
    LpropsHasYi = odtl.LpropsHasYi;
    iPropsYi    = odtl.iPropsYi;
    LhasVel     = odtl.LhasVel;
    LhasTemp    = odtl.LhasTemp;
    LhasRxn     = odtl.LhasRxn;

    LhasEta     = odtl.LhasEta;
    LhasMom     = odtl.LhasMom;
    Lprxn       = odtl.Lprxn;
    nmom        = odtl.nmom;
    neta        = odtl.neta;
    eta         = odtl.eta;
    mom         = odtl.mom;

    LhasScl     = odtl.LhasScl;
    nscl        = odtl.nscl;
    scl         = odtl.scl;

    uvel    = odtl.uvel;
    vvel    = odtl.vvel;
    wvel    = odtl.wvel;
    uptr    = &uvel;
    enth    = odtl.enth;
    pres    = odtl.pres;
    spNames = odtl.spNames;
    temp    = odtl.temp;
    mixf    = odtl.mixf;
    chi     = odtl.chi;  
    gradZ   = odtl.gradZ;  
    diffZ   = odtl.diffZ;  

    gas     = odtl.gas;                            // note: points to same object
    tran    = odtl.tran;                           // note: points to same object
    strm    = odtl.strm;                           // note: points to same object
    odtP    = odtl.odtP;                           // note: points to same object

    adptVar = odtl.adptVar;

    yspc.resize(nspc);
    for(int k=0; k<nspc; k++)
        yspc[k] = odtl.yspc[k];

    if(LhasEta) {
        eta.resize(neta);
        for(int k=0; k<neta; k++)
            eta[k] = odtl.eta[k];
    }

    if(LhasMom) {
        mom.resize(nmom);
        for(int k=0; k<nmom; k++)
            mom[k] = odtl.mom[k];
    }

    if(LhasScl) {
        scl.resize(nscl);
        for(int k=0; k<nscl; k++)
            scl[k] = odtl.scl[k];
    }

    props.resize(nprops);

    int idmb = 0;
    if(LhasVel) {
        props[0] = &uvel;
        props[1] = &vvel; 
        props[2] = &wvel;
        idmb=3;
    }
    if(LhasTemp){
        props[0+idmb] = &temp;
        idmb++;
    }
    if(LhasRxn) {
        props[0+idmb] = &enth;
        idmb++;
        for(int k=0; k<nspc; k++)
            props[k+idmb] = &yspc[k];
        idmb += nspc;
    }
    if(LhasEta) {
        for(int k=0; k<neta; k++)
            props[k+idmb] = &eta[k];
        idmb += neta;
    }
    if(LhasMom) {
        for(int k=0; k<nmom; k++)
            props[k+idmb] = &mom[k];
        idmb += nmom;
    }

    if(LhasScl) {
        for(int k=0; k<nscl; k++)
            props[k+idmb] = &scl[k];
        idmb += nscl;
    }

    upBprops  = odtl.upBprops;
    loBprops  = odtl.loBprops;
    upBrho    = odtl.upBrho;
    loBrho    = odtl.loBrho;
    upBmolec  = odtl.upBmolec;
    loBmolec  = odtl.loBmolec;
    upBlambda = odtl.upBlambda;
    loBlambda = odtl.loBlambda;
    upBuptr   = odtl.upBuptr;
    loBuptr   = odtl.loBuptr;

    if (odtl.odtP->Iparticles)
        voidFrac = odtl.voidFrac;

    etaTools = odtl.etaTools;
    momTools = odtl.momTools;
    heatLoss = odtl.heatLoss;
    lambda   = odtl.lambda;

    //---------- setup the mesh adapter variables

    setMesher(true);


}

///////////////////////////////////////////////////////////////////////////////

/**Constructor builds a line from a section (or sections) of an existing line).
 *
 * @param odtl  \input odtline object to copy.
 * @param i1    \input beginning index of the construction.
 * @param i2    \input ending index of the construction.
 * @param Lwrap \input true if the region wraps around the domain (periodic).
 * @param Ladpt \input false if no mesh adaption is wanted (e.g. for eddies)
 *
 * Used to create eddy lines from existing lines.  The Lwrap flag handles
 * the periodic case of an eddy that wraps around the domain.  See below.
 */

odtline::odtline(const odtline *odtl, int &i1, int &i2, bool Lwrap, bool Ladpt) {
    nprops    = 0;      
    this->initFromCopy(odtl, i1, i2, Lwrap, Ladpt);
}

void odtline::initFromCopy(const odtline *odtl, 
            int &i1, int &i2, bool Lwrap, bool Ladpt) {
    // Note Lwrap is optional and defaults to false
    
    upBprops  = odtl->upBprops;
    loBprops  = odtl->loBprops;
    upBrho    = odtl->upBrho;
    loBrho    = odtl->loBrho;
    upBmolec  = odtl->upBmolec;
    loBmolec  = odtl->loBmolec;
    upBlambda = odtl->upBlambda;
    loBlambda = odtl->loBlambda;
    upBuptr   = odtl->upBuptr;
    loBuptr   = odtl->loBuptr;
    
    pres      = odtl->pres;
    
    gas       = odtl->gas;              // note, same object
    tran      = odtl->tran;             // note, same object
    strm      = odtl->strm;             // note, same object
    odtP      = odtl->odtP;             // note, same object
    nspc      = odtl->nspc;
    
    nmom      = odtl->nmom;
    neta      = odtl->neta;
    
    nscl      = odtl->nscl;
    
    LpropsHasYi = odtl->LpropsHasYi;
    iPropsYi    = odtl->iPropsYi;
    LhasVel     = odtl->LhasVel;
    LhasTemp    = odtl->LhasTemp;
    LhasRxn     = odtl->LhasRxn;
    Lprxn       = odtl->Lprxn;
    LhasEta     = odtl->LhasEta;
    LhasMom     = odtl->LhasMom;
    
    LhasScl     = odtl->LhasScl;
    
    
    spNames = odtl->spNames;
    
    yspc.resize(nspc);
    eta.resize(neta);
    mom.resize(nmom);
    
    scl.resize(nscl);
    
    etaTools   = odtl->etaTools;
    momTools   = odtl->momTools;
    voidFrac   = odtl->voidFrac;
    heatLoss   = odtl->heatLoss;
    lambda     = odtl->lambda;
    
    probType = odtl->probType;
    
    /**\vc{
     * ---------------------------------------------------------------------------------
     *
     * nonwrap: |   | * | * | * | * | * | * |   |   |
     *               i1                  i2
     * new line consists of *'d cells
     * }
     */
    
    
    if(!Lwrap) {          // new vector filled with elements it1 to it2-1
        
        int i2b = i2+1;   // for the vector constructor
        
        ngrd    = i2-i1+1;
        ngrdf   = ngrd+1;
        Ldomain = odtl->posf[i2+1] - odtl->posf[i1]; 
        
        pos.assign(&odtl->pos[i1],   &odtl->pos[i2b]);
        posf.assign(&odtl->posf[i1],  &odtl->posf[i2b]);
        posf.push_back(odtl->posf[i2+1]);
        rho.assign(&odtl->rho[i1],   &odtl->rho[i2b]);
        molec.assign(&odtl->molec[i1], &odtl->molec[i2b]);
        phase.assign(&odtl->phase[i1], &odtl->phase[i2b]);
        lambda.assign(&odtl->lambda[i1], &odtl->lambda[i2b]);
        if (odtP->Iparticles)
            voidFrac.assign(&odtl->voidFrac[i1], &odtl->voidFrac[i2b]);
        if(LhasRxn) {
            enth.assign(&odtl->enth[i1],  &odtl->enth[i2b]);
            temp.assign(&odtl->temp[i1],  &odtl->temp[i2b]);
            mixf.assign(&odtl->mixf[i1],  &odtl->mixf[i2b]);
            chi.assign(&odtl->chi[i1],  &odtl->chi[i2b]);
            gradZ = vector<double> (odtl->gradZ.begin()+i1, odtl->gradZ.begin()+i2b);
            diffZ = vector<double> (odtl->diffZ.begin()+i1, odtl->diffZ.begin()+i2b);
            for(int k=0; k<nspc; k++)
                yspc[k].assign(&odtl->yspc[k][i1], &odtl->yspc[k][i2b]);
        }
        if(LhasVel) {
            uvel.assign(&odtl->uvel[i1],  &odtl->uvel[i2b]);
            vvel.assign(&odtl->vvel[i1],  &odtl->vvel[i2b]);
            wvel.assign(&odtl->wvel[i1],  &odtl->wvel[i2b]);
        }
        if(LhasTemp) {
            temp.assign(&odtl->temp[i1],  &odtl->temp[i2b]);
        }
        if(LhasEta) {
            for(int k=0; k<neta; k++)
                eta[k].assign(&odtl->eta[k][i1], &odtl->eta[k][i2b]);
        }
        if(LhasMom) {
            for(int k=0; k<nmom; k++)
                mom[k].assign(&odtl->mom[k][i1], &odtl->mom[k][i2b]);
        }
        
        if(LhasScl) {
            for(int k=0; k<nscl; k++)
                scl[k].assign(&odtl->scl[k][i1], &odtl->scl[k][i2b]);
        }
    }
    
    /**
     *
     * \vc{---------------------------------------------------------------------------------
     *
     * Wrap: | 4 | 5 |   |   |   |   | 1 | 2 | 3 |
     *            i2                  i1
     * New line consists of #'d cells: 1 2 3 4 5}
     * Meant for periodic domains with eddies that wrap around
     * This creates the eddy, which is then triplet mapped, then inserted.
     * Apply the jump periodic conditions here.\n
     * \b Note: the position of the eddy in this case becomes:
     * \vc{
     *           | 1 | 2 | 3 | 4 | 5 |
     * }
     *  That is, it extends beyond the original region.
     *
     */
    
    else {               // new vectors filled with elements from [i2:ngrd-1, 0:i1]
        
        //int i1b = i1+1;  // for the insert function  //  !!!!!  currently unused variable
        int i2b = i2+1;  // for the insert function
        
        ngrd    = (odtl->ngrd-i1) + (i2+1);
        ngrdf   = ngrd+1; 
        Ldomain = (odtl->posf[odtl->ngrd] - odtl->posf[i1]) + (odtl->posf[i2b]
                    - odtl->posf[0]);
        
        pos.assign(odtl->pos.begin()+i1,   odtl->pos.end());
        int idmb = pos.size();
        pos.insert(  pos.end(),   odtl->pos.begin(),   odtl->pos.begin()+i2b   ); 
        for(int i=idmb; i<(int)pos.size(); i++)
            pos[i]+=odtl->Ldomain;
        
        posf.assign(odtl->posf.begin()+i1,  odtl->posf.end()-1);
        idmb  = posf.size();
        posf.insert( posf.end(),  odtl->posf.begin(),  odtl->posf.begin()+i2b+1);
        for(int i=idmb; i<(int)posf.size(); i++)
            posf[i]+=odtl->Ldomain;
        
        rho.assign(odtl->rho.begin()+i1,   odtl->rho.end());
        rho.insert(  rho.end(),   odtl->rho.begin(),   odtl->rho.begin()+i2b   );
        
        molec.assign(odtl->molec.begin()+i1, odtl->molec.end());
        molec.insert(molec.end(), odtl->molec.begin(), odtl->molec.begin()+i2b );
        
        phase.assign(odtl->phase.begin()+i1, odtl->phase.end());
        phase.insert(phase.end(), odtl->phase.begin(), odtl->phase.begin()+i2b );
        
        lambda.assign(odtl->lambda.begin()+i1, odtl->lambda.end());
        lambda.insert(lambda.end(), odtl->lambda.begin(), odtl->lambda.begin()+i2b );
        
        if (odtP->Iparticles) {
            voidFrac.assign(odtl->voidFrac.begin() + i1, odtl->voidFrac.end());
            voidFrac.insert(voidFrac.end(), odtl->voidFrac.begin(), odtl->voidFrac.begin() + i2b);
        }
        
        if(LhasRxn) {
            
            temp.assign(odtl->temp.begin()+i1, odtl->temp.end());
            temp.insert(temp.end(), odtl->temp.begin(), odtl->temp.begin()+i2b );
            
            mixf.assign(odtl->mixf.begin()+i1, odtl->mixf.end());
            mixf.insert(mixf.end(), odtl->mixf.begin(), odtl->mixf.begin()+i2b );
            
            chi.assign(odtl->chi.begin()+i1, odtl->chi.end());
            chi.insert(chi.end(), odtl->chi.begin(), odtl->chi.begin()+i2b );

            gradZ = vector<double>(odtl->gradZ.begin()+i1, odtl->gradZ.end());
            gradZ.insert(gradZ.end(), odtl->gradZ.begin(), odtl->gradZ.begin()+i2b );

            diffZ = vector<double>(odtl->diffZ.begin()+i1, odtl->diffZ.end());
            diffZ.insert(diffZ.end(), odtl->diffZ.begin(), odtl->diffZ.begin()+i2b );
            
            enth.assign(odtl->enth.begin()+i1, odtl->enth.end());
            enth.insert(enth.end(), odtl->enth.begin(), odtl->enth.begin()+i2b );
            
            for(int k=0; k<nspc; k++) {
                yspc[k].assign(odtl->yspc[k].begin()+i1, odtl->yspc[k].end());
                yspc[k].insert(yspc[k].end(), odtl->yspc[k].begin(), odtl->yspc[k].begin()+i2b );
            }
        }
        
        if(LhasVel) {
            uvel.assign(odtl->uvel.begin()+i1, odtl->uvel.end());
            uvel.insert(uvel.end(), odtl->uvel.begin(), odtl->uvel.begin()+i2b );
            
            vvel.assign(odtl->vvel.begin()+i1, odtl->vvel.end());
            vvel.insert(vvel.end(), odtl->vvel.begin(), odtl->vvel.begin()+i2b );
            
            wvel.assign(odtl->wvel.begin()+i1, odtl->wvel.end());
            wvel.insert(wvel.end(), odtl->wvel.begin(), odtl->wvel.begin()+i2b );
        }
        
        if(LhasTemp) {
            temp.assign(odtl->temp.begin()+i1, odtl->temp.end());
            temp.insert(temp.end(), odtl->temp.begin(), odtl->temp.begin()+i2b );
        }
        if(LhasEta) {
            for(int k=0; k<neta; k++) {
                eta[k].assign(odtl->eta[k].begin()+i1, odtl->eta[k].end());
                eta[k].insert(eta[k].end(), odtl->eta[k].begin(), odtl->eta[k].begin()+i2b );
            }
        }
        
        if(LhasMom) {
            for(int k=0; k<nmom; k++) {
                mom[k].assign(odtl->mom[k].begin()+i1, odtl->mom[k].end());
                mom[k].insert(mom[k].end(), odtl->mom[k].begin(), odtl->mom[k].begin()+i2b );
            }
        }
        
        if(LhasScl) {
            for(int k=0; k<nscl; k++) {
                scl[k].assign(odtl->scl[k].begin()+i1, odtl->scl[k].end());
                scl[k].insert(scl[k].end(), odtl->scl[k].begin(), odtl->scl[k].begin()+i2b );
            }
        }
        
    }
    
    if (nprops != odtl->nprops) {
    // assume that this will not change later
    propNames = odtl->propNames;
    nprops    = odtl->nprops;         
    bcprops   = odtl->bcprops;
    
    props.resize(nprops);
    
    int idmb = 0;
    if(LhasVel) {
        props[0] = &uvel;
        props[1] = &vvel; 
        props[2] = &wvel;
        idmb=3;
    }
    if(LhasTemp) {
        props[0+idmb] = &temp;
        idmb++;
    }
    if(LhasRxn) {
        props[0+idmb] = &enth;
        idmb++;
        for(int k=0; k<nspc; k++)
            props[k+idmb] = &yspc[k];
        idmb += nspc;
    }
    if(LhasEta) {
        for(int k=0; k<neta; k++)
            props[k+idmb] = &eta[k];
        idmb += neta;
    }
    if(LhasMom) {
        for(int k=0; k<nmom; k++)
            props[k+idmb] = &mom[k];
        idmb += nmom;
    }
    
    if(LhasScl) {
        for(int k=0; k<nscl; k++)
            props[k+idmb] = &scl[k];
        idmb += nscl;
    }
    
    // Set BCs for props
    this->setBCprops();
   
    uptr = &uvel;
    }
    
    //---------- setup the mesh adapter variables
    
    setMesher(Ladpt);
}


///////////////////////////////////////////////////////////////////////////////
/** This function sets the boundary conditions depending on the input parameter
 *  bcType of the variable odtP.
 *  
 */
void odtline::setBCprops() {
    
    bcprops = vector<vector<double> >(nprops, vector<double>(8,-10.0));
    // for each prob:   0     1       2       3       4       5      6       7
    //                bcl   bcla    bclb    bclc     bcu    bcua   bcub    bcuc
    // bcl  =   0.0 => periodic     | Dirichlet: bcla =  f(x_bcl)
    //      =   1.0 => Dirichlet    | Neumann:   bcla =  df(x_bcl)/dx
    //      =   2.0 => Neumann      | Cauchy:    bcla =  f(x_bcl) 
    //      =   3.0 => Cauchy       |            bclb =  df(x_bcl)/dx
    //      =   4.0 => Robin        | Robin:     bclc =  bcla * f(x_bcl)
    //      = -10.0 => not set      |                   +bclb * df(x_bcl)/dx
    
    if(odtP->bcType == 0){ // Dirichlet
        for(int i=0; i<nprops; i++){
            bcprops[i][0] = 1.0; bcprops[i][4] = 1.0;
            bcprops[i][1] = 0.0; bcprops[i][5] = 0.0; // currently wall; need to be read
            bcprops[i][2] = 0.0; bcprops[i][6] = 0.0;
            bcprops[i][3] = 0.0; bcprops[i][7] = 0.0;
        }
    }
    else if(odtP->bcType == 1){ // periodic
        for(int i=0; i<nprops; i++){
            bcprops[i][0] = 0.0; bcprops[i][4] = 0.0;
            bcprops[i][1] = 0.0; bcprops[i][5] = 0.0;
            bcprops[i][2] = 0.0; bcprops[i][6] = 0.0;
            bcprops[i][3] = 0.0; bcprops[i][7] = 0.0;
        }
    }
    else if(odtP->bcType == 2){ // wall
        for(int i=0; i<nprops; i++){
            bcprops[i][0] = 1.0; bcprops[i][4] = 1.0;
            bcprops[i][1] = 0.0; bcprops[i][5] = 0.0;
            bcprops[i][2] = 0.0; bcprops[i][6] = 0.0;
            bcprops[i][3] = 0.0; bcprops[i][7] = 0.0;
        }
    }
    else if(odtP->bcType == 3){ // outflow
        for(int i=0; i<nprops; i++){
            bcprops[i][0] = 2.0; bcprops[i][4] = 2.0;
            bcprops[i][1] = 0.0; bcprops[i][5] = 0.0;
            bcprops[i][2] = 0.0; bcprops[i][6] = 0.0;
            bcprops[i][3] = 0.0; bcprops[i][7] = 0.0;
        }
    }
    else if(odtP->bcType == 4){ // inlet
        // TBS
    }
    else if(odtP->bcType == 5){ // wall outflow
        // TBS
    }
    else if(odtP->bcType == 70){ // wall heated channel
        for(int i=0; i<4; i++){
            bcprops[i][0] = 1.0; bcprops[i][4] = 1.0;
            bcprops[i][1] = 0.0; bcprops[i][5] = 0.0;
            bcprops[i][2] = 0.0; bcprops[i][6] = 0.0;
            bcprops[i][3] = 0.0; bcprops[i][7] = 0.0;
        }

        //// not needed boundary conditions for lambda
        //bcprops[4][0] = 2.0; bcprops[4][4] = 2.0;
        //bcprops[4][1] = 0.0; bcprops[4][5] = 0.0;
        //bcprops[4][2] = 0.0; bcprops[4][6] = 0.0;
        //bcprops[4][3] = 0.0; bcprops[4][7] = 0.0;
    }
    else{
        cout << endl << "ERROR:";
        cout << endl << "Unknown boundary type. bcType = " << odtP->bcType;
        cout << endl;
        exit(0);
    }
    
    
    // function to read boundary conditions
    
    
}


///////////////////////////////////////////////////////////////////////////////

/**Output the anyline properties.  
 * Redefined from anyline base class.
 *
 * @param fname \input output file name.
 */

void odtline::outputProperties(string fname) {

    setMixfVec();
    setTempVec();
    setChiVec();
    setDiffZVec(); 
    
    if(odtP->ItableLookup) etaTools->updateOdtLineVecs();
    
    string       s1;
    stringstream ss1;
    int idmb = 0;
    
    ofstream ofile(fname.c_str()); 
    if(!ofile) 
        *proc.ostrm << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
    
    ofile << "# grid points = " << ngrd;
    ofile << "\n# Domain Size = " << Ldomain;
    if(fabs(posf[ngrd] -posf[0] - Ldomain) > 1.0E-6) {
        ofile << "\n# last posf-first posf != Ldomain, last posf = " << posf[ngrd];
        *proc.ostrm << "\n # WARNINGWARNINGWARNING file written with last last face != Ldomain" << endl;
    }
    ofile << "\n# pressure (Pa) = " << pres;
    ofile << endl;
    ofile << setw(19) << "#  1_pos"
          << setw(19) << "2_posf"
          << setw(19) << "3_rho"
          << setw(19) << "4_molec"
          << setw(19) << "5_phase";
    if(LhasTemp) {
        ofile << setw(19) << "6_lambda";
        idmb = 1;
    }
    if(LhasRxn)  {
        ofile << setw(19) << "6_mixf"
              << setw(19) << "7_chi"
              << setw(19) << "8_temp"
             << setw(19) << "9_gradZ"
             << setw(19) << "10_diffZ";
        idmb = 3;
    }
    
    if(odtP->ItableLookup & !LhasRxn){
        ofile<< setw(19) << "6_temp";
        ofile << setw(19) << "7_heatLoss";
        idmb = 2;
    }
    
    for(int k=0; k<nprops; k++) {
        ss1.clear();
        ss1 << k+6+idmb << "_" << propNames[k];
        ss1 >> s1;
        ofile << setw(19) << s1;
    }
    
    
    if(probType==4) //for opposedJets add Pr, Le_s and Sc_s number
    {
        // Pr
        ss1.clear(); ss1 << (6+idmb+nprops) << "_" << "Pr";
        ss1 >> s1;       ofile << setw(19) << s1;
        
        // Le_s
        for(int k=0; k<nspc; k++) {
            ss1.clear(); ss1 << (6+idmb+nprops)+1+k << "_" << "Le_"+spNames.at(k);
            ss1 >> s1;
            ofile << setw(19) << s1;
        }
        
        // Sc_s
        for(int k=0; k<nspc; k++) {
            ss1.clear(); ss1 << (6+idmb+nprops)+1+nspc+k << "_" << "Sc_"+spNames.at(k);
            ss1 >> s1;
            ofile << setw(19) << s1;
        }
    }
    
    
    ofile << scientific;
    ofile << setprecision(10);
    
    for(int i=0; i<ngrd; i++) {
        ofile << endl;
        ofile << setw(19) << pos[i] 
              << setw(19) << posf[i]
              << setw(19) << rho[i]
              << setw(19) << molec[i]
              << setw(19) << phase[i];
        if(LhasTemp) // are part of the props-vector
            ofile << setw(19) << lambda[i];
        if(LhasRxn) 
            ofile << setw(19) << mixf[i]
                  << setw(19) << chi[i]
                  << setw(19) << temp[i]
                   << setw(19) << gradZ[i]
                   << setw(19) << diffZ[i];
        
        if (odtP->ItableLookup && !LhasRxn){
            ofile << setw(19) << temp[i];
            ofile << setw(19) << heatLoss[i];
        }
        
        for(int k=0; k<nprops; k++) 
            ofile << setw(19) << (*props[k])[i];
        
        if(probType==4) //for opposedJets add Pr, Le_s and Sc_s number
        {
            //check size of vectors
            if((int)Pr.size()   != ngrd) Pr.resize(ngrd);
            if((int)Le_s.size() != nspc) Le_s.resize(nspc);
            if((int)Sc_s.size() != nspc) Sc_s.resize(nspc);
            
            // Pr
            ofile << setw(19) << Pr.at(i);
            
            // Le_s
            for(int k=0; k<nspc; k++) {
                if((int)Le_s[k].size() != ngrd)
                    Le_s[k].resize(ngrd, 0.0);
                ofile << setw(19) << (Le_s.at(k)).at(i);
            }
            
            // Sc_s
            for(int k=0; k<nspc; k++) {
                if((int)Sc_s[k].size() != ngrd)
                    Sc_s[k].resize(ngrd, 0.0);
                ofile << setw(19) << (Sc_s.at(k)).at(i);
            }
            
        }
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

void odtline::readProperties(string fname) {

    *proc.ostrm << endl << "# Reading odtline property file " << fname << endl;
    
    ifstream ifile(fname.c_str()); 
    if(!ifile) {
        *proc.ostrm << "\n\n# ***************** ERROR OPENING FILE " << fname << endl << endl;
        exit(0);
    }
    
    string       s1;
    stringstream ss1;
    
    getline(ifile, s1);                        // read line "# grid point = 100"
    ss1.str(s1);
    ss1 >> s1 >> s1 >> s1 >> s1 >> ngrd;
    ss1.clear();
    
    getline(ifile, s1);                        // read line "# Domain Size = 2"
    ss1.str(s1);
    ss1 >> s1 >> s1 >> s1 >> s1 >> Ldomain;
    
    getline(ifile, s1);                        // read line "# pressure (Pa) = 101325 Pa"
    ss1.str(s1);
    ss1 >> s1 >> s1 >> s1 >> s1 >> pres;
    ss1.clear();
    getline(ifile, s1);
    
    
    ngrdf = ngrd + 1;
    
    pos.resize(ngrd);            // this should not be needed ?
    posf.resize(ngrdf);
    rho.resize(ngrd);
    molec.resize(ngrd);
    phase.resize(ngrd);
    lambda.resize(ngrd);
    if (odtP->Iparticles)
        voidFrac.resize(ngrd);
    if (LhasTemp)
        temp.resize(ngrd);
    if(LhasRxn) {
        temp.resize(ngrd);
        enth.resize(ngrd);
        mixf.resize(ngrd);
        chi.resize(ngrd);
        gradZ.resize(ngrd);
        diffZ.resize(ngrd);
        for(int k=0; k<nspc; k++)
            yspc[k].resize(ngrd);
        normalizeSpecies();
    }
    if(LhasVel) {
        uvel.resize(ngrd);    // this should not be needed ?
        vvel.resize(ngrd);
        wvel.resize(ngrd);
    }
    if(LhasEta) {
        for(int k=0; k<neta; k++)
            eta[k].resize(ngrd);
    }
    if(LhasMom) {
        for(int k=0; k<nmom; k++)
            mom[k].resize(ngrd);
    }
    
    if(LhasScl) {
        for(int k=0; k<nscl; k++)
            scl[k].resize(ngrd);
    }
    
    
    for(int i=0; i<ngrd; i++) {
        ifile >> pos[i]
              >> posf[i]
              >> rho[i]
              >> molec[i] 
              >> phase[i];
        if(LhasTemp)
            ifile >> lambda[i];
        if(LhasRxn)
            ifile >> mixf[i]
                  >> chi[i]
                  >> temp[i]
                  >> gradZ[i]
                  >> diffZ[i];
        if(LhasVel)
            ifile >> uvel[i]
                  >> vvel[i]
                  >> wvel[i];
        if(LhasTemp)
            ifile >> temp[i];
        if(LhasRxn){
            ifile >> enth[i];
            for(int k=0; k<nspc; k++)
                ifile >> yspc[k][i];
        }
        if(LhasEta){
            for(int k=0; k<neta; k++)
                ifile >> eta[k][i];
        }
        if(LhasMom){
            for(int k=0; k<nmom; k++)
                ifile >> mom[k][i];
        }
        
        if(LhasScl){
            for(int k=0; k<nscl; k++)
                ifile >> scl[k][i];
        }
        
    }
    posf[ngrd] = Ldomain + posf[0];
    
    ifile.close();

}

///////////////////////////////////////////////////////////////////////////////

/** Return the species vector at a given point.  This is needed
 *  because of the funny structure of the species organization on the grid.
 * 
 *  @param i \input grid point to return species vector at.
 *  @param yi \inout species vector at point i.
 */

void odtline::getYspVecAtPt(const int &i, vector<double> &yi) {
    
    if(!LhasRxn && !Lprxn) return;
    
    for(int k=0; k<nspc; k++)
        yi[k] = yspc[k][i];
    
}

///////////////////////////////////////////////////////////////////////////////

/** Set the temperature vector (auxiliary) from the odtline state. */

void odtline::setTempVec() {
    
    if(!LhasRxn) return;
    
    if((int)temp.size() != ngrd)
        temp.resize(ngrd);
    vector<double> yi(nspc);
    for(int i=0; i<ngrd; i++) {
        getYspVecAtPt(i, yi);
        gas->setState_PY(pres, &yi[0]);
        gas->setState_HP(enth[i], pres, 1.E-10);
        //gas->setState_HP(enth[i], pres);
        temp[i] = gas->temperature();
    }
    
}

///////////////////////////////////////////////////////////////////////////////

/** Reads an input file to initialize a scalar field. 
 *
 *  @param fname \input name of odtline setup input file.
 */

void odtline::setOdtline(string fname)  {
    
    caseInput_.setFile(fname);
    
    double delta_mixf;
    double fyc1;
    double fyc2;
    double delta_vel;          // velocity transition width
    double vel_diff;           // velocity difference
    double vel_min;            // minimum velocity
    double vyc1;               // first transition position
    double vyc2;               // second transition position
    
    double initUnburntZone;    //fraction of the domain with unburnt mixture (initial profile is a step function for premixed combustion)
    
    caseInput_.getParameter("probType",           &probType);
    caseInput_.getParameter("odtline_delta_mixf", &delta_mixf, 0.0);
    caseInput_.getParameter("odtline_fyc1",       &fyc1, 0.0);
    caseInput_.getParameter("odtline_fyc2",       &fyc2, 0.0);
    caseInput_.getParameter("odtline_delta_vel",  &delta_vel, 0.0);
    caseInput_.getParameter("odtline_vel_diff",   &vel_diff, 0.0);
    caseInput_.getParameter("odtline_vel_min",    &vel_min, 0.0);
    caseInput_.getParameter("odtline_vyc1",       &vyc1, 0.0);
    caseInput_.getParameter("odtline_vyc2",       &vyc2, 0.0);
    
    //  caseInput_.setFile(caseSetupFile);
    caseInput_.getParameter("odtline_unifCompress",	&unifCompress, 0);			// input for uniform compression or velocity based: 1 or 0
    caseInput_.getParameter("odtline_mixf_premixed",	&mixf_premixed, 0.);		// for pre-mixed flow mixture fraction: 0<mixf_premixed<1
    caseInput_.getParameter("odtline_nonPremixed",	&nonPremixed, 1);			// for pre-mixed flow mixture fraction: 1 or 0
    
    //---------- for premixed combustion
    caseInput_.getParameter("initUnburntZone", &initUnburntZone, 1.0);
    //---------- set the mixture fraction and velocity profiles
    
    mixf.resize(ngrd);            // not normally used for !LhasRxn, but size here for generic initialization
    if(LhasVel) uvel.resize(ngrd);
    
    //---------- a mixing layer with a tanh transition
    
    if(probType == 1 || probType==5) {      
    
        fyc1 += 0.5*(posf[0]+posf[ngrd]); // fyc1=dist from center --> dist from left
        for(int i=0; i<ngrd; i++)
            mixf[i] = 0.5*(1.0+tanh(2.0/delta_mixf*(pos[i]-fyc1)));
        
        vyc1 += 0.5*(posf[0]+posf[ngrd]);  // vyc1=dist from center --> dist from left
        for(int i=0; i<ngrd; i++){
            uvel[i] = vel_diff*0.5*(1.0+tanh(2.0/delta_vel*(pos[i]-vyc1))) + vel_min;
            if (probType==5) 
                uvel[i] *= tanh(2.0/delta_vel*pos[i]) * tanh(2.0/delta_vel*(Ldomain-pos[i]));
        }
    }
    
    //---------- a jet with a tanh transition
    
    else if(probType == 2) {
    
        fyc1 += 0.5*(posf[0]+posf[ngrd]); // fyc1=dist from center --> dist from left
        fyc2 += 0.5*(posf[0]+posf[ngrd]); // fyc2=dist from center --> dist from left
        for(int i=0; i<ngrd; i++){
            mixf[i] = 0.5*(1.0+tanh(2.0/delta_mixf*(pos[i]-fyc1))) * 
                        0.5*(1.0+tanh(2.0/delta_mixf*(fyc2-pos[i])));
        }
        
        vyc1 += 0.5*(posf[0]+posf[ngrd]);  // vyc1=dist from center --> dist from left
        vyc2 += 0.5*(posf[0]+posf[ngrd]);  // vyc2=dist from center --> dist from left
        if (LhasVel) {
            
            for(int i=0; i<ngrd; i++){
                uvel[i] = vel_diff * 0.5 * (1.0 + tanh(2.0 / delta_vel * (pos[i] - vyc1))) *
                          0.5 * (1.0 + tanh(2.0 / delta_vel * (vyc2 - pos[i]))) + vel_min;
            }
        }
    }
    
    //---------- flame wall case tanh mixture fraction
    //----------    laminar velocity, boundary layer solution: 
    //----------    Fit Blasius solution: u/U = fifth order polynomial in eta, where
    //----------    eta = y*sqrt(U/nu*x) and x is the downstream dist, taken from vyc1; 
    //----------    U is vel_diff
    
    else if(probType == 3) {             // flame wall case tanh mixture fraction
    
        // fyc1 is distance from left
        for(int i=0; i<ngrd; i++)
            mixf[i] = 1.0-0.5*(1.0+tanh(2.0/delta_mixf*(pos[i]-fyc1))); //eim check for isothermal wall
        
        double nu = 4.765E-5;           // kinematic viscosity m2/s air at 1 atm, 300 K
        double EtaFac = sqrt(vel_diff/nu/vyc1);
        double Eta;
        
        for(int i=0; i<ngrd; i++) {
            Eta = EtaFac * pos[i];
            if(Eta < 7.4)
                uvel[i] = vel_diff*( -1.55926E-4*pow(Eta,5.0) +
                               3.70924E-3*pow(Eta,4.0) -
                               2.84820E-2*pow(Eta,3.0) +
                               4.77151E-2*pow(Eta,2.0) +
                               3.05872E-1*Eta +
                               2.02326E-3 );
            else
                uvel[i] = vel_diff;
        }
    }
    
    else if(probType == 4) {    // opposed jets set up for spatial initial condition
        if(LhasRxn) {           // set spatial mixture fraction profile  for non-premixed case
            for(int i=0; i<ngrd; i++) {
                if(nonPremixed)
                {
                    mixf[i] = 1-( pos[i]-posf[0] ) / ( posf[ngrd]-posf[0] );	// set linear profile: fuel ----- air
                    // double distance=posf.at(ngrd)-posf.at(0);
                    // double temp= 20.0 / distance *pos[i]-10.0;
                    // mixf[i] = 0.5-0.5*tanh( temp ); 
                    // mixf[i]=0;
                }
                else
                    mixf[i] = mixf_premixed;    // set constant profile
                
                temp[i] = strm->T1 + ( (pos[i]-posf[0]) / (posf[ngrd]-posf[0]) ) * ( strm->T0 - strm->T1 );	// set linear temperature profile
            }
        }
        if(LhasVel){
            for(int i=0; i<ngrd; i++)
            {
                uvel[i] = odtP->uBClo + ( (pos[i]-posf[0]) / (posf[ngrd]-posf[0]) ) * ( odtP->uBChi - odtP->uBClo );	// set linear velocity profile
                //uvel[i] = 0.01 + ( (pos[i]-posf[0]) / (posf[ngrd]-posf[0]) ) * ( -0.01 - 0.01 );			// set linear velocity profile
                //if( pos[i] >=0.25 && pos[i]<= 0.35 )
                //uvel[i]=0.1;
            }
        }
    }
    
    else if(probType == 6) { // freely propagating premixed flame; Only working for CH4 at the moment
        
        if(!LhasRxn){
            *proc.ostrm << "\n\nError setting the odtline: LhasRxn is false, but premixed must have reaction" << endl;
            exit(0);
        }
#ifdef DOCANTERA
        vector<double> Y1(nspc,0.0);
        vector<double> Y2(nspc,0.0);
        for(int i=0; i<ngrd; i++) {
        
            temp[i] = odtP->inletTemp; 
            if(odtP->nInletFuel != 0 && odtP->nInletOxidizer != 0){
                gas->setMoleFractionsByName("C2H4:0.33, O2:1, N2:3.76");
                gas->getMassFractions(&Y1[0]);
                Y1[gas->speciesIndex("C2H4")] *= odtP->inletEquRatio;
            }
            else {
                *proc.ostrm << "\n\n Not available yet";
                exit(0);
            }
            
            gas->setMassFractions(&Y1[0]);
            gas->setTemperature(temp[i]);
            gas->setPressure(pres);
            //equilibrate(*g,"TP");
            gas->getMassFractions(&Y1[0]);
            rho[i]  = gas->density();
            enth[i] = gas->enthalpy_mass(); 
            for(int k=0;k<nspc;k++){
                yspc[k][i] = Y1[k];
            }
        }
        int nInitUnburntZone = int(initUnburntZone*ngrd);
        for(int i=nInitUnburntZone; i<ngrd; i++) {
            gas->setMassFractions(&Y1[0]);
            gas->setTemperature(temp[i]);
            gas->setPressure(pres);
            equilibrate(*gas,"HP");
            gas->getMassFractions(&Y2[0]);
            temp[i] = gas->temperature(); 
            enth[i] = gas->enthalpy_mass();  
            rho[i] = gas->density();
            for(int k=0;k<nspc;k++){
                yspc[k][i] = Y2[k];
            }
        }
#endif   // DOCANTERA
        
        //
    }

    else if (probType == 7) { // homogeneous turbulence with initial sine wave u velocity
        for(int i=0; i<ngrd; i++) {
            double r = static_cast<double>(proc.myid)/static_cast<double>(proc.nproc);
            uvel[i] = 6.55*sin(2.0*M_PI*(pos[i]+r*0.0254)/0.0254);
        }
    }

    //---------- Error condition
    
    else {
        *proc.ostrm << "\n\nError setting the odtline: option " << probType 
                    << " not recognized " << endl;
        exit(0);
    }
    
    //----------- set the gas state along the line
    
    if(LhasRxn && probType!=6) {
        pres = strm->pres;                      // this is redundant
        
        vector<double> yspecies(nspc);          // dummy storage, 		corrected from ngrd-->nspc	ZJ
        
        for(int i=0; i<ngrd; i++) {
            strm->getProdOfCompleteComb( mixf[i], yspecies, enth[i], temp[i], probType ); 
            for(int k=0; k<nspc; k++)
                yspc[k][i] = yspecies[k];
        }
    }
    
    //----------- set anyline objects
    
    
    if(odtP->IetaType == odtP->IETATYPE_TABLELOOKUP ){
        
        pres = strm->pres;                      // this is redundant
        eta[0] = mixf;
        
        enth.resize(ngrd);                      
        for(int i=0; i<ngrd; i++)
            enth[i] = (strm->h0 * (1 - mixf[i]) + strm->h1 * mixf[i]);
        eta[1] = enth;
    }
    if(!odtP->ItableLookup){
        setRhoVec();
        setViscosity();     // molec
    }
    
}

///////////////////////////////////////////////////////////////////////////////

/**Set density based on enthalpy, pressure, Yi. */

void odtline::setRhoVec() {
    
    if(!LhasRxn) {
        if(odtP->ItableLookup){
            etaTools->updateOdtLineVecs(true);
            return;
        }
        else
            rho = vector<double>(ngrd,odtP->rho_0);
        return;
    }
    
    if((int)rho.size() != ngrd)
        rho.resize(ngrd);
    
    vector<double> yi(nspc);
    for(int i=0; i<ngrd; i++) {
        getYspVecAtPt(i, yi);
        gas->setState_PY(pres, &yi[0]);
        gas->setState_HP(enth[i], pres, 1.E-10);
        //gas->setState_HP(enth[i], pres);
        rho[i] = gas->density();
    }
    
}

///////////////////////////////////////////////////////////////////////////////

/** Set member molec from the current state on the points of the line. */

void odtline::setViscosity() {
    
    if(!LhasRxn) {
        molec = vector<double>(ngrd,odtP->visc_0);
        return;
    }
    
    vector<double> yi(nspc);           // working array
    
    for(int i=0; i<ngrd; i++) {
        getYspVecAtPt(i, yi);
        gas->setState_PY(pres, &yi[0]);
        gas->setState_HP(enth[i], pres, 1.E-10);
        //gas->setState_HP(enth[i], pres);
        molec[i]   = tran->viscosity();
    }
}

///////////////////////////////////////////////////////////////////////////////

/** Set member mixf (mixture fraction) from the composition on the line. */ 

void odtline::setMixfVec() {
    
    if(!LhasRxn) return;
    
    mixf.resize(ngrd);
    
    vector<double> yi(nspc);
    for(int i=0; i<ngrd; i++) {
        getYspVecAtPt(i, yi);
        mixf[i] = strm->getMixtureFraction(&yi[0]);
    }
    
}

///////////////////////////////////////////////////////////////////////////////

/** Set scalar dissipation rate (chi) \cond
 *  Chi = 2D(df/dx)^2                    
 *  Compute the derivative as follows:
 *        *       *                 *
 *       i-1      i                i+1
 *            a            b 
 *  Linearly interpolate the derivatives at a and b to point i
 *  f' = ( d2*fa' + d1*fb' ) / (d1+d2) where d1 is dist between i and a
 *  and d2 is distance between b and i \endcond
 */

void odtline::setChiVec() {

    if(!LhasRxn) return;
    
    chi.resize(ngrd);
    gradZ.resize(ngrd);
    
    vector<double> D(ngrd);   
    
    //------------- Get diffusivity = thermal diffusivity
    
    for(int i=0; i<ngrd; i++) {
        D[i] = getThermalDiff(i);
    }
    
    //------------- Compute chi

    gradZ[0] = (mixf[1]-mixf[0])/(pos[1]-pos[0]); 
    double d1, d2;
    for(int i=1; i<ngrd-1; i++) {
        d1 = 0.5*(pos[i]-pos[i-1]);
        d2 = 0.5*(pos[i+1]-pos[i]);
        gradZ[i] = (d2*d2*(mixf[i]-mixf[i-1]) + d1*d1*(mixf[i+1]-mixf[i]))/(2*d1*d2*(d1+d2));
    }
    gradZ[ngrd-1] = (mixf[ngrd-1]-mixf[ngrd-2])/(pos[ngrd-1]-pos[ngrd-2]);

    for(int i=0; i<ngrd; i++)
        chi[i] = 2.0*D[i]*gradZ[i]*gradZ[i];
    
//    vector<double> dfdx(ngrd);
//    dfdx[0] = (mixf[1]-mixf[0])/(pos[1]-pos[0]); 
//    double d1, d2;
//    for(int i=1; i<ngrd-1; i++) {
//        d1 = 0.5*(pos[i]-pos[i-1]);
//        d2 = 0.5*(pos[i+1]-pos[i]);
//        dfdx[i] = (d2*d2*(mixf[i]-mixf[i-1]) + d1*d1*(mixf[i+1]-mixf[i]))/(2*d1*d2*(d1+d2));
//    }
//    dfdx[ngrd-1] = (mixf[ngrd-1]-mixf[ngrd-2])/(pos[ngrd-1]-pos[ngrd-2]);
//    
//    for(int i=0; i<ngrd; i++)
//        chi[i] = 2.0*D[i]*dfdx[i]*dfdx[i];
    
}

///////////////////////////////////////////////////////////////////////////////

/** Set mixture fraction diffusion  (diffZ) \cond
 *  DiffZ = (d/dx)(D(df/dx))                    

 *  Compute the derivative as follows:
 *        *       *                 *
 *       i-1      i                i+1
 *            a            b 
 *  Linearly interpolate the derivatives at a and b to point i
 *  f' = ( d2*fa' + d1*fb' ) / (d1+d2) where d1 is dist between i and a
 *  and d2 is distance between b and i \endcond
 */

void odtline::setDiffZVec() {

    if(!LhasRxn) return;
 
    diffZ.resize(ngrd);

    vector<double> D(ngrd);   

    //------------- Get diffusivity = thermal diffusivity

    for(int i=0; i<ngrd; i++) {
        D[i] = getThermalDiff(i);
    }

    //------------- Compute diffZ

    vector<double> dfdx(ngrd);
    vector<double> dx(ngrd);
    vector<double> Dcoeff(ngrd);


    dfdx[0] = (mixf[1]-mixf[0])/(pos[1]-pos[0]); 

    for(int i=0; i<ngrd-1; i++) {
        dx[i] = ( pos[i+1] - pos[i] );
        Dcoeff[i] = 0.5 * ( D[i+1] + D[i] );
        dfdx[i] = (mixf[i+1]-mixf[i])/dx[i];
    }
    dfdx[ngrd-1] = (mixf[ngrd-1]-mixf[ngrd-2])/(pos[ngrd-1]-pos[ngrd-2]);

    diffZ[0] = ( Dcoeff[1]*dfdx[1] - Dcoeff[0]*dfdx[0] ) / ( pos[1] - pos[0] );
    for(int i=1; i<ngrd-1; i++)
      diffZ[i] = ( Dcoeff[i+1]*dfdx[i+1] - Dcoeff[i]*dfdx[i] ) * 2. / ( pos[i+2] - pos[i] );
    diffZ[ngrd-1] = ( Dcoeff[ngrd-1]*dfdx[ngrd-1] - Dcoeff[ngrd-2]*dfdx[ngrd-2] ) / ( pos[ngrd-1] - pos[ngrd-2] );
}

///////////////////////////////////////////////////////////////////////////////

/** Return the mean molecular weight array from the composition on the line. 
 *
 *  @return vector of mean molecular weights on the line.
 */

vector<double> odtline::getMMW() {
    
    if(!LhasRxn) return vector<double>();
    
    vector<double> mmw(ngrd);
    vector<double> yi(nspc);
    for(int i=0; i<ngrd; i++) {
        getYspVecAtPt(i, yi);
        gas->setState_PY(pres, &yi[0]);
        mmw[i] = gas->meanMolecularWeight();
    }
    return mmw;
    
}

///////////////////////////////////////////////////////////////////////////////

/** Just checks for negative species. */
    
void odtline::enforceYsBounds(vector<vector<double> > &Y_old)  {

    if(!LhasRxn) return;

    bool eflag;
    if(odtP->LlimitMassFrac!=0){
        for(int i=0; i<ngrd; i++) {
            for(int k=0; k < nspc; k++) {
                eflag = false;
                if(yspc[k][i] < 0.0 && yspc[k][i] > -1.0e-15) {
                    yspc[k][i]=0.0;
                    eflag = true;
                }
                if(yspc[k][i] < 0.0 || yspc[k][i] > 1.0) { 
                    if(odtP->LlimitMassFrac==2)
                        *proc.ostrm << endl << "Species < 0 or > 1 : " << propNames[k+1] << " = " << yspc[k][i];
                    for(int k=0; k < nspc; k++) yspc[k][i]=Y_old[k][i];
                    break;
                }
            }
            
            if(eflag) 
                normalizeSpecies(i);
        }
        return;
    }
    
    for(int i=0; i<ngrd; i++) {
    
        eflag = false;
        for(int k=0; k<nspc; k++) 
            if(yspc[k][i] < 0.0) {
                yspc[k][i] = 0.0;
                eflag = true;
            }
        if(eflag) 
            normalizeSpecies(i);
        
    }
    
    
    
}

///////////////////////////////////////////////////////////////////////////////

/** Normalize species at grid point i.
 * without a parameter, use the default -1 which gives all grid points
 *
 * @param ipos \input grid point to normalize species on.
 */

void odtline::normalizeSpecies(int ipos) {
    
    if(!LhasRxn) return;
    
    int ilo = (ipos<0) ? 0      : ipos;
    int ihi = (ipos<0) ? ngrd-1 : ipos;
    
    double sum;
    for(int i=ilo; i<=ihi; i++) {
        sum = 0.0;
        for(int k=0; k<nspc; k++) 
            sum += yspc[k][i];
        if(sum > 0.0)
            for(int k=0; k<nspc; k++)
                yspc[k][i] /=sum;
    }
    
}

///////////////////////////////////////////////////////////////////////////////

/** Sets the variable to do mesh adaption on (sum of two or more normalized vars)
 *
 */

void odtline::setAdptVar() {                       
    
    adptVar.resize(ngrd);
    
    if(LhasRxn)
        adptVar = rho;
    if(LhasVel)
        adptVar = uvel;
    
}

///////////////////////////////////////////////////////////////////////////////
/** Sets the mesher array of pointers to vectors to adapt on.
 *  Can have any number for a given line. 
 *  If you use something like temperature, then you might need to make sure temperature
 *     is set before the adaption happens.
 *  If the oldmesher is being used, then adaption only happens on the first 
 *     variable in the list.
 *
 *  @param Ladpt \input false if no mesh adaption wanted (e.g. for eddies)
*/

void odtline::setMesher(bool Ladpt) {
    
    meshAdapter.phi.resize(0);
    
    if(LhasVel && (LhasRxn || odtP->ItableLookup))  {
        meshAdapter.phi.push_back(&temp);
        //meshAdapter.phi.push_back(&rho);
        meshAdapter.phi.push_back(&uvel);
        
        //  for(int i=0; i< yspc.size(); i++)
        //  meshAdapter.phi.push_back(&yspc[i]);		//for opposed flow also check mass fractions
        
    }
    else if(LhasVel && LhasTemp) {
        meshAdapter.phi.push_back(&uvel);
        meshAdapter.phi.push_back(&temp);
    }
    else if(LhasVel) {
        meshAdapter.phi.push_back(&uvel);
    }
    else if(LhasRxn) {
        meshAdapter.phi.push_back(&rho);
    }
    else if(odtP->Llem)
        meshAdapter.phi.push_back(&scl[0]);
    else {
       *proc.ostrm << "\n\n***************** ERROR no LhasRxn or LhasVel\n";
       exit(0);
    }
    
    if(LhasEta) {
        meshAdapter.phi.push_back(&eta[0]);
    }
    
#ifndef COMPSCI
    meshAdapter = adaptMesh(this, odtP, meshAdapter.phi, Ladpt);
#else
    // performance notes:
    // * meshAdapter is already a valid instance, no need to create a temporary one
    // and copy over the contents. (compiler might fix this though)
    // * every time a copy of mesherPhi must be created for the call to the adaptMesh constuctor,
    // but the contents of mesherPhi do not actually change between iterations
    meshAdapter.init(this, odtP, Ladpt);
#endif
}

///////////////////////////////////////////////////////////////////////////////

/** Sets the cell void fraction for particle cases.                                                        \cond
 *  Void frac = (Vcell - Vpart_in_cell) / (Vcell) = 1-sum(4/3 * pi * r_i^3*n_i)/dx                         \endcond
 *     Sum is over each particle in a given grid cell.
 *                                                                                                         
 *  @param part \inout particles object to use
 */

void odtline::setVoidFrac(particles *part) {
    
    if (!odtP->Iparticles) return;
    part->set_iyPos();
    
    voidFrac.resize(ngrd);
    for(int i = 0; i<ngrd; i++){
        voidFrac[i] = 1;
    }
    return;
}

///////////////////////////////////////////////////////////////////////////////

/**
 * Computes the thermal diffusivity at index i
 * @param i \input Index in the domain at which thermal diffusivity is computed.
 */

double odtline::getThermalDiff(int& i) {
    vector<double> yi(nspc);
    double lambda = 1;
    double cp = 0;
    getYspVecAtPt(i, yi);
    gas->setState_PY(pres, &yi[0]);
    gas->setState_HP(enth[i], pres, 1.E-10);
    //gas->setState_HP(enth[i], pres);
    temp[i] = gas->temperature();
    
    lambda = tran->thermalConductivity();
    cp = gas->cp_mass();
    return lambda / rho[i] / cp;
}

////////////////////////////////////////////////////////////////////////////

/**
 * Method to get the total gas enthalpy.
 * @return the total gas enthalpy (J)
 */
double odtline::getTotalGasEtaEnthalpy(){
    double totalGasEnth = 0;
    vector<double> myVoidFrac(ngrd,1.0);
    if(odtP->Iparticles){
        myVoidFrac = voidFrac;
    }
    for(int i=0; i<ngrd;i++){
       totalGasEnth += eta[1][i] * rho[i] * (posf[i+1] - posf[i]) * myVoidFrac[i]; 
    }
    return totalGasEnth;
}

////////////////////////////////////////////////////////////////////////////

/**
 * Method to get the total mass of the gas.
 * @return the total mass of the gas (kg)
 */
double odtline::getTotalGasMass(){
    double totalGasMass = 0.0;
    vector<double> myVoidFrac(ngrd, 1.0);
    if (odtP->Iparticles) {
        myVoidFrac = voidFrac;
    }
    for(int i=0;i<ngrd;i++){
        totalGasMass += rho[i] * (posf[i+1] - posf[i]) * myVoidFrac[i];
    }
    return totalGasMass;
}

///////////////////////////////////////////////////////////////////////////////
/** This function increases the domain size of the odtline with one or two 
 *  cells depending on the input LUC. There will be add one cell at the lower
 *  boundary and/or one cell at the upper boundary. The values of the new
 *  cell/cells is/are not set.
 *  @Param  delz    Size which is added to the current domain size
 *  @Param  LUC     Switch for expandung at lower boundary 'L', at upper 
 *                  boundary 'U', or the old line centered 'C'.
 */
void odtline::expand(const double delz_in, const char LUC){
    
    double delz = delz_in;
    if (LUC == 'C')
        delz = delz_in/2.0;
    *proc.ostrm << endl << "###############################";
    *proc.ostrm << endl << "# delz_in = " << delz_in;
    *proc.ostrm << "  # delz = " << delz;
    *proc.ostrm << "  # char = " << LUC;
    
    if (LUC == 'L' || LUC == 'C'){
        *proc.ostrm << endl << "# expand at lower bc";
        // moving positions of current cells
        for(int i=0; i<ngrd; i++){
            pos.at(i) += delz;
            posf.at(i) += delz;
        }
        posf.at(ngrd) += delz;
        
        // inserting one cell at the beginning
        pos.insert(pos.begin(), 0.5*delz);
        posf.insert(posf.begin(), 0.0);
        ngrd++;
        ngrdf++;
        
        // changing domain size
        odtP->domainLength += delz;
        Ldomain += delz;
        
        // inserting cell in property arrays with default value
        uvel.insert(uvel.begin(), 0.0);
        vvel.insert(vvel.begin(), 0.0);
        wvel.insert(wvel.begin(), 0.0);
        rho.insert(rho.begin(), 0.0);
        molec.insert(molec.begin(), 0.0);
        phase.insert(phase.begin(), 0.0);
        lambda.insert(lambda.begin(), 0.0);
        if(LhasTemp)
            temp.insert(temp.begin(), 0.0);
    }
    if (LUC == 'U' || LUC == 'C'){
        *proc.ostrm << endl << "# expand at upper bc";
        // inserting one cell at the end
        pos.push_back(Ldomain + 0.5*delz);
        posf.push_back(Ldomain + delz);
        ngrd++;
        ngrdf++;
        
        // changing domain size
        odtP->domainLength += delz;
        Ldomain += delz;
        
        // inserting cell in property arrays with default value
        uvel.push_back(0.0);
        vvel.push_back(0.0);
        wvel.push_back(0.0);
        rho.push_back(0.0);
        molec.push_back(0.0);
        phase.push_back(0.0);
        lambda.push_back(0.0);
        if(LhasTemp)
            temp.push_back(0.0);
    }
    *proc.ostrm << endl << "# expanding domain finished";
    *proc.ostrm << endl << "###############################";
}


///////////////////////////////DOXYGEN DOCUMENTATION/////////////////////////

/*! \fn void odtline::setChiVec()
 * \f[
 *  \chi = 2D ( \frac{df}{dx} )^2
 * \f]
 *  Compute the derivative as follows:
 * \vc{
 *        *       *                 *
 *       i-1      i                i+1
 *            a            b
 * }
 *  Linearly interpolate the derivatives at \c a and \c b to point <code>i</code>.
 *  \f[
 *      f' = \frac{ d_2 f'_a + d_1 f'_b }{ d_1+d_2 }
 *  \f]
 *  where \fun{d_1} is distance between \c i and \c a
 *  and   \fun{d_2} is distance between \c b and \c i
 */
 
/*! \fn void odtline::setVoidFrac(particles *part) 
 *  \f[
 *     \text{Void frac} = \frac{V_{cell} - V_{part in cell} }{ V_{cell}} = 1 - \frac{ \displaystyle\sum\limits_i{ \frac{4}{3}  \pi  {r_i}^3 n_i}}{dx} 
 *  \f]
 */
///////////////////////////END DOXYGEN DOCUMENTATION//////////////////////////
