/**
 * @file MOM_soot.cc
 * Source file for class MOM_soot
 */

#include "MOM_soot.h"
#include "processor.h"
#include "diffuser.h"
#include "table.h"
#include <cmath>
#include <cstdlib>

using namespace std;

extern processor proc;

///////////////////////////////////////////////////////////////////////////////
/** Constructor function.
 *
*/

MOM_soot::MOM_soot(odtline *odtlp, std::string ftable){
    
    odtl         = odtlp;
    // *diff is set externally
    
    mom_source   = std::vector<std::vector<double> > (odtl->nmom, std::vector<double> (odtl->ngrd, 0.0)); 
    //gasSp_source = std::vector<std::vector<double> > (odtl->nspc, std::vector<double> (odtl->ngrd, 0.0)); 

    
    //----------- set default moments

    rho_s       = 1850.0;  //solid soot density, kg/m^3 
    kb          = 1.38066E-23;         // J/K
    Na          = 6.02257E26;          // particles/kmol
    Ca          = 9.0;
    Cmin        = 100.0;
    if(!odtl->odtP->ItableLookup) 
    {
        iO2         = odtl->gas->speciesIndex("O2");
        iO2         = (iO2 >= 0) ? iO2 : odtl->gas->speciesIndex("o2");
    } else {
        iO2         = odtl->etaTools->luTable->index.species_mass_index;
    }
    MWc     = odtl->gas->atomicWeight(odtl->gas->elementIndex("C"));

    //--------- init table properties
    
    readYc2h2Table(ftable);

}

///////////////////////////////////////////////////////////////////////////////


/** Update fluxes at cell faces
 *
 *  @param dd    - input  - the grad part of grad phi (multiplier
 */

void MOM_soot::computeFluxes(std::vector<double> &dd){

    int i, im, k;

    resizeVar();

    std::vector<double> viscf(odtl->ngrd);
    std::vector<double> rhof(odtl->ngrd);
    for(i=0; i<odtl->ngrd; i++){
        viscf[i] = diff->linearInterpToFace(i, *odtl, odtl->molec);
        rhof[i]  = diff->linearInterpToFace(i, *odtl, odtl->rho);
    }

    //////////////////// interior faces ////////////////////

    for(i=1, im=0; i<odtl->ngrd; i++, im++) {

        //-------------- Thermophoretic transport is "hyperbolic", flux looks like f=v*phi, so upwind it for stability

        double vel = -0.55415*viscf[i]*(log(odtl->temp[i])-log(odtl->temp[im]))*dd[i];
        double ii = (vel > 0) ? im : i;
        for(k=0; k<odtl->odtP->Imom; k++) {
            diff->flxProp[diff->iptMom+k][i] = vel * odtl->mom[k][ii];
        }
    } 

    //////////////////// boundary faces ////////////////////

    //-------------- periodic 

    if (diff->odtP->bcType == odtParam::BCTYPE_PERIODIC) { 

        i = 0;
        im = odtl->ngrd - 1;

        double vel = -0.554*viscf[i]*(log(odtl->temp[i])-log(odtl->temp[im]))*dd[0];
        double ii = (vel > 0) ? im : i;

        for(k=0; k<odtl->odtP->Imom; k++) {
            diff->flxProp[diff->iptMom+k][i]            = vel * odtl->mom[k][ii];
            diff->flxProp[diff->iptMom+k][odtl->ngrd]   = diff->flxProp[diff->iptMom+k][i];
        }
    }

    //----------------- wall on left outflow on right or outflow, etc.

    else if (diff->odtP->bcType == odtParam::BCTYPE_WALL    ||
             diff->odtP->bcType == odtParam::BCTYPE_OUTFLOW ||
             diff->odtP->bcType == odtParam::BCTYPE_INLET   ||
             diff->odtP->bcType == odtParam::BCTYPE_WALLOUTFLOW) { 

        for (int k = 0; k < odtl->odtP->Imom; k++) {
            diff->flxProp[k][0] = 0.0;
            diff->flxProp[k][odtl->ngrd] = 0.0;
        }

    }

    //------------------ unknown bc

    else { 

        *proc.ostrm << endl << "****** ERROR IN MOM_soot, BCTYPE = " << diff->odtP->bcType << " unknown" << endl;
        exit(0);

    }


} // end function

///////////////////////////////////////////////////////////////////////////////

/** Leung and Lindstedt Soot Model
    Calculates the nucleation, growth, oxidation and coagulation rates.
    This calculates for the entire line.
    Then from those the source terms for number density and mass fraction 
    calculated
 */

void MOM_soot::computeSourceTerms() {
    
    resizeVar();
    odtl->setMixfVec();                 // needed for scaling flmlt c2h2

    for(int i=0; i<odtl->ngrd; i++){

        if(odtl->odtP->Imom == 1) {
            getSootSourceTermsPt1mom(diff->rhsSrc[diff->iptMom][i], i);
        }
        else if(odtl->odtP->Imom == 2) {
            getSootSourceTermsPt2mom(diff->rhsSrc[diff->iptMom][i], diff->rhsSrc[diff->iptMom+1][i], i);
        }
        else if(odtl->odtP->Imom == 3) {
            getSootSourceTermsPt3mom(diff->rhsSrc[diff->iptMom][i], diff->rhsSrc[diff->iptMom+1][i], 
                                     diff->rhsSrc[diff->iptMom+2][i], i);
        }
        else {
            cout << "\n\n***************** ERROR: wrong Imom in MOM_soot" << endl;
            exit(0);
        }


                                    //momsource and gSp_source are outputs
        //for(int k = 0; k<odtl->nspc; k++)
        //    gasSp_source[k][i] = 0.0;
        //gasSp_source[iO2][i] = gSp_source[iO2]; 

    }
            //exit(0);//doldb
    //--------Add soot source term from the particles
    if(odtl->odtP->Iparticles && odtl->odtP->Lprxn) {
        for(int i = 0; i < odtl->ngrd; i++) 
            for(int k = 0; k < odtl->nmom; k++) 
                diff->rhsSrc[diff->iptMom + k][i] += (*(diff->gSource))[diff->part->iPtMom + k][i]; //add the source term for eta
    }

    return;
}

///////////////////////////////////////////////////////////////////////////////

/** Single moment source terms at a point
 *
 *  Used for the fire cases.  
 *  Used for table lookup, (or detailed chemistry)
 *  Growth is based on particle reaction rate using a constant soot yield,
 *       assumed small enough to ignore mass loss from the gas.
 *  Oxidation is using O2, OH, O rates as in Lignell 2007 2-D soot paper.
 * 
 */

void MOM_soot::getSootSourceTermsPt1mom(double &Ys_source, const int ipos) {

    if(!odtl->odtP->Iparticles) { 
        *proc.ostrm << endl << "****** ERROR IN MOM_soot::getSootSourceTermsPt1mom, must use particles here" << endl;
        exit(0);
    }

    //---------- Soot growth 

    double mdotPartRxn = (*(diff->gSource))[diff->part->iPtMass][ipos]; // kg/m3*s
    double Ys_yield = 0.02;

    Ys_source = mdotPartRxn * Ys_yield;                            // kg/m3*s, divide by rho below

    //---------- Soot oxidation

    double xO2, xOH, xO, cO2;

    { //------------ get species concentrations for oxidation

        //double mixf = odtl->odtP->ItableLookup ? odtl->eta[0][ipos] : odtl->mixf[ipos]; //  !!!!!     currently unused variable

        int    iO2, iOH, iO;
        iO2 = odtl->gas->speciesIndex("O2");
        iO2 = (iO2 >= 0) ? iO2 : odtl->gas->speciesIndex("o2");
        iOH = odtl->gas->speciesIndex("OH");
        iOH = (iOH >= 0) ? iOH : odtl->gas->speciesIndex("oh");
        iO  = odtl->gas->speciesIndex("O");
        iO  = (iO  >= 0) ? iO  : odtl->gas->speciesIndex("o");

        double mmw  = odtl->gas->meanMolecularWeight();

        if(!odtl->odtP->ItableLookup) { 
            double yO2, yOH, yO;
            yO2 = odtl->yspc[iO2][ipos];
            yOH = odtl->yspc[iOH][ipos];
            yO  = odtl->yspc[iO ][ipos];
            double mwO2 = odtl->gas->molecularWeight(iO2);
            double mwOH = odtl->gas->molecularWeight(iOH);
            double mwO  = odtl->gas->molecularWeight(iO);
            xO2 = yO2 * mmw/mwO2;
            xOH = yOH * mmw/mwOH;
            xO  = yO  * mmw/mwO;
        }
        else {
            xO2 = odtl->etaTools->luTable->
                  getValAtGridPoint(ipos,odtl->etaTools->luTable->index.species_mole_index+iO2) ;
            xOH = odtl->etaTools->luTable->
                  getValAtGridPoint(ipos,odtl->etaTools->luTable->index.species_mole_index+iOH) ;
            xO  = odtl->etaTools->luTable->
                  getValAtGridPoint(ipos,odtl->etaTools->luTable->index.species_mole_index+iO) ;
        }

        cO2 = odtl->rho[ipos]/mmw * xO2;

    } //------------ end get species concentrations

    double dp_soot = 40.E-9;                                       // m
    double Ys = odtl->mom[0][ipos];       
    double Asoot = 6.0*odtl->rho[ipos]*Ys/rho_s/dp_soot;           // m2 soot / volume
    double T = odtl->temp[ipos];
    double rO2, rOH, rO;


    rO2 = 10000.0 * exp(-19640/T) * sqrt(T) * Asoot * cO2;         // kmol/m3*s
    rOH = 106.0 * 0.13 / sqrt(T) * Asoot * xOH;                    // kmol/m3*s
    rO  = 55.4  * 0.5  / sqrt(T) * Asoot * xO;                     // kmol/m3*s

    Ys_source -= (rO2 + rOH + rO) * MWc;                           // kg/m3*s

    //----------- finish up: 

    Ys_source /= odtl->rho[ipos];                                  // 1/s

}


///////////////////////////////////////////////////////////////////////////////

/** Leung and Lindstedt Soot Model
 *  Calculates the nucleation, growth, oxidation and coagulation rates.
 *  This calculates for the entire line.
 *  Then from those the source terms for number density and mass fraction 
 *  calculated
 *
 *  d[n]/dt      = S_n   (=) part/m3*s
 *  d[rho*Ys]/dt = S_Ys  (=) kg/m3*s
 *  ODT form:
 *  d[rho*(n/rho)]/dt = S_n,  where (n/rho) = Mhat_0 = M_0/rho, where M_0 = n
 *  d[rho*(Ys)]/dt    = S_Ys, where (Ys)    = Mhat_1 = M_1/rho, where M_1 = rho*Ys
 *  So, the natural source terms are the same
 */

void MOM_soot::getSootSourceTermsPt2mom(double &mhat0_source,
                                        double &mhat1_source,
                                        //vector<double> &gSp_source,
                                        const int ipos) {

    //set working variables
    double nd = odtl->mom[0][ipos] * odtl->rho[ipos];
    
    double Ys = odtl->mom[1][ipos];
    double mixf = odtl->odtP->ItableLookup ? odtl->eta[0][ipos] : odtl->mixf[ipos];

    //get acetylene concentration

    double yC2H2;
    int iC2H2;
    iC2H2   = odtl->gas->speciesIndex("C2H2");
    iC2H2   = (iC2H2 >= 0) ? iC2H2 : odtl->gas->speciesIndex("c2h2");

    if(odtl->nspc > 6 && !odtl->odtP->ItableLookup) { 
        yC2H2   = odtl->yspc[iC2H2][ipos];
    }
    else{
        if(!odtl->odtP->ItableLookup) {
            getC2H2(mixf, odtl->enth[ipos], yC2H2, ipos);
        }
        else {
            //getC2H2(odtl->eta[0][ipos], odtl->eta[1][ipos], yC2H2, ipos);

            yC2H2 = odtl->etaTools->luTable->
                        getValAtGridPoint(ipos,odtl->etaTools->luTable->index.species_mass_index+iC2H2) ;
        }
    }

    //---------- scale the acetylene to "account" for reduced local concentrations with soot present

    //double Ycmax = 0.8563 * odtl->mixf[ipos];
    //yC2H2 *= (1-odtl->mom[1][ipos]/Ycmax);



    //get xMassSP
    std::vector<double> yi(odtl->nspc);           // working array

    if(!odtl->odtP->ItableLookup)
        odtl->getYspVecAtPt(ipos, yi);
    
    //Calculate reaction rates for nucleation, growth, oxidation, and coagulation
    double MWC2H2 = 2*odtl->gas->atomicWeight(odtl->gas->elementIndex("C")) 
                    + 2*odtl->gas->atomicWeight(odtl->gas->elementIndex("H"));

    double cC2H2 = odtl->rho[ipos]*yC2H2/MWC2H2;        // kmol/m3
    double Asurf = M_PI*pow(6.0*odtl->rho[ipos]*Ys/M_PI/rho_s,2.0/3.0)*pow(nd,1.0/3.0);

    double r1 = 0.1E5 * exp(-21100.0/odtl->temp[ipos]) * cC2H2;             //kmol/m^3/s 
    double r2 = 0.6E4 * exp(-12100.0/odtl->temp[ipos]) * cC2H2 * sqrt(Asurf); //kmol/m3/s

    double r3;
    if(!odtl->odtP->ItableLookup)
        r3 = 0.1E5  * sqrt(odtl->temp[ipos]) * exp(-19680.0/odtl->temp[ipos]) * 
                (odtl->rho[ipos]*yi[iO2]/odtl->gas->molecularWeight(iO2)) * Asurf;   //kmol/m3*s
    else{
        r3 = 0.1E5  * sqrt(odtl->temp[ipos]) * exp(-19680.0/odtl->temp[ipos]) * 
                (odtl->rho[ipos]*odtl->etaTools->luTable->getValAtGridPoint(ipos, "O2") / 31.9989) * Asurf;
    }
    double r4 = 2.0 * Ca * pow(6.0*MWc/M_PI/rho_s, 1.0/6.0) * sqrt(6.0*kb*odtl->temp[ipos]/rho_s)
                * pow(odtl->rho[ipos]*Ys/MWc, 1.0/6.0) * pow(nd, 11.0/6.0);                 //#/m^3/s 

    //Define number density and soot mass fraction source terms

    mhat0_source   = (2.0*Na*r1/Cmin - r4)/odtl->rho[ipos];         // = S_n/rho  in Lignell dissertation
    mhat1_source   = MWc * (2.0*r1 + 2.0*r2 - r3)/odtl->rho[ipos];  // = S_Ys/rho in Lignell dissertation

    if (!isnan(mhat0_source) && !isnan(mhat1_source) && (mhat0_source > 1e100 || mhat0_source < -1e100))
        cout << "i - " << ipos << " nd - " << mhat0_source << " ys - " << mhat1_source << endl;
    
    //Define source terms for gas phase species

    //    for(int k = 0; k<odtl->nspc; k++)    //set all the unaffected species source terms to 0
    //        gSp_source[k] = 0.0;             //define the terms for affected species
    //gSp_source[iC2H2] = -odtl->gas->molecularWeight(iC2H2) * (r1 + r2);
    //gSp_source[iO2]   = -odtl->gas->molecularWeight(iO2) * r3 / 2;


    //------------ Beji 2008 model:

//    mhat0_source = 0.0;
//    double mixf = odtl->odtP->ItableLookup ? odtl->eta[0][ipos] : odtl->mixf[ipos];
//
//    mhat1_source = 3.5E-5 * odtl->rho[ipos]*odtl->rho[ipos]*(mixf-0.063)/(1.0-0.063)*
//                pow(mixf, 2.2)*exp(-2000/odtl->temp[ipos]) *
//                exp(-5.0*(mixf-0.08)/0.08*(mixf-0.08)/0.08);


    return;
}

///////////////////////////////////////////////////////////////////////////////

/** 
 * 3 moment lognormal using Leung and Lindstead, from dns code soot_m.f90
 */

void MOM_soot::getSootSourceTermsPt3mom(double &mhat0_source,
                                        double &mhat1_source,
                                        double &mhat2_source,
                                        const int ipos) {


    double         M0, M1, M2;
    double         M23, M16, M13, Mm12, Mm16, M76, M43, M56, M53, M12;

    double         b0, b2, Kf, mmin;
    double         f1;//, f2, f3; // unused variables
    double         kk;

    vector<double> Coa(3); 
    vector<double> Nuc(3);
    vector<double> Grw(3);

    //const double   exp1   = 1.0/3.0; // unused variable
    const double   exp2   = 1.0/6.0;
    //const double   exp3   = 11.0/6.0; // unused variable
    const double   exp4   = 2.0/3.0;
    //const double   exp5   = 3.0/2.0; // unused variable
    const double   pi     = 3.141592654;
    const double   mw_cs  = 12.01115;           // carbon molecular weight kg/kmol

    double         temp = odtl->temp[ipos];

    //int            i,j,k,L; // unused variables
    
    //---------------------------

    vector<double> soot_rr(4,0.0);
    get_soot_rr(soot_rr, odtl->temp[ipos], ipos);

    //---------------------------

    M0 = odtl->mom[0][ipos] * odtl->rho[ipos];      // ndens
    M1 = odtl->mom[1][ipos] * odtl->rho[ipos];      // rho*Ys
    M2 = odtl->mom[2][ipos] * odtl->rho[ipos];      // 

    //M0 = 5.0E16;      // doldb
    //M1 = 3.0E-5;      // doldb
    //M2 = 5.0E-27;     // doldb

    //---------------------------

    if(M0 > 0.0 && M1 > 0.0 && M2 > 0.0) {
        kk = 2.0/3.0;    M23  = pow(M0, 0.5*kk*kk-1.5*kk+1) * pow(M1,2*kk-kk*kk) * pow(M2, 0.5*kk*kk-0.5*kk);
        kk = 5.0/3.0;    M53  = pow(M0, 0.5*kk*kk-1.5*kk+1) * pow(M1,2*kk-kk*kk) * pow(M2, 0.5*kk*kk-0.5*kk);
        kk = 1.0/6.0;    M16  = pow(M0, 0.5*kk*kk-1.5*kk+1) * pow(M1,2*kk-kk*kk) * pow(M2, 0.5*kk*kk-0.5*kk);
        kk = 1.0/3.0;    M13  = pow(M0, 0.5*kk*kk-1.5*kk+1) * pow(M1,2*kk-kk*kk) * pow(M2, 0.5*kk*kk-0.5*kk);
        kk = -0.5;       Mm12 = pow(M0, 0.5*kk*kk-1.5*kk+1) * pow(M1,2*kk-kk*kk) * pow(M2, 0.5*kk*kk-0.5*kk);
        kk = -1.0/6.0;   Mm16 = pow(M0, 0.5*kk*kk-1.5*kk+1) * pow(M1,2*kk-kk*kk) * pow(M2, 0.5*kk*kk-0.5*kk);
        kk = 7.0/6.0;    M76  = pow(M0, 0.5*kk*kk-1.5*kk+1) * pow(M1,2*kk-kk*kk) * pow(M2, 0.5*kk*kk-0.5*kk);
        kk = 4.0/3.0;    M43  = pow(M0, 0.5*kk*kk-1.5*kk+1) * pow(M1,2*kk-kk*kk) * pow(M2, 0.5*kk*kk-0.5*kk);
        kk = 5.0/6.0;    M56  = pow(M0, 0.5*kk*kk-1.5*kk+1) * pow(M1,2*kk-kk*kk) * pow(M2, 0.5*kk*kk-0.5*kk);
        kk = 1.0/2.0;    M12  = pow(M0, 0.5*kk*kk-1.5*kk+1) * pow(M1,2*kk-kk*kk) * pow(M2, 0.5*kk*kk-0.5*kk);
    }
    else {
        M23 = 0.0;
        M53 = 0.0;
        M16 = 0.0;
        M13 = 0.0;
        Mm12= 0.0;
        Mm16= 0.0;
        M76 = 0.0;
        M43 = 0.0;
        M56 = 0.0;
        M12 = 0.0;
    }

    //Kf = sqrt(kb*temp/rho_s)*pow(0.75/pi/rho_s,exp2);
    Kf = sqrt(kb*temp/120./rho_s)*pow(0.75/pi/rho_s,exp2);    // 120 is DNS bug, should be temp, not temp/120, but include 120 for consistency with DNS (doesn't matter at short times anyway since coagulation timescales too long)
    b0 = 0.85;
    b2 = 0.85;

    //------ Nucleation

    mmin = mw_cs*Cmin/Na;                    
    Nuc[0] = soot_rr[0]*Na/Cmin/mw_cs;
    Nuc[1] = Nuc[0] * mmin;
    Nuc[2] = Nuc[0] * mmin*mmin;

    //------  Growth

    f1     = pi*pow(6.0/rho_s/pi, exp4);           // geometric factor for SA
    Grw[0] = 0.0;
    Grw[1] = f1 * (soot_rr[1] - soot_rr[2]) * M23;
    Grw[2] = f1 * (soot_rr[1] - soot_rr[2]) * M53 * 2.0;

    //------ Coagulation

    Coa[0] = -Kf*b0  * ( M0*M16 + 2*M13*Mm16 + M23*Mm12 );
    Coa[1] = 0.0;
    Coa[2] = 2*Kf*b2 * ( M1*M76 + 2*M43*M56  + M53*M12  );

    //Coa[0] = Coa[1] = Coa[2] = 0.0; //doldb

    //------ Moment sources

    mhat0_source = ( Nuc[0] + Grw[0] + Coa[0] ) / odtl->rho[ipos];
    mhat1_source = ( Nuc[1] + Grw[1] + Coa[1] ) / odtl->rho[ipos];
    mhat2_source = ( Nuc[2] + Grw[2] + Coa[2] ) / odtl->rho[ipos];

//    cout << endl << odtl->pos[ipos] << " " << mhat0_source*odtl->rho[ipos] << " "  //doldb
//                                           << mhat1_source*odtl->rho[ipos] << " "  //doldb
//                                           << mhat2_source*odtl->rho[ipos] << " " 
//                                           << soot_rr[0] << " " 
//                                           << soot_rr[1] << " " 
//                                           << soot_rr[2] << " " 
//                                           << Nuc[0] << " " 
//                                           << Nuc[1] << " " 
//                                           << Nuc[2] << " " 
//                                           << Grw[0] << " " 
//                                           << Grw[1] << " " 
//                                           << Grw[2] << " " 
//                                           << Coa[0] << " " 
//                                           << Coa[1] << " " 
//                                           << Coa[2] << " " ;
//

}

///////////////////////////////////////////////////////////////////////////////

/**
 *-------------------------------------------------------------------------------------
 *  Compute the reaction rates for soot reactions (kmol/m3*s)
 *  Reaction:
 *			1. Nucleation	        C2H2   ---R1---> 2C(s) + H2
 *			2. Growth	    nC(s) + C2H2   ---R2---> (n+2)C(s) + H2
 *			3. Oxidation     C(s) + 0.5 O2 ---R3---> CO
 *			4. Coagulation          n/Na   ---R5---> n/Na - 1
 *  Note: Everything is done dimensionally
 *-------------------------------------------------------------------------------------
 */

void MOM_soot::get_soot_rr(vector<double> &soot_rr, double temp, int ipos) {

    //int i; // unused variable

    //double rr_ys_conv, rr_yn_conv; // unused variables

    const double  exp1   = 1.0/3.0;
    const double  exp2   = 1.0/6.0;
    const double  exp3   = 11.0/6.0;
    const double  exp4   = 2.0/3.0;
    const double  pi     = 3.141592654;
    const double  mw_cs  = 12.01115;          // carbon molecular weight kg/kmol

    double M0 = odtl->mom[0][ipos] * odtl->rho[ipos];      // ndens
    double M1 = odtl->mom[1][ipos] * odtl->rho[ipos];      // rho*Ys
    //double M2 = odtl->mom[2][ipos] * odtl->rho[ipos];      // unused variable

    //M0 = 5.0E16;      // doldb
    //M1 = 3.0E-5;      // doldb
    //M2 = 5.0E-27;     // doldb

    
    //---------- Concentrations

    int iC2H2, iO2;  
    iC2H2     = odtl->gas->speciesIndex("C2H2");
    iC2H2     = (iC2H2 >= 0) ? iC2H2 : odtl->gas->speciesIndex("c2h2");
    iO2       = odtl->gas->speciesIndex("O2");
    iO2       = (iO2 >= 0) ? iO2 : odtl->gas->speciesIndex("o2");
    if(iC2H2 <=0 || iO2 <=0) {
        cout << "\n\n***************** ERROR: missing C2H2 or O2 in MOM_soot::get_soot_rr " << endl;
        exit(0);
    }
    double c_c2h2 = odtl->yspc[iC2H2][ipos] * odtl->rho[ipos] / odtl->gas->molecularWeight(iC2H2); // kmol/m3
    double c_o2   = odtl->yspc[iO2][ipos]   * odtl->rho[ipos] / odtl->gas->molecularWeight(iO2);   // kmol/m3

    //--------------------------------------------------------------------
    //Leung and Lindstedt rates 
    //--------------------------------------------------------------------

    // Nucleation rate: kg/m3*s    (kg of carbon/m3*s)
   
    soot_rr[0] = 0.1E5*exp(-21100./temp) * c_c2h2 * 2.0 * mw_cs;                           
   
    // Growth Rate:     kg/m2*s
   
    soot_rr[1] = 0.6E4*exp(-12100./temp) * c_c2h2 *
                       2.0 * mw_cs *       // 1/sqrt(Surf area: m2/m3)
                       pow( pi*pow(6.0/pi/rho_s*(abs(M1)+1.1766E-10), exp4) * 
                       pow(abs(M0)+1.E-5,exp1), -0.5 );

    // Oxidation Rate:   kg/m2*s

    soot_rr[2] = 0.1E5*sqrt(temp)*exp(-19680./temp) * c_o2 *mw_cs;

    // Coagulation Rate:  #/m3*s
   
    soot_rr[3] = 18.0*pow(abs(6.0*M1/pi/rho_s),exp2) *
                 sqrt(6.0*kb*temp/rho_s) *
                 pow(abs(M0), exp3);

}

///////////////////////////////////////////////////////////////////////////////

/** Reset var size for mixture fraction and mom_source
 */

void MOM_soot::resizeVar(){

    //resize needed variables
    odtl->setMixfVec();
    //odtl->hLoss.resize(odtl->ngrd);
    for(int i=0; i<odtl->nmom; i++){
       mom_source[i].resize(odtl->ngrd);
    }

    for(int i=0; i<odtl->ngrd; i++)
        for(int j=0; j<odtl->nmom; j++)
            if(odtl->mom[j][i] < 0.0)
                odtl->mom[j][i] = 0.0;

}
///////////////////////////////////////////////////////////////////////////////

/** Get acetylene mass fraction from mixture fraction and heat loss from
* look up table
 * @param mf input: mixture fraction at pt given to calculate heatloss and Y_C2H2.
 * @param h input: enthalpy at pt given to calculate heatloss and Y_C2H2.
 * @param YC2H2 output: Mass fraction of Acetylene calculated with mf and h.
 */

void MOM_soot::getC2H2(const double mf,const double h, double &YC2H2, int posi) {
   int iu, id;         //indecies in table for mixf just greater and less than mf 
   int iup, idn;       //indecies in table for heat loss just greater and less than hloss 
   double hsens;            //Sensible heat for given mf
   double had;              //Adiabatic enthalpy for given mf
   double hloss;            //Heat loss used in table
   double var_up, var_dn;   //interpolation variables

    getIndices(mf,mixF,iu,id);
    //Calculation of Heat loss at mf
    //Sensible heat
    if(mf <= 0.001)                //find mixture fraction indecies
        hloss = 0.0;
    else{
        if(iu == id)
            hsens = Y_C2H2_Table[1][iu];
        else
            interpOnePt(Y_C2H2_Table[0][iu],Y_C2H2_Table[1][iu],Y_C2H2_Table[0][id],
                        Y_C2H2_Table[1][id],mf,hsens);
        //Adiabatic heat
        had   = odtl->strm->h0 * (1-mf) + odtl->strm->h1 * mf;
        hloss = (had - h) / hsens;
    }

    if (hloss >1){
        hloss = 0.99;
    }

    //odtl->hLoss[posi] = hloss;
    getIndices(hloss,hl,iup,idn);
    iup += 2;
    idn += 2;

    if(hloss < hl[0]){                //find YC2H2
        //heat loss lower than chart
        //extrapolate
        interpOnePt(hl[1],Y_C2H2_Table[iup][iu],hl[0],Y_C2H2_Table[idn][iu],hloss, var_up);
            if(var_up < 0)
                var_up = 0.0;
        if(iu == id){ //if mixture fraction is 0 or 1
            YC2H2 = var_up;
        }
        else{ //otherwise extrapolate with lower mixture fraction
            interpOnePt(hl[1],Y_C2H2_Table[iup][id],hl[0],Y_C2H2_Table[idn][id],hloss, var_dn);
            if(var_dn < 0)
                var_dn = 0.0;
            //then interpolate with correct mixture fraction
        interpOnePt(Y_C2H2_Table[0][iu],var_up,Y_C2H2_Table[0][id],var_dn,mf, YC2H2);
        }
    }
    else if(hloss == hl[0]){                //find YC2H2
        //heat loss is equal to lowest on chart
        if (iu == id)
            YC2H2 = Y_C2H2_Table[2][iu];
        else{
            //then interpolate with correct mixture fraction
            interpOnePt(Y_C2H2_Table[0][iu],Y_C2H2_Table[2][iu],Y_C2H2_Table[0][id],
                        Y_C2H2_Table[2][id],mf, YC2H2);
          }
    }
    else if(hloss == hl[nhl-1]){                //find mixture fraction indecies
        //heat loss is equal to highest on chart
        if (iu == id)
            YC2H2 = Y_C2H2_Table[nhl+1][iu];
        else{
            //then interpolate with correct mixture fraction
            interpOnePt(Y_C2H2_Table[0][iu],Y_C2H2_Table[nhl+1][iu],Y_C2H2_Table[0][id],
                        Y_C2H2_Table[nhl+1][id],mf, YC2H2);
        }
    }
    else if(hloss > hl[nhl-1]){                //find mixture fraction indecies
        //interpolate with YC2H2 = 0 at hl =1.000
        interpOnePt(1.000,0.000,hl[nhl-1],Y_C2H2_Table[iup][iu],
                    hloss, var_up);
        if(var_up < 0.0)
            var_up = 0.0;
        if(iu == id){ //if mixture fraction is 0 or 1
            YC2H2 = var_up;
        }
        else{ //otherwise interpolate with lower mixture fraction
            interpOnePt(1.000,0.000,hl[nhl-1],Y_C2H2_Table[iup][id],
                    hloss, var_dn);
            if(var_dn < 0.0)
                var_dn = 0.0;
            //then interpolate with correct mixture fraction
            interpOnePt(Y_C2H2_Table[0][iu],var_up,Y_C2H2_Table[0][id],var_dn,mf, YC2H2);
        }

    }
    else{
        //heat loss is in chart range
        //Y_C2H2 at hloss and mf(iu)
       interpOnePt(hl[iup-2],Y_C2H2_Table[iup][iu],hl[idn-2],Y_C2H2_Table[idn][iu],
                  hloss, var_up);
       if (iu == id)       
           YC2H2 = var_up;
       else {
           //Y_C2H2 at hloss and mf(id)
            interpOnePt(hl[iup-2],Y_C2H2_Table[iup][id],hl[idn-2],Y_C2H2_Table[idn][id],
                        hloss, var_dn);
                                       
           //Interpolation over mf using var_up and var_dn
            interpOnePt(Y_C2H2_Table[0][iu],var_up,Y_C2H2_Table[0][id],var_dn,mf, YC2H2);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

/**Read Y_C2H2 look up table from a file.
 *
 * @param fname input: file name to write.
 */
void MOM_soot::readYc2h2Table(std::string fname) {
     cout << endl << "Reading Y_C2H2 look up table " << fname << endl;

     ifstream ifile(fname.c_str());
     if(!ifile) {
         cout << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
         exit(0);
     }

     string       s1;
     stringstream ss1;

     getline(ifile, s1);            //read line "# nmixf = 200"
     ss1.str(s1);
     ss1 >> s1 >> s1 >> s1 >> nmixf;
     ss1.clear();
     getline(ifile, s1);            //read line "# nhl = 15"
     ss1.str(s1);
     ss1 >> s1 >> s1 >> s1 >> nhl;
     ss1.clear();
     getline(ifile, s1);            //read line "# row has h1 value"

     hl.resize(nhl);

     for(int i=0; i<nhl; i++){       //get heat loss values
         ifile >> hl[i];
     }

     //getline(ifile, s1);            //read line of column headings

    mixF.resize(nmixf);
    Y_C2H2_Table.resize(nhl+2);
    for(int k=0; k<(nhl+2); k++)
        Y_C2H2_Table[k].resize(nmixf);
    
    for(int i=0; i<nmixf; i++){     //loops down column
        for(int k=0; k<(nhl + 2); k++) //loops across the row
            ifile >> Y_C2H2_Table[k][i];
    }
    
    for(int i=0; i<nmixf; i++){     //loops down column
        mixF[i] = Y_C2H2_Table[0][i];
    }
}
///////////////////////////////////////////////////////////////////////////////

void MOM_soot::interpOnePt(double x1, 
                           double y1,
                           double x2, 
                           double y2, 
                           double actualx, 
                           double &ansy) {
    ansy = y1 + (actualx - x1)*(y2 - y1)/(x2 - x1);

}

///////////////////////////////////////////////////////////////////////////////


/** Searches a vector to find the indices of the values that the given
 *      value is between
 *
 *  @param value    input: The value to be searched for
 *  @param vec      input: The vector which the value correlates with
 *  @return an indices object containing the upper and lower bound indices
 */

 void MOM_soot::getIndices(double value, std::vector<double> vec, int &up, int &dn){

   //find the indices

    if (value <= vec[0]){                        //First case, Below the lowest given Value

        up = dn = 0;

    } 
    else if (value >= vec[vec.size() - 1]){ //Second Case, Above the Highest Given Value

        up = dn = vec.size() - 1;

    } 
    else {                                        //Last Case, Between given values
        double dVal;            //Difference between value and the first element in the vector
        double di;              //Change between value in indices
        double i;               //How many steps to go from first value to the desired value           
        //Calculate values
        dVal = value - vec[0];
        di   = vec[1] - vec[0];
        i = dVal / di;

        // indices are above and below the found value
        dn  = (int)floor(i);
        up  = dn + 1;
               //Verify Values (Close values)
        if(vec[dn] > value){
                   --dn;
                   --up;
               }
               if(vec[up] < value){
                   ++up;
                   ++dn;
               }

        }
}

///////////////////////////////////////////////////////////////////////////////


/** Set soot volume fraction
 */

 void MOM_soot::setVolumeFraction(std::vector<double> &volFrac){
     
     if(odtl->odtP->Imom == 1) 
         for(int i=0; i<odtl->ngrd; i++)
             volFrac[i] = odtl->rho[i] * odtl->mom[0][i] / rho_s;
     else
         for(int i=0; i<odtl->ngrd; i++)
            volFrac[i] = odtl->rho[i] * odtl->mom[1][i] / rho_s;

 }







