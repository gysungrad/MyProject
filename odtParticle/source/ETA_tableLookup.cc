/**
 * @file ETA_tableLookup.cc
 * @author Ryan Hintze and Liz Monson
 * Source file for class ETA_tableLookup
 */

#include "odtline.h"
#include "table.h"
#include "diffuser.h"
#include "ETA_tableLookup.h"
#include <math.h>
#include <assert.h>
#include <cstdlib>

using namespace std;

extern processor proc;

///////////////////////////////////////////////////////////////////////////////

/** Constructor used normally by the TableReader class
 * @param _odtl     \input Pointer to the odtline
 * @param luTable_  \input Pointer to table
 */

ETA_tableLookup::ETA_tableLookup(odtline *_odtl, table *luTable_) {

    luTable = luTable_;
    odtl = _odtl;
}

///////////////////////////////////////////////////////////////////////////////


/** Update fluxes at cell faces
 *
 *  Note: This requires a more complete lookup table. The following must be included
 *  - rho
 *  - Dmixf
 *  - lambda
 *  - species in mole fraction
 *  - kabs (absorption coefficient)
 *  
 *  @param dd \input undocumented
 */

void ETA_tableLookup::computeFluxes(vector<double> &dd){

    int i;
    int im;
    //------------- compute the diffusivities
    
    vector<vector<double> > Deta_f(odtl->odtP->Ieta);
    for(int i=0; i<odtl->odtP->Ieta; i++)
        Deta_f[i].resize(odtl->ngrd);

    for(int i =0; i < odtl->ngrd; i++) {
        Deta_f[0][i] = luTable->getValAtGridPoint(i, luTable->index.Dmixf_index);
//        Deta_f[0][i] = odtl->tran->thermalConductivity()/(odtl->rho[i]*luTable->getValAtGridPoint(i, luTable->index.cp_index));
        Deta_f[1][i] = Deta_f[0][i];
        /***********NOTE: diffusion coefficient for enthalpy is equal to the diffusion
        coefficient of the mixture fraction because of unity Lewis number.***********/
        
        if(odtl->odtP->Ieta == 3)
            Deta_f[2][i]  = Deta_f[0][i];
    }
    //////////////////// boundary faces ////////////////////


    if (odtl->odtP->bcType == odtParam::BCTYPE_PERIODIC) { ////////////////// periodic

        i = 0;
        im = odtl->ngrd - 1;

        //todo in source function  + etaSource[k][i]/getValAtGridPoint(i,index.rho_index)

        for (int k = 0; k < odtl->odtP->Ieta; k++) {

            diff->flxProp[diff->iptEta+k][i]          = -dd[i] * Deta_f[k][i] *odtl->rho[i] * (odtl->eta[k][i] - odtl->eta[k][im]);
            diff->flxProp[diff->iptEta+k][odtl->ngrd] = diff->flxProp[k][i];
        }

    }
    
    else if (odtl->odtP->bcType == odtParam::BCTYPE_WALL    ||
             odtl->odtP->bcType == odtParam::BCTYPE_OUTFLOW ||
             odtl->odtP->bcType == odtParam::BCTYPE_INLET   ||
             odtl->odtP->bcType == odtParam::BCTYPE_WALLOUTFLOW) { 

        for (int k = 0; k < odtl->odtP->Ieta; k++) {
            diff->flxProp[diff->iptEta+k][0] = 0.0;
            diff->flxProp[diff->iptEta+k][odtl->ngrd] = 0.0;
        }

    }

    else { //---------------- unknown bc

        *proc.ostrm << endl << "****** ERROR IN DIFFUSER, BCTYPE = " 
                    << odtl->odtP->bcType << " unknown" << endl;
        exit(0);

    }

    //////////////////// interior faces ////////////////////

    for (int k = 0; k < odtl->odtP->Ieta; k++) {
        for (i = 1, im = 0; i < odtl->ngrd; i++, im++) {
            diff->flxProp[diff->iptEta+k][i] = -dd[i] * odtl->rho[i] * Deta_f[k][i] * (odtl->eta[k][i] - odtl->eta[k][im]);
        }
    }

}


///////////////////////////////////////////////////////////////////////////////

void ETA_tableLookup::computeSourceTerms(){
    //-------Add radiative source term 
    if(odtl->odtP->Iradiation) {
        for (int i = 0; i < odtl->ngrd; i++){ 
            diff->rhsSrc[diff->iptEta + 1][i] += diff->radSource_G[i] / odtl->rho[i];
        }
    }
    //------Add particle source term
    if(odtl->odtP->Iparticles) {
        for (int i = 0; i < odtl->ngrd; i++) {
            for (int k = 0; k < odtl->neta; k++) {
                diff->rhsSrc[diff->iptEta + k][i] += (*(diff->gSource))[diff->part->iPtEta + k][i]; //add the source term for eta
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

/** Updates all the odtl->yspc values as well as odtl->rho and odtl->temp by using
 *  the table lookup based on the mixture fraction and heat loss.
 * @param updateRho \input the flag that indicates if density should be updated in line using 
 * lookup table. If no argument is provided to the method, the default value is false. 
 * Density (rho) in line is updated only when updating grid in the diffusion.
 */


void ETA_tableLookup::updateOdtLineVecs(bool updateRho){

    double  mf;                     //The mixture fraction at the current position
    double  h;                      //The enthalpy at the current position
    table::indices mfi;             //Indices of the mixture fractions bounding mf
    table::indices hli;             //Indices of the mixture fraction bounding the hloss
    double  hloss;                  //Heat loss used in table
    vector<double> mixf = luTable->mixf;
    vector<double> hlosses = luTable->hlosses;

    // --------------------------- Resize key odt variables.

    if (updateRho)
        odtl->rho.resize(odtl->ngrd);
    odtl->molec.resize(odtl->ngrd);
    odtl->temp.resize(odtl->ngrd);
    odtl->lambda.resize(odtl->ngrd);
    odtl->heatLoss.resize(odtl->ngrd,0.0);
    
    //---------------------------
    
    for(int i=0; i < odtl->ngrd; i++) {

        h     = odtl->eta[1][i];
        mf    = odtl->eta[0][i];
        if(mf < 0.0){
            *proc.ostrm << endl << "Mixture fraction is negative (Caused probably by adaptMesh::adaptEddyRegionMesh())."
                    "Correcting from "<< mf << " to 0.0 "<< endl;
            odtl->eta[0][i] = 0.0;
            mf = odtl->eta[0][i];
        }
        mfi   = luTable->getIndices(mf, mixf);
        hloss = luTable->calcHLossFromIndices(mf, h, mfi);
        hli   = luTable->getIndices(hloss, hlosses);

        //OK Now we have both the heat loss and mixf,
        //All that's left to do is interpolate all the
        //values between the 4 hloss,mixf values
        //Update rho and molec

        if(luTable->index.rho_index < 0)
            cout << "INVALID RHO INDEX" << endl;
        if(updateRho){
                odtl->rho[i]   = luTable->bilinearInterpolate
                                 (mfi.d, mfi.u, hli.d, hli.u, mf, hloss, luTable->index.rho_index);
        }
        
        odtl->heatLoss[i] = hloss;
        if(luTable->index.visc_index < 0)
            cout << "INVALID VISCOSITY INDEX" << endl;
        odtl->molec[i] = luTable->bilinearInterpolate
                (mfi.d, mfi.u, hli.d, hli.u, mf, hloss, luTable->index.visc_index);
        
        if(luTable->index.temp_index < 0)
            cout << "INVALID TEMP INDEX" << endl;
        odtl->temp[i]  = luTable->bilinearInterpolate
                (mfi.d, mfi.u, hli.d, hli.u, mf, hloss, luTable->index.temp_index);
        
        if(luTable->index.lambda_index < 0)
            cout << endl << "INVALID LAMBDA INDEX" << endl;
        odtl->lambda[i] = luTable->bilinearInterpolate
                (mfi.d, mfi.u, hli.d, hli.u, mf, hloss, luTable->index.lambda_index);
        
    }
}

void ETA_tableLookup::updateOdtLineVecs(int & i,bool updateRho){

    double  mf;                     //The mixture fraction at the current position
    double  h;                      //The enthalpy at the current position
    table::indices mfi;             //Indices of the mixture fractions bounding mf
    table::indices hli;             //Indices of the mixture fraction bounding the hloss
    double  hloss;                  //Heat loss used in table
    vector<double> mixf = luTable->mixf;
    vector<double> hlosses = luTable->hlosses;

    // --------------------------- Resize key odt variables (make sure it is resized just once).

    if(i==0){
        if (updateRho)
            odtl->rho.resize(odtl->ngrd);
        odtl->molec.resize(odtl->ngrd);
        odtl->temp.resize(odtl->ngrd);
    }

    //---------------------------

        h     = odtl->eta[1][i];
        mf    = odtl->eta[0][i];
        if (mf < 0.0) {
            *proc.ostrm << endl << "Mixture fraction is negative (Caused probably by adaptMesh::adaptEddyRegionMesh())."
                    "Correcting from " << mf << " to 0 " << endl;
            odtl->eta[0][i] = 0.0;
            mf = odtl->eta[0][i];
        }
        mfi   = luTable->getIndices(mf, mixf);
        hloss = luTable->calcHLossFromIndices(mf, h, mfi);
        hli   = luTable->getIndices(hloss, hlosses);

        //OK Now we have both the heat loss and mixf,
        //All that's left to do is interpolate all the
        //values between the 4 hloss,mixf values
        //Update rho and molec

        if(luTable->index.rho_index < 0)
            cout << "INVALID RHO INDEX" << endl;
        if(updateRho)
           odtl->rho[i]   = luTable->bilinearInterpolate
                            (mfi.d, mfi.u, hli.d, hli.u, mf, hloss, luTable->index.rho_index);
        
        if(luTable->index.visc_index < 0)
            cout << "INVALID VISCOSITY INDEX" << endl;
        odtl->molec[i] = luTable->bilinearInterpolate
                (mfi.d, mfi.u, hli.d, hli.u, mf, hloss, luTable->index.visc_index);
        
        if(luTable->index.temp_index < 0)
            cout << "INVALID TEMP INDEX" << endl;
        odtl->temp[i]  = luTable->bilinearInterpolate
                (mfi.d, mfi.u, hli.d, hli.u, mf, hloss, luTable->index.temp_index);
        
    
}
