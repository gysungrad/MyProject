/**
 * @file ETA_Aq.cc
 * Source file for class ETA_Aq
 */
#include "ETA_Aq.h"
#include "processor.h"
#include <cmath>
#include "QR.h"
#include <iomanip>
#include<ctime>
#include<limits>
#include<cstdlib>

using namespace std;

extern processor proc;

///////////////////////////////////////////////////////////////////////////////

/** Constructor function
 *
 * 
 */

ETA_Aq::ETA_Aq(odtline *odtl){ 

    Llewis = false; // False = Differential Diffusion ON
    nElem = 6;  
    nEqSp = 10;
    nAqSp = 12;
    nEqSolnVars = 1;
    nEqRxns = 6;
    

    aqSpNames = vector<string>(nAqSp,"");
    aqSpNames[0] = "H2O";
    aqSpNames[1] = "H+";
    aqSpNames[2] = "OH-";
    aqSpNames[3] = "Ca++";
    aqSpNames[4] = "H2CO3*";
    aqSpNames[5] = "HCO3-";
    aqSpNames[6] = "CO3--";
    aqSpNames[7] = "CaHCO3+";
    aqSpNames[8] = "CaCO3o";
    aqSpNames[9] = "CaOH+";
    aqSpNames[10] = "Na+";
    aqSpNames[11] = "Cl-";
    
    z_charges = vector<int>(nAqSp);
    z_charges[0] = 0;
    z_charges[1] = 1;
    z_charges[2] = -1;
    z_charges[3] = 2;
    z_charges[4] = 0;
    z_charges[5] = -1;
    z_charges[6] = -2;
    z_charges[7] = 1;
    z_charges[8] = 0;
    z_charges[9] = 1;
    z_charges[10]= 1;
    z_charges[11]= -1;

    gamma = vector<double>(nAqSp,1.0);
    Klog10 = vector<double>(nEqRxns,0.0);
    Conc   = vector<double>(nAqSp,0.0);

    vector<vector<double> > A(nEqSp, vector<double>(nEqSp, 0.0));   // dummy assigned to A_rxnCoef below
    A[0][1] =  1; A[0][4] = -1; A[0][5] =  1;
    A[1][1] =  1; A[1][5] = -1; A[1][6] =  1;
    A[2][3] = -1; A[2][5] = -1; A[2][7] =  1;
    A[3][3] = -1; A[3][6] = -1; A[3][8] =  1;
    A[4][2] = -1; A[4][3] = -1; A[4][9] =  1;
    A[5][1] =  1; A[5][2] =  1; 
    A[6][0] =  1; 
    A[7][1] =  1; 
    A[8][6] =  1; 
    A[9][3] =  1; 

    Q_rxnCoef = vector<vector<double> > (nEqSp,vector<double>(nEqSp,0.0));
    J = vector<vector<double> > (4,vector<double>(4,0.0));
    A_rxnCoef = A;
    R_rxnCoef = A;

    QRdecomp(R_rxnCoef, Q_rxnCoef);

    mwel.resize(6,0.0);
    mwel[0] = 12.011;  // C
    mwel[1] = 1.0079;  // H
    mwel[2] = 16;      // O
    mwel[3] = 40.078;  // Ca  // lingel had it at 40.079
    mwel[4] = 22.99;   // Na
    mwel[5] = 35.453;  // Cl

    mwsp.resize(nAqSp,0.0);
    mwsp[0]  = mwel[1]*2 + mwel[2]; 
    mwsp[1]  = mwel[1];
    mwsp[2]  = mwel[1] + mwel[2];
    mwsp[3]  = mwel[3];
    mwsp[4]  = mwel[1]*2 + mwel[0] + mwel[2]*3;
    mwsp[5]  = mwel[1]   + mwel[0] + mwel[2]*3;
    mwsp[6]  =             mwel[0] + mwel[2]*3;
    mwsp[7]  = mwel[1]   + mwel[0] + mwel[2]*3 + mwel[3];
    mwsp[8]  =           + mwel[0] + mwel[2]*3 + mwel[3];
    mwsp[9]  = mwel[1]   +           mwel[2]   + mwel[3];
    mwsp[10] = mwel[4];
    mwsp[11] = mwel[5];

    moleElSp = vector<vector<int> >(nAqSp, vector<int>(nElem,0));
    moleElSp[0][1]  = 2; moleElSp[0][2] = 1;
    moleElSp[1][1]  = 1;
    moleElSp[2][1]  = 1; moleElSp[2][2] = 1;
    moleElSp[3][3]  = 1;
    moleElSp[4][1]  = 2; moleElSp[4][0] = 1; moleElSp[4][2] = 3;
    moleElSp[5][1]  = 1; moleElSp[5][0] = 1; moleElSp[5][2] = 3;
    moleElSp[6][0]  = 1; moleElSp[6][2] = 3; 
    moleElSp[7][1]  = 1; moleElSp[7][0] = 1; moleElSp[7][2] = 3; moleElSp[7][3] = 1;
    moleElSp[8][0]  = 1; moleElSp[8][2] = 3; moleElSp[8][3] = 1;
    moleElSp[9][1]  = 1; moleElSp[9][2] = 1; moleElSp[9][3] = 1;
    moleElSp[10][4] = 1;
    moleElSp[11][5] = 1;

    Yspel = vector<vector<double> >(nAqSp, vector<double>(nElem,0));
    for(int i=0; i<nAqSp; i++)
        for(int k=0; k<nElem; k++)
            Yspel[i][k] = moleElSp[i][k]*mwel[k]/mwsp[i];

    cp = 4180;      
    xOld.resize(nEqSolnVars,0.0);
    xOldi.resize(nEqSolnVars,0.0);

        xOldi[0] = 1.0E-11;                  // molar concentration of H+
     //   xOldi[1] = .0666757;                      // molar concentration of Ca++

      // Conc[1]=1e-11;
      // Conc[3]=.0666757;
    
    b  = vector<double>(nEqSp,0.0);
    ISc = vector<double>(7,1.0);  // Coefficients for gamma fit
    ISc[0] =    0.5345;
    ISc[1] =   -1.2157;
    ISc[2] =   -2.7348;
    ISc[3] =  -14.933 ;
    ISc[4] =  284.74  ;
    ISc[5] = -779.45  ;
    ISc[6] =  636.77  ;


   if (odtl->neta<=6)
       cout << "\n********ERROR NOT enought eta slots allocated for species****** \n";


   vector<vector<double> > eta(odtl->neta, vector<double> (odtl->ngrd,0.0));

    setEta(odtl->mixf,eta  );
    for(int i=0; i<odtl->ngrd; i++)
        for(int k=0; k<odtl->neta; k++)
                         odtl->eta [k][i] = eta[k][i];




}

///////////////////////////////////////////////////////////////////////////////

/** Set the initial condition for the eta variable (elemental mass fractions).
 *
 *  @param eta input odtline object
 */

void ETA_Aq::setEta( vector<double> &mixf, vector<vector<double> > &eta  ) {



    int nssp=7;         // number of stream species

    //------------ set streams

    T0 = 273.15 +25;
    h0 = 00.0;
    double h1 = h0+cp*(0.0);

    vector<double> y0(nssp,0.0);
    vector<double> y1(nssp,0.0);
    vector<double> mwsp2(nssp,0.0);
    mwsp2[0] = 2*mwel[1]+mwel[2];   // Water
    mwsp2[1] = 2*mwel[5]+mwel[3];   // CaCl2
    mwsp2[2] = 2*mwel[4]+mwel[0]+mwel[2]*3; // Na2Co3
    mwsp2[3] =   mwel[1]+mwel[5];           // HCL
    mwsp2[4] =   mwel[1]+mwel[2] + mwel[4] ; // NaOH
    mwsp2[5] =   mwel[4]+mwel[5];            // NaCl
    mwsp2[6] =   mwel[4]+mwel[0]+mwel[2]*3+mwel[1] ; // NaHCO3
 
 
    vector<vector<int> > moleNumbers(nElem, vector<int>(nssp,0.0));
    moleNumbers[0][0] = 0; moleNumbers[0][1] = 0; moleNumbers[0][2] = 1; moleNumbers[0][3] = 0;   moleNumbers[0][4] = 0;  moleNumbers[0][5] = 0;// C in H2O, CaCl2, Na2CO3, HCl, NaOH, and NaCl
    moleNumbers[1][0] = 2; moleNumbers[1][1] = 0; moleNumbers[1][2] = 0; moleNumbers[1][3] = 1;   moleNumbers[1][4] = 1;  moleNumbers[1][5] = 0;// H
    moleNumbers[2][0] = 1; moleNumbers[2][1] = 0; moleNumbers[2][2] = 3; moleNumbers[2][3] = 0;   moleNumbers[2][4] = 1;  moleNumbers[2][5] = 0;// O
    moleNumbers[3][0] = 0; moleNumbers[3][1] = 1; moleNumbers[3][2] = 0; moleNumbers[3][3] = 0;   moleNumbers[3][4] = 0;  moleNumbers[3][5] = 0;// Ca
    moleNumbers[4][0] = 0; moleNumbers[4][1] = 0; moleNumbers[4][2] = 2; moleNumbers[4][3] = 0;   moleNumbers[4][4] = 1;  moleNumbers[4][5] = 1;// Na
    moleNumbers[5][0] = 0; moleNumbers[5][1] = 2; moleNumbers[5][2] = 0; moleNumbers[5][3] = 1;   moleNumbers[5][4] = 0;  moleNumbers[5][5] = 1;// Cl

    moleNumbers[0][6] = 1;
    moleNumbers[1][6] = 1;
    moleNumbers[2][6] = 3;
    moleNumbers[3][6] = 0;
    moleNumbers[4][6] = 1;
    moleNumbers[5][6] = 0;

//------------------------GEBAEUR Directly Concentrations-------------------//
/*
    T0 = 273.15 +24;
    double V1, V2,C_CaCl2, C_Na2CO3,C_NaHCO3,C_NaOH;

    vector<double> tm(6,0.0); // total mass assuming a liter basis
V1 = 0.025 ;  // Volume of input 1 [L]
//-----------pH of about 10----------------//
C_Na2CO3 =4.52*1e-3;  // Concentration of Na2CO3 in input 2 [mol/L]
C_NaHCO3 =5.48*1e-3;  // Concentration of NaHCO3 in input 2 [mol/L]
//-----------------------------------------//
//-----------pH of about 9.0----------------//
//C_Na2CO3 =0.655*1e-3;  // Concentration of Na2CO3 in input 2 [mol/L]
//C_NaHCO3 =9.345*1e-3;  // Concentration of NaHCO3 in input 2 [mol/L]
//-----------------------------------------//
C_NaOH =  0;  // Concentration of Na2CO3 in input 2 [mol/L]

 tm[0]= 1000/18.015  ;   //         Concentration of H2O    [mol/L]
 tm[1]= 1e-10;          //           // Concentration of CaCl2  [mol/L]
 tm[4]= C_NaOH;
 tm[2]= C_Na2CO3;              // Concentration of Na2CO3 [mol/L]
 tm[5]= 0;           //      Concentration of NaHCO3 NaCL [mol/L]
 tm[6]= C_NaHCO3;           //      Concentration of NaHCO3 NaCL [mol/L]
 Total= 0;
 for (int i=0; i<nssp; i++)
     Total+=tm[i]*mwsp2[i];

// cout <<  tm[0]<< "  " << tm[1]<< "  " << tm[2]<< "  " << tm[5] <<  " \n";

    y0[0] = tm[0]/Total*mwsp2[0];       // water
    y0[1] = tm[1]/Total*mwsp2[1];       // CaCl2
    y0[2] = tm[2]/Total*mwsp2[2];        // Na2CO3
    y0[3] = 0;                         // HCL
    y0[4] = tm[4]/Total*mwsp2[4];                         // NaOH
    y0[5]=0;
    y0[6]=  tm[6]/Total*mwsp2[6];
cout << " \n";
//    for (int ik=0; ik<6 ; ik++)
//cout << y0[ik] << "     ";
//cout <<   Total <<"  \n";
*/
//--------------------------------Sawada Stream Variables----------------------//
    double V1, V2,C_CaCl2, C_Na2CO3,C_NaCl_mix, Total_1, Total_2;
 double tm0, tm1;  // Total Mass in stream 0 and 1
    vector<double> tm(nssp,0.0); // total mass assuming a liter basis
//------------------------SAWADA PREMIXED Concentrations-------------------// 


V1 = 0.25 ;  // Volume of input 1 [L]
V2 = 0.62 ;  // Volume of input 2 [L]
C_CaCl2  = 6.7e-2 ;  // Concentration of CaCl2  in input 1 [mol/L]
C_Na2CO3 = 1.3e-2 ;  // Concentration of Na2CO3 in input 2 [mol/L]
C_NaCl_mix = 0.1 ;     // Concentration of NaCl added to mixture [mol/L] to make ionic strength = 0.2 (approximate)

 tm[0]= 1000/18.015  ;   //         Concentration of H2O    [mol/L]
 tm[1]= C_CaCl2 *V1/(V1+V2);                     // Concentration of CaCl2  [mol/L]
 tm[2]= C_Na2CO3*V2/(V1+V2);              // Concentration of Na2CO3 [mol/L]
 tm[5]= C_NaCl_mix;           //      Concentration of NaCL [mol/L]
 tm0= tm[0]*mwsp2[0]+C_CaCl2*mwsp2[1]+tm[5]*mwsp2[5];   //Total mass is = sum( mol(i)*mw(i) )
 tm1= tm[0]*mwsp2[0]+C_Na2CO3*mwsp2[2]+tm[5]*mwsp2[5];  //Total mass is = sum( mol(i)*mw(i) ) 
Total= tm[0]*mwsp2[0] + tm[1]*mwsp2[1] + tm[2]*mwsp2[2] + tm[5]*mwsp2[5]; // Total Mass  Stream 1 + stream 2 + Salt // THIS IS USED LATER IN ETA FUNCTION 
// cout <<  tm[0]<< "  " << tm[1]<< "  " << tm[2]<< "  " << tm[5] <<  " \n";

    y0[0] = tm[0]/Total*mwsp2[0];       // water
    y0[1] = tm[1]/Total*mwsp2[1];       // CaCl2
    y0[2] = tm[2]/Total*mwsp2[2];        // Na2CO3
    y0[5] = tm[5]/Total*mwsp2[5];        // Nacl  // 


//    y0[0]+=-y0[3]-y0[2]-y0[1]-y0[4]; // water
 //-----------------NOT PREMIXED Sawada----------------------//  For 1-D cases, run this in combination with Premixed Sawada.

    // stream Index //
    //  0  -  Water
    //  1  -  CaCl2
    //  2  -  Na2CO3
    //  3  -  HCl
    //  4  -  NaOH
    //  5  -  NaCl

V1 = 0.25 ;  // Volume of input 1 [L]
V2 = 0.62 ;  // Volume of input 2 [L]
C_CaCl2  = 6.7e-2 ;  // Concentration of CaCl2  in input 1 [mol/L]
C_Na2CO3 = 1.3e-2 ;  // Concentration of Na2CO3 in input 2 [mol/L]
C_NaCl_mix = 0.1 ;     // Concentration of NaCl added to mixture [mol/L] to make ionic strength = 0.2 (approximate)

  tm[0]= 1000/18.015  ;   //         Concentration of H2O    [mol/L]
  tm[1]= C_CaCl2 ;                     // Concentration of CaCl2  [mol/L]
  tm[2]= C_Na2CO3;                     // Concentration of Na2CO3 [mol/L]
  tm[5]= C_NaCl_mix;           //      Concentration of NaCL [mol/L]

    tm0= tm[0]*mwsp2[0]+C_CaCl2*mwsp2[1]+tm[5]*mwsp2[5];   //Total mass of species i  = sum( mol(i)*mw(i) )   for stream 1 
    tm1= tm[0]*mwsp2[0]+C_Na2CO3*mwsp2[2]+tm[5]*mwsp2[5];  //Total mass of species i  = sum( mol(i)*mw(i) )   for stream 2 

    Total_1   = tm[0]*mwsp2[0] + tm[1]*mwsp2[1] + 5e-6*mwsp2[2] + tm[5]*mwsp2[5]; // Total Mass  Stream 1 + stream 2 + Salt // THIS IS USED LATER IN ETA FUNCTION 
    Total_2   = tm[0]*mwsp2[0]   + 1e-8*mwsp2[1] + tm[2]*mwsp2[2] + tm[5]*mwsp2[5]; // Total Mass  Stream 1 + stream 2 + Salt // THIS IS USED LATER IN ETA FUNCTION 
    Total     = (Total_1*V1+Total*V2) /(V1+V2);

    y0[0] = tm[0]/Total_1*mwsp2[0];       // water
    y0[1] = tm[1]/Total_1*mwsp2[1];       // CaCl2
    y0[2] = 5e-6/Total_1*mwsp2[2];        // Na2CO3
    y0[5] = tm[5]/Total_1*mwsp2[5];        // Nacl  // 

    y1[0] = tm[0]/Total_2*mwsp2[0];       // water
    y1[1] = 1e-8/Total_2*mwsp2[1];       // CaCl2
    y1[2] = tm[2]/Total_2*mwsp2[2];        // Na2CO3
    y1[5] = tm[5]/Total_2*mwsp2[5];        // Nacl  // 


//    y0[0]+=-y0[3]-y0[2]-y0[1]-y0[4]; // water


//---------------------Rodriguez-Blanco Concentrations--------------//
/*
    y1[0] = 1.0;       // water
    y1[1] = 0.0;        // CaCl2
    y1[2] = 0.106;      // Na2CO3
    y1[3] = 0.00000;        // HCl
    y1[4] = 0.00;        // NaOH

    y0[0] = 1.0;       // water
    y0[1] = 0.1109;       // CaCl2
    y0[2] = 0.000;        // Na2CO3
    y0[3] = 0.00001;        // HCl
    y0[4] = 0.00;        // NaOH
    y1[0]+=-y1[3]-y1[2]-y1[1] -y1[4]; //water


     h1=h0;  // uniform temp
*/
//-----------------------------------------------------------------//
//--------------------OLI -------------------//
/*
    y0[0] = 1.0;       // water
    y0[1] = 0.0099965;      //CaCl2 
   // y0[2] = 0.0127250;       //Na2Co3 
    y0[2] = 0.0130157;       //Na2Co3 

    y0[3] = 0.0;        // HCl
    y0[4] = 0.0;        // NaOH
    y0[5] = .0000000000     ; // NaCl  // found through trial and error using acheivein I=.2

    y0[0]+=-y0[3]-y0[2]-y0[1]-y0[4] -y0[5]; // water

    y1[0] = 1.0;       // water
    y1[1] = 0.0099965;      //CaCl2 
    y1[2] = 0.0130157;       //Na2Co3 

    y1[3] = 0.0;        // HCl
    y1[4] = 0.000000;        // NaOH
    y1[5] = .0000000000     ; // NaCl  // found through trial and error using acheivein I=.2


    y1[0]+=-y1[3]-y1[2]-y1[1] -y1[4] -y1[5]; //water

     h1=h0;  // uniform temp
     */
//--------------------Wolf -------------------//
/*
    y0[0] = 1.0;       // water
    y0[1] = 0.1e-10  ;       //CaCl2 
    y0[2] = 0.1e-10   ;       //Na2Co3 

    y0[3] = 0.0;        // HCl
    y0[4] = 0.0;        // NaOH
    y0[5] = 0.0     ; // NaCl  // foun

    y0[0]+=-y0[3]-y0[2]-y0[1]-y0[4] -y0[5]; // water

    y1[0] = 1.0;       // water
    y1[1] = 0.1e-10   ;       //CaCl2 
    y1[2] = 0.1e-10   ;       //Na2Co3 

    y1[3] = 0.00000;        // HCl
    y1[4] = 0.000000;        // NaOH
    y1[5] = .4090000000     ; // NaCl  // ionic strength of 7
    y1[0]+=-y1[3]-y1[2]-y1[1] -y1[4] -y1[5]; //water

     h1=h0;  // uniform temp
     */
//-----------------------------------------------------//
//-----------------------------------------------------//
/*    
    y0[0] = 0.97;       // water
    y0[1] = 0.015;       // CaCl2
    y0[2] = 0.015;        // Na2CO3

    y1[0] = 0.97;      // water
    y1[1] = 0.015;        // CaCl2
    y1[2] = 0.015;      // Na2CO3
     h1=h0;
*/


    int ngrd = eta[0].size();
    vector<vector<double> > ys(nssp,vector<double>(ngrd,0.0));
    for(int k=0; k<nssp; k++)
        for(int i=0; i<ngrd; i++){
            ys[k][i] = y0[k]*(1-mixf[i]) + y1[k]*(mixf[i]);// Hyperbolic tangent profile
  //          ys[k][i] = y0[k]*((double) ngrd-1-( double)  i)/((double) ngrd-1) + y1[k]*((double) i)/((double) ngrd-1);// Linear Profile
    //             ys[k][i] = y0[k];  // UNIFORM profile for DEBUGGINg  or 0-D
/*
if (i>=2*ngrd/3){
           ys[k][i] = y0[k]*((double) ngrd/3-1-( double)  i+200.0)/((double) ngrd/3-1) + y1[k]*((double) i-200.0)/((double) ngrd/3-1);// 3 linear profiles for WOLF
continue;
}
if (i>=ngrd/3){
           ys[k][i] = y0[k]*((double) ngrd/3-1-( double)  i+100.0)/((double) ngrd/3-1) + y1[k]*((double) i-100.0)/((double) ngrd/3-1);// 3 linear profiles for WOLF
continue;
}
else
           ys[k][i] = y0[k]*((double) ngrd/3-1-( double)  i)/((double) ngrd/3-1) + y1[k]*((double) i)/((double) ngrd/3-1);// 3 linear profiles for WOLF
*/
    }
//  0   1  2  3  4  5
//  C   H  O  Ca Na Cl

    //------------ set elemental mass fractions in eta

    for(int j=0; j<ngrd; j++){
        for(int k=0; k<nElem; k++){
            for(int i=0; i<nssp; i++){
                eta[k][j] += ys[i][j]*moleNumbers[k][i]*mwel[k]/mwsp2[i];  // mass fraction of species!!!

          //    eta[i][j]   += moleNumbers[k][i] * mwel[k]/mwsp2[i];


            }
        }
        /*
*/
    }


 //-------------For Comparison of equilibrium to Matlab scripts ----------------//
 /*
    for(int j=0; j<ngrd; j++){
        double mixf;
        mixf=((double)j + .5)/(double)ngrd;// mixture fraction from .005 to .995

 eta[4][j]=  0.0    *(1.-mixf)  +  0.01128 *(mixf);
 eta[5][j]=  0.01278*(1.-mixf)  +  0.0    *(mixf);
 eta[0][j]=  0.00000*(1.-mixf)  +  0.00295*(mixf);
 eta[1][j]=  0.10965*(1.-mixf)  +  0.10898*(mixf);
 eta[2][j]=  0.87035*(1.-mixf)  +  0.87679*(mixf);
 eta[3][j]=  0.00722*(1.-mixf)  +  0.     *(mixf);


    for(int jj=0; jj<6; jj++)
        cout << eta[jj][j] << "  ";
        cout << " \n";

}
*/
 //--------------------------------------------------------------//
/*
    for(int j=0; j<ngrd; j++){
    for(int jj=0; jj<6; jj++)
        cout << eta[jj][j] << "  ";
        cout << " \n";
}
*/

        //--------------Wolf-----------------//
        /*
    T0 = 273.15 ;
    h0=0;
     h1 = h0+cp*(10.0);
    double h2 = h0+cp*(25.0);
    double h3 = h0+cp*(60.0);
    for(int i=0; i<ngrd; i++){
if (i>=2*ngrd/3){
         eta[nElem][i]  = h3;;// 3 linear profiles for WOLF
continue;
}
if (i>=ngrd/3){
         eta[nElem][i]  = h2 ;// 3 linear profiles for WOLF
continue;
}
else
         eta[nElem][i]  = h1;// 3 linear profiles for WOLF
    }
    */
    //-------------------------------------------//


    //------------ enthalpy is last entry of eta
    for(int i=0; i<ngrd; i++)
   //    eta[nElem][i] = h0*((double) ngrd-1-( double)  i)/((double) ngrd-1) + h1*((double) i)/((double) ngrd-1);// Linear Profile
       eta[nElem][i] = mixf[i]*h1 + (1.0-mixf[i])*h0;
      //--------------------------------//



/*
double zero;
vector<double> z2(6,0.0);
z2[0]=4   ; // C
z2[1]= 1  ;  // H
z2[2]=  -2 ; //O
z2[3]=  2 ; // Ca
z2[4]= 1  ; // Na
z2[5]= -1  ; // Cl


zero=0;


    for(int ii=0; ii<ngrd; ii++){
    for(int i=0; i<6; i++)
 zero+=   eta[i][ii]*1000/mwel[i]*z2[i];
cout <<mixf[ii]<<"   " << zero << "  \n";
zero=0;
}

*/


}

///////////////////////////////////////////////////////////////////////////////

/** Computes the activity coefficients gamma on a molar scale from the Davies 
 *  equation.
 *
 *  @param T input, temperature in Kelvin
 *  @param C input, vector of all aqueous species concentrations in moles/L.
 */

void ETA_Aq::setActivityCoefficients(double T, double rho) {

   double eps = 78.46; 
//   eps=32;   This is a fitting parameter the dielectric constant?
 //  eps=80.3;
 //  wolf
//  int datasize=3;


/*
 vector<double> D2(datasize,0.0);
 vector<double> E2(datasize,0.0);
 vector<double> T2(datasize,0.0);
 double D=0.3;

 D2[0] = .0985;     D2[1] = .105;     D2[2] = .100;
 E2[0] = 88.00;   E2[1] = 80.00;  E2[2] = 75.46;
 T2[0] = 283.15; T2[1] = 298.15; T2[2] = 333.15;

if (T>= T2[datasize-1]){
      int i= datasize-2;
      D  = D2[i]+(T-T2[i])*(D2[i+1]-D2[i])/(T2[i+1]-T2[i]);
      eps= E2[i]+(T-T2[i])*(E2[i+1]-E2[i])/(T2[i+1]-T2[i]);
}//extrapolate
else{
 for (int i=0; i<datasize-2; i++)
    if (T <= T2[i+1]){
      D  = D2[i]+(T-T2[i])*(D2[i+1]-D2[i])/(T2[i+1]-T2[i]);
      eps= E2[i]+(T-T2[i])*(E2[i+1]-E2[i])/(T2[i+1]-T2[i]);
      break;
    }

}
*/
  eps = T*(T*.0005959864603  -.7112982350097) +  237.6464504261370; // Temperature dependant dielectric constant or permitivity of space
  double D=.4;

   double A = 1.8248E6 * pow(eps*T,-1.5);


   double I = 0.0;                    // ionic strength
   for(int i=0; i<nAqSp; i++)
       I += 0.5*Conc[i]*(double)z_charges[i]*(double)z_charges[i];

   for(int i=0; i<nAqSp; i++){
       if (z_charges[i] == 0)
           gamma[i] = 1.0; // pow (10.0,.03*I);   // The Value of .03 from Dolejs_2008  This is the extended Debye-Huckel model simplified for neutral species.
       else
       gamma[i] = pow( 10.0, -A * (double) z_charges[i]* (double) z_charges[i] * ( sqrt(I)/  (1.0+sqrt(I)) - D*I ) );
    } 

}

void ETA_Aq::setActivityCoefficients2(double T,double x) {

   for(int i=0; i<nAqSp; i++)
       gamma[i]=1.0;  // Assume ideal for ALL 
/*
I=ISc[6]; // Old mixture fraction guess approach for the ionic strength
for (int j=5; j>=0; j--) 
    I=I*x+ISc[j];
for(int i=0; i<nAqSp; i++)
       gamma[i] = pow( 10.0, -A*z_charges[i]*z_charges[i]*( sqrt(I)/(1.0+sqrt(I)) - 0.3*I ) );
       */
}




///////////////////////////////////////////////////////////////////////////////
/** Compute the RHS function constants for the equilibrium solve
 */
 void ETA_Aq::Compute_RHS_Constants(){  // Only functions of gamma, not even for the 2x2 case 


    b[0] = Klog10[0] / ( gamma[1]*gamma[5]/gamma[4]  ); 
    b[1] = Klog10[1] / ( gamma[1]*gamma[6]/gamma[5]  );
    b[2] = Klog10[2] / ( gamma[7]/gamma[3]/gamma[5]  );
    b[3] = Klog10[3] / ( gamma[8]/gamma[3]/gamma[6]  );
    b[4] = Klog10[4] / ( gamma[9]/gamma[2]/gamma[3]  );
    b[5] = Klog10[5] / ( gamma[1]*gamma[2]           );

    a=1/b[0]/b[1];
    B=(1/b[1]);
    c=(b[2]/b[1]);
    d=(b[3]);
    e=(b[5]);
    f=(b[4]*b[5]);

    return;
 }




///////////////////////////////////////////////////////////////////////////////
/** Compute the RHS function for the equilibrium solve.
Solves nonlinear system of 10 equations and 10 unknowns by writing them
all in terms of H+ concentration and then using a newtonian solve to 
solve the H+ concentration.
 */
void ETA_Aq::equilFrhs(vector<double> &x, vector<double> &y) {
    /*
    //-----------------------------------//
        b[0]-b[5]   - >  equilibrium Constants
             b[6]   - >  H2O   Concentration
             b[7]   - >  H+    Concentration
             b[8]   - >  CO3-- Concentration
             b[9]   - >  Ca++  Concentration
    //-----------------------------------//
    */

     
      
      b[7] =  x[0];  // H+ Concentration
     // b[9] =  x[1] ; // for 2x2 formulation




//    b[8]=Cc/(a*b[7]*b[7]   +B*b[7] + 1.0   + c*b[7]*b[9]   + d*b[9]);  //old 2x2 way  //old Carbon mass balance

//-------- iteratively solves for b[8] gives  nearly the same as analytical, within 1e-16---------------------//
/*

for (int i=0; i<7; i++){
    double b9;
    b9=(1.0   + b[2]/b[1]*b[7]*b[8] +  b[3]*b[8] +   b[4]*b[5]/b[7] );
    b[8]=Cc* b9/(a*b[7]*b[7]* b9  +B*b[7]*b9 + b9)  + c*b[7]*Cca    + d*Cca  );
    cout << b[8] << " guess 2 \n ";
}

*/
// --------------------------Mass balance on Carbon -----------------------//

//--This Carbon balance is simplifeid into a quadratic form with Coefficients A B and C --//

//   zero=  -Cc* (1.0   + b[2]/b[1]*b[7]*b[8] +  b[3]*b[8] +   b[4]*b[5] / 
//   b[7] )+b[8]*(a*b[7]*b[7]* (1.0   + b[2]/b[1]*b[7]*b[8] +  b[3]*b[8] +   b[4]*b[5]/b[7] )  +
//   B*b[7]*(1.0   + b[2]/b[1]*b[7]*b[8] +  b[3]*b[8] +   b[4]*b[5]/b[7] ) + 
//   1.0* (1.0   + b[2]/b[1]*b[7]*b[8] +  b[3]*b[8] +   b[4]*b[5]/b[7] )  + c*b[7]*Cca    + d*Cca  );

//---------------------Mass balance on Carbon simplified-------------------//
    double CoefC, CoefB, CoefA;
    CoefC =-Cc*(1.0+f/b[7]);
    CoefB=-Cc*(c*b[7] +d)  + (a*b[7]*b[7]+B*b[7]+1.0)*(1.0+f/b[7])+c*b[7]*Cca    + d*Cca;
    CoefA=(a*b[7]*b[7]+B*b[7]+1.0)*(c*b[7]+d);
    // use qaudratic formula for roots 
    b[8]=(-CoefB + sqrt(CoefB*CoefB-4.0*CoefA*CoefC))/(2.0*CoefA); 
  //  cout << (-CofB -sqrt(CofB*CofB-4.0*CofA*CofC))/(2.0*CofA) << "    " << b[8] << "  \n";
   // use + , rather than ( -b -( b^2 - 4ac )^.5 ) / 2a  
   // Negative  b[8] generated by negative, maybe always, not verified.  May have multiple physical solutions.

// --------------------------Mass balance on Calcium -----------------------//
    b[9] =   Cca / (1.0 + c*b[7]*b[8] +  d*b[8] + f/b[7]);
// -------------------------Mass balance on Hydrogen -----------------------//
    b[6]=( b[7] + e/b[7]   + 2*a*b[7]*b[7]*b[8] + B*b[7]*b[8]  + c*b[7]*b[8]*b[9] + f/b[7]*b[9]     - Ch ) /-2.0;

    // This constitutes the linear solver, redone analytically
    Conc[0]=b[6];                 //   
    Conc[1]=b[7];
    Conc[2]=b[5]/b[7];            //
    Conc[3]=b[9];                 
    Conc[4]=b[7]*b[7]*b[8]*a;     //  
    Conc[5]=b[7]*b[8]*B;          //
    Conc[6]=b[8];                 //
    Conc[7]=b[9]*b[7]*b[8]*c;     //
    Conc[8]=b[3]*b[9]*b[8];        //
    Conc[9]=f*b[9]/b[7];             //

// -------------------------Mass balance on Oxygen ------------------------//
      y[0] =   Conc[0]   + Conc[2] + 3*Conc[4] + 3*Conc[5] + 3*Conc[6] + 3*Conc[7] + 3*Conc[8] + Conc[9] - Co  ;

 //     y[1] =   Conc[3]   + Conc[7] +   Conc[8] +   Conc[9]                                               - Cca;  // for 2x2

}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
/** Compute the RHS function for the equilibrium solve used in comparison with WOLF.
//ASSUMES in equilibrium with solid phase DOES NOT CONSERVE mass of O, Ca, or C!!!!
 */

void ETA_Aq::equilFrhsw(vector<double> &x, vector<double> &y) {

      b[7] =  x[0] ;
      b[9] =  x[1] ;



double co2;
co2=922./101325.*pow(10.0,108.3865 +.01985076*TA-6919.53/TA-40.45154*log10(TA)+669365/TA/TA)/gamma[4];
    
//double kc2=pow(10.0,108.3865 +.01985076*TA-6919.53/TA-40.45154*log10(TA)+669365/TA/TA); //  !!!!!     currently unused variable
double kc=pow(10.,-171.9065 - 0.077993  *TA + 2839.319/TA + 71.595  *log10(TA));
b[6]=( b[7] + e/b[7]   + 2*a*b[7]*b[7]*b[8] + B*b[7]*b[8]  + c*b[7]*b[8]*b[9] + f/b[7]*b[9]     - Ch ) /-2.; // using this because it reports how much water is in the system on a hydrogen balance in the system
b[8]=co2/b[7]/b[7]*b[1]*b[0];
b[9]=kc/gamma[6]/gamma[3]/b[8];

Conc[0]=b[6];
Conc[1]=b[7];
Conc[2]=b[5]/b[7];
Conc[3]=b[9];
Conc[4]=co2;
Conc[5]=b[7]*b[8]/b[1];
Conc[6]=b[8];
Conc[7]=b[2]*b[9]/b[1]*b[7]*b[8];
Conc[8]=b[3]*b[9]*b[8];
Conc[9]=b[4]*b[9]*b[5]/b[7];

y[1]=kc/gamma[6]/gamma[3]/b[8]-b[9];





y[0]=0;
for (int i=1; i<12; i++)
   y[0]+=Conc[i]*z_charges[i]; 
/*

   cout <<x[0]<<"  "  <<y[0] << "    \n";
for (int i=1; i<12; i++)
 cout <<   aqSpNames[i]<< "       ";
   cout << "\n";
for (int i=1; i<12; i++)
   cout << Conc[i] << "  ";
   cout << "\n";
int delme;
cin >> delme;
*/
}
/** Compute the Jacobian for the equilibrium solve
 */

void ETA_Aq::equilJac(vector<double> x1, vector<vector<double> > &Jac) {
 // Could be optimized for the 1D case
  vector<double> F1(nEqSolnVars,0.0);
  vector<double> F2(nEqSolnVars,0.0);
 // vector<double> F3(nEqSolnVars,0.0);
  vector<double> x2(nEqSolnVars,0.0);
 // vector<double> x3(nEqSolnVars,0.0);
  double h;
  equilFrhs(x1,F1);
  for(int j=0; j<nEqSolnVars; j++) {
    
    h = x1[j]*1.0E-4;
    x2 = x1;
    x2[j] = x2[j] + h;
    equilFrhs(x2,F2);

    for(int i=0; i<nEqSolnVars; i++){
        Jac[i][j] = (F2[i]- F1[i])  / h;
    //    Jac[i][j] = (F2[i]- F3[i])  / (h*2.0);
  

}

  }

}

///////////////////////////////////////////////////////////////////////////////

/** Set log10 (Keq)
 */

void ETA_Aq::setKlog10(double T) {

    Klog10[0] = pow(10., -356.3094 - 0.06091964*T + 21834.37/T + 126.8339*log10(T) - 1684915. /T/T);
    Klog10[1] = pow(10., -107.8871 - 0.03252849*T + 5151.79 /T + 38.92561*log10(T) - 563713.9/T/T);
    Klog10[2] = pow(10., 1209.120 + 0.31294   *T - 34765.05/T - 478.782 *log10(T));
    Klog10[3] = pow(10.,-1228.732 - 0.299444  *T + 35512.75/T + 485.818 *log10(T));
    Klog10[4] = pow(10., 2.272282                - 260.0688/T);
    Klog10[5] = pow(10.,6.0875   - 0.01706   *T - 4470.99 /T);
 //   cout<< Klog10[3] << " \n ";
//  int   delme;
//  cin >> delme;

}

///////////////////////////////////////////////////////////////////////////////

/** Do liquid equilibrium, sets vector Conc
 * @param Yelem input elemental mass fractions
 * @param rho   input density (kg/m3)
 */

void ETA_Aq::equilibrateAll(vector<vector<double> > &eta, double rho ) {

    int ngrd   = eta[0].size();
    CaCO3.resize(  ngrd  );
    CaCO3_gamma.resize(  ngrd  );
    vector<double> Yelem(nElem,0.0);
    vector<double> Temp(ngrd, 298.15);
    getTempVecFromH(eta, Temp);
    allConc = vector<vector<double> >(nAqSp, vector<double>(ngrd,0.0));

    for(int i=0; i<ngrd; i++) {

        for(int k=0; k<nElem; k++)
            Yelem[k] = eta[k][i];

            TA=Temp[i]; /// Derek
        if(i==0){
            equilibrate(Yelem, rho, Temp[i], true);
        }
        else {
           equilibrate(Yelem, rho, Temp[i], false);
        }

        for(int k=0; k<nAqSp; k++)
            allConc[k][i] = Conc[k];

        CaCO3[i]=Conc[8];
        CaCO3_gamma[i]=gamma[8];
        
    }


}


///////////////////////////////////////////////////////////////////////////////
/** Do liquid equilibrium, sets vector Conc
// Mass is Conserved to machine precision for all elements, except for possibly oxygen.
Oxygen is also to machine precision for small tolerences, however the 
convergence criteria does not work for small tolerences <1e-10.  This can be improved.
 * @param Yelem input elemental mass fractions
 * @param rho   input density (kg/m3)
 * @param T     input temperature (K)
 */

bool ETA_Aq::equilibrate(vector<double> &Yelem, double rho, double T, bool Lfirst) {



    // Correct_rho(rho,Yelem[4], Yelem[5]);

    int    maxit  = 100;
    int    maxit2 = 50;
    double tol    = 1.0E-6;
    double tol2   = 1.0E-15;  
    //double normo = 0.0; //  !!!!!     currently unused variable
    vector<vector<double> > Jac(nEqSolnVars, vector<double>(nEqSolnVars,0.0));
    vector<double> F0(nEqSolnVars,0.0);
    vector<double> dx(nEqSolnVars,0.0);
    vector<double> x(nEqSolnVars,0.0);

    vector<double> xOld2(nEqSolnVars,0.0);
    setKlog10(T);
    elemMassFrac_to_elemConc(Yelem,rho,T);  // get elemental concentrations


    if(Lfirst) {
        xOld[0]=1e-13;
    }
    else {

        xOld[0] =    Conc[1];

    }
    double norm;
    for (int it2=1; it2<=maxit2; it2++){

        //--------------- Newton solve

        vector<double> gamma_old = gamma;
        double relax = 1.0;
        int id=0;

        Compute_RHS_Constants();

        for(int it=1; it<=maxit+1; it++) {    

            equilFrhs(xOld, F0);

            for(int i=0; i<nEqSolnVars; i++){
                F0[i] *= -1.0;
            }


            equilJac(xOld, Jac);
            //  QR(Jac, F0, dx);       // For 2x2 system
            //      if (Jac[0][0]==0){
            //          xOld[0]=1e-7;
            //          continue;
            //     }

            dx[0]= F0[0]/Jac[0][0]; // use for 1D solvers othwise call QR


            id=0;  // This variable is the index in the array of varioubles solved for.
            while (id<nEqSolnVars){
                x[id] = xOld[id]+relax*dx[id];
                if (x[id]<0 ){  //will Generate NAN if this "if" statement isn't included 19-9-2012
                    id+=-1;
                    relax=relax*.8;
                }
                id+=1 ;
            }


            if(it>maxit){ 
                *proc.ostrm << "\n#*********** WARNING NO EQUILIBRIUM CONVERGENCE IN " << maxit-1 << " ITERATIONS" << "  " << Cc<< "  "  <<Ch<< "  "  << Co<< " " << Cca<< "  " << Cna<< "  " << Ccl <<"  "<<x[0] << "  " << F0[0]<<" "  << endl;
                return false;
            }
            //------------old  Convergence check
            norm = 0.0;
            for(int i=0; i<nEqSolnVars; i++)
                norm += abs(dx[i]/relax/x[i]); // This is a sum and should be dived by nEqSolnVars for expanson
            if(norm < tol){
                break;
            }

            //------------ Update x

            relax=1.0;
            xOld = x;


        }   // Inner loop, solve for optimal H+


        gamma_old=gamma;
        setActivityCoefficients( T,  rho); 
        norm=0.0;
        for(int i=0; i<nEqSp; i++)
            norm += abs(gamma[i]-gamma_old[i])/gamma[i]; // This is a sum and should be dived by nEqSolnVars for expanson


        if(norm < tol2){
            break;
        }



    } // Outer loop converge Activity Coefficients


    return true;
}









///////////////////////////////////////////////////////////////////////////////
/*
//----------------Density is a function of ionic strength---------------//
// this function assumes that salt is the primary contributor to density changes
// and uses a small 1-D lookup table to approximate the density of the fluid
*/
void ETA_Aq::Correct_rho(double &rho, double Yna,double Ycl){
    int datasize=5;
    vector<double> rhoI(datasize,0.0);
    vector<double> salI(datasize,0.0);
    vector<double> IS(datasize,0.0);
    double apI =0.0;

    rhoI[0]=1000.0 ;rhoI[1]=1039.6 ;rhoI[2]=1080 ;rhoI[3]=1165.628 ;rhoI[4]=1257.9 ;
    salI[0]=0      ;salI[1]=50.0   ;salI[2]=100.0;salI[3]=200.0    ;salI[4]=300.0  ; 
    // data fro http://www.csgnetwork.com   Water density calculator;
    for (int i=0; i<datasize; i++)
        IS[i]= salI[i]/(mwel[4]+mwel[5])*rhoI[i]*.001;


    for (int ii=0; ii<3; ii++){

        Cna      = Yna*rho/mwel[4];    // Na
        Ccl      = Ycl*rho/mwel[5];    // Cl
        apI= Cna+Ccl;

        for (int i=0; i<datasize; i++){
            if (apI > IS[i]){
                rho= rhoI[i]+(apI-IS[i])*(rhoI[i+1]-rhoI[i])/(IS[i+1]-IS[i]);
                return;
            }
            if (i==datasize-1)
                rho= rhoI[i]+(apI-IS[i])*(rhoI[i]-rhoI[i-1])/(IS[i]-IS[i-1]);
        }
    }
    return;
}

///////////////////////////////////////////////////////////////////////////////

/** Convert elemental mass fractions to concentration
Mass fraction notation
-----------------------------
Summ of species + sum of third moments of each polymorph = 1.0
 */

void ETA_Aq::elemMassFrac_to_elemConc(vector<double> &Yelem, double rho, double T ){

    //    double sum=0.0;
    //    for(int i=0; i<Yelem.size(); i++)
    //        sum += Yelem[i];               
    //cout << sum << "\n";
    //    for(int i=0; i<Yelem.size(); i++)   // Equilibrate normalizes elements here!   // DO NOT normalize for moment cases. 
    //        Yelem[i] /= sum;               // should I still normalize with the moment term included?
    rho=Total;  // To be consistent with sean.

    //rho=1.008963000528736e3;

    Cc       = Yelem[0]*rho/mwel[0];    // C
    Ch       = Yelem[1]*rho/mwel[1];    // H
    Co       = Yelem[2]*rho/mwel[2];    // O
    Cca      = Yelem[3]*rho/mwel[3];    // Ca
    Cna      = Yelem[4]*rho/mwel[4];    // Na
    Ccl      = Yelem[5]*rho/mwel[5];    // Cl
    Conc[10] = Cna; 
    Conc[11] = Ccl;


    /*

       Cc       = 0.009264367816092        ;    // C
       Ch       =   1.110185956147655e+02        ;    // H
       Co       =   55.537090910831012     ;    // O
       Cca      =   0.019252873563218    ;    // Ca
       Cna      =   0.118528735632184     ;    // Na
       Ccl      =   0.138505747126437       ;    // Cl
       Conc[10] = Cna; 
       Conc[11] = Ccl;
     */
    //Ch=111.02;
    //Co=55.537;
    // Correct_rho(rho,Yelem[4], Yelem[5]);
    //    cout << rho << "  ";
    //cout << " \n";
    //     cout <<setprecision(17) << " \n";
    //     cout << Cc << "   " << Ch   << "   " << Co   << "   " <<  Cca  << "   " <<Cna   << "   " <<Ccl   << "  \n";
    //    int delme;
    //    cin >> delme;
    /*
       Cc       = Yelem[0]*rho/mwel[0];    // C
       Ch       = Yelem[1]*rho/mwel[1];    // H
       Co       = Yelem[2]*rho/mwel[2];    // O
       Cca      = Yelem[3]*rho/mwel[3];    // Ca
       Cna      = Yelem[4]*rho/mwel[4];    // Na
       Ccl      = Yelem[5]*rho/mwel[5];    // Cl
       Conc[10] = Cna; 
       Conc[11] = Ccl;
     */

    //int delme;
    //cin >> delme;
}

///////////////////////////////////////////////////////////////////////////////

/** Fill temperature vector from enthalpy vector
 *  @param eta input vector of vectors of eta [ieta][pos] 
 *  @param T output vector of temperatures 
 */
void ETA_Aq::getTempVecFromH(vector<vector<double> > &eta, vector<double> &T) {

    int ngrd = eta[0].size();
    int ienth = nElem;            // enth index 
    for(int i=0; i<ngrd; i++)
        T[i] = getTgivenH(eta[ienth][i]);
    
}

///////////////////////////////////////////////////////////////////////////////

/** Get temperature from enthalpy 
 *  (note, enthalpy is with enthalpy of formation=0 at T0)
 * @param h input enthalpy (J/kg)
 * @return tempeature (K)
 */

double ETA_Aq::getTgivenH(double h) {

    return T0 + (h-h0)/cp;  

}

///////////////////////////////////////////////////////////////////////////////

/** Get species diffusivities
 * @param T input tempeature (K)
 * @param Di output diffusivities (m2/s) 
 */

void ETA_Aq::getDiffusivities(vector<double> &T, vector<vector<double> > &Di) {

    int ngrd = T.size();

    for(int i=0; i<ngrd; i++) {

        Di[0][i]  = 2.2108E-9;     // H2O
        Di[1][i]  = 8.27258E-9;    // H+
        Di[2][i]  = 4.72208E-9;    // OH-
        Di[3][i]  = 7.50157E-10;   // Ca++
        Di[4][i]  = 1.84956E-9;    // H2CO3* = CO2 (aq) and H2CO3, value taken as CO2 (aq)
        Di[5][i]  = 1.10511E-9;    // HCO3-
        Di[6][i]  = 9.21507E-10;   // CO3-
        Di[7][i]  = 6.85158E-10;   // CaHCO3+
        Di[8][i]  = 6.62867E-10;   // CaCO3o
        Di[9][i]  = 7.47279E-10;   // CaOH+
        Di[10][i] = 1.24578E-9;    // Na+
        Di[11][i] = 1.87421E-9;    // Cl-

        for(int j=0; j<nAqSp; j++)
            Di[j][i] = Di[j][i] * T[i] / 298.15 * 0.00089 / getViscosity(T[i]);

    }

}

///////////////////////////////////////////////////////////////////////////////

/** Get Viscosity
 * @param T input tempeature (K)
 * @return viscosity Pa*s
 */

double ETA_Aq::getViscosity(double T) {
    return 3.2509E12*pow(T, -6.28614);
}

///////////////////////////////////////////////////////////////////////////////

/** Get Thermal Conductivity
 * @param T input tempeature (K)
 * @return viscosity Pa*s
 */

double ETA_Aq::getThermalConductivity(double T) {
    return 8.68732E-3*pow(T,0.745028);
}

///////////////////////////////////////////////////////////////////////////////

/** Get diffusive fluxes of elemental mass fractions
 * Using the Nernst-Planck Equation, with an assumed ideal fluid
 * @param eta    input transported elements and enthalpy
 * @param rho    input density
 * @param dx     input cell sizes
 * @param flxEta output fluxes
 */

void ETA_Aq::getEtaDiffusiveFluxes(vector<vector<double> > &eta, double rho,
                                vector<double> &dx,
                                vector<vector<double> > &flxEta,
                                bool Lmom) {
    //Derek:  Remove variable bool Ladaptive from this function
    int ngrd = eta[0].size();
    int neta = eta.size();

    vector<double> T(ngrd,0.0);
    getTempVecFromH(eta, T);

/*    double one=0;            // MASS BALANCE CHECK
   for ( int i=0; i<eta[0].size(); i++){
   for ( int j=0; j<eta.size(); j++)
      one +=  eta[j][i];
      if (one <(1.0-1e-6))
      cout << one << "  SMALL  \n";
      if (one >(1.0 +1e-6))
      cout << one << "   LARGE\n";
      one =0;
   }
   */
    //////////////////////// For a simple scaling approach (constant Le)
    
    if(Llewis) {
        if (Lmom)
            equilibrateAll(eta, rho);


        double Sc = 500;
        double Pr = 1;

        int i, im;
        for(i=1, im=0; i<ngrd; i++, im++) {

            double dd = 2.0*getViscosity(T[i])/(dx[im]+dx[i]) / Sc;

            for(int k=0; k<neta-1; k++) 
                flxEta[k][i]  = -dd*(eta[k][i]-eta[k][im]);  // Yelem
            for(int k=0; k<neta-1; k++)  // FLUXES ZERO  DEBUG
                flxEta[k][i]  = 0.0;  // Yelem


            dd = 2.0*getViscosity(T[i])/(dx[im]+dx[i]) / Pr;
            flxEta[neta-1][i] = -dd*(eta[neta-1][i]-eta[neta-1][im]);  // enthalpy
        }

        for(int k=0; k<neta; k++) {            // boundary fluxes
            flxEta[k][0]    = 0.0;
            flxEta[k][ngrd] = 0.0;
        }
        //----------------//
     // equilibrateAll(eta, rho);
//---------------------------// to be deleted, after confirmed that this isn't needed.

        return;
    }
    equilibrateAll(eta, rho);
    //////////////////////// For a more general approach using the Nernst-Planck equation

    vector<vector<double> > Di(nAqSp, vector<double>(ngrd,0.0));
    getDiffusivities(T, Di);

    vector<vector<double> > Ji(nAqSp, vector<double>(ngrd+1,0.0));
   // vector<vector<double> > Ji2(nAqSp, vector<double>(ngrd+1,0.0));
    int i, im;

    // Derek's Version------------------------------------
    //Maxwell_Stefan_Diffusion(dx,Ji2,T);
    // Currently configured to output reformulated M-S diffusion with Non-idealities included could use clean up
    //--------- compute species fluxes (kg/m2*s)

    for(i=1, im=0; i<ngrd; i++, im++) {        // interior cells

        double dd = 2.0/(dx[im]+dx[i]);
        double Cz2D_sum    = 0;
        double DgradCz_sum = 0;
        double s1 = 0.0;

        for(int k=0; k<nAqSp; k++) {
            Cz2D_sum    += allConc[k][i]*z_charges[k]*z_charges[k]*Di[k][i];  
            DgradCz_sum += Di[k][i] * (allConc[k][i]-allConc[k][im]) * dd * z_charges[k];
        }

        for(int k=0; k<nAqSp; k++) {
            Ji[k][i] = ( -Di[k][i]*(allConc[k][i]-allConc[k][im])*dd + 
                    allConc[k][i]*z_charges[k]*Di[k][i]/Cz2D_sum*DgradCz_sum ) * mwsp[k];
            s1 += Ji[k][i];


        }


        Ji[0][i] -= s1;   // sum j_i = 0.0, enforced through j_h2o (which is neutral!)

        //----------- Now convert to elemental mass fluxes

        for(int k=0; k<nElem; k++) {
            flxEta[k][i] = 0.0;
            for(int j=0; j<nAqSp; j++)
                flxEta[k][i] += Ji[j][i] * Yspel[j][k];
        } 
        //---------- Do enthalpy, which is assumed only a function of temperature

        flxEta[neta-1][i] = -getThermalConductivity(T[i])*(T[i]-T[im])*dd;
    }
    for(int k=0; k<neta; k++) {            // boundary fluxes
        flxEta[k][0]    = 0.0;
        flxEta[k][ngrd] = 0.0;
    }


return;

}


void ETA_Aq::Analytical_Jac( vector<double> x, double T){  // 4x4 formulation

vector<double> Keq2(6,0.0);
vector<double> b(6,0.0);

    Keq2[0] = pow(10,-356.3094 - 0.06091964*T + 21834.37/T + 126.8339*log10(T) - 1684915 /T/T);
    Keq2[1] = pow(10,-107.8871 - 0.03252849*T + 5151.79 /T + 38.92561*log10(T) - 563713.9/T/T);
    Keq2[2] = pow(10, 1209.120 + 0.31294   *T - 34765.05/T - 478.782 *log10(T));
    Keq2[3] = pow(10, -1228.732 - 0.299444  *T + 35512.75/T + 485.818 *log10(T));
    Keq2[4] = pow(10, 2.272282                - 260.0688/T);
    Keq2[5] = pow(10, 6.0875   - 0.01706   *T - 4470.99 /T);

    b[0] = Keq2[0] /(    gamma[1] * gamma[5] / gamma[4] ) ;
    b[1] = Keq2[1] / (    gamma[1] * gamma[6] / gamma[5] ) ;
    b[2] = Keq2[2] /  (    gamma[7] / gamma[3] / gamma[5] ) ;
    b[3] = Keq2[3] /   (    gamma[8] / gamma[3] / gamma[6] ) ;
    b[4] = Keq2[4] /    (   gamma[9] / gamma[2] / gamma[3])  ;
    b[5] = Keq2[5] /    (   gamma[1] * gamma[2]         )  ;


J[0][0]=0;
J[1][0]=2;
J[2][0]=1;
J[3][0]=0;

J[0][1]=x[2]/(b[1]) + (x[2]*x[3]*(b[2]))/(b[1]) + (2*x[1]*x[2])/((b[0])*(b[1]));
J[1][1]=x[2]/(b[1]) - (b[5])/x[1]/x[1] + (x[2]*x[3]*(b[2]))/(b[1]) - (x[3]*(b[4])*(b[5]))/x[1]/x[1] + (4*x[1]*x[2])/((b[0])*(b[1])) + 1;
J[2][1]=(3*x[2])/(b[1]) - (b[5])/x[1]/x[1] + (3*x[2]*x[3]*(b[2]))/(b[1]) - (x[3]*(b[4])*(b[5]))/x[1]/x[1] + (6*x[1]*x[2])/((b[0])*(b[1]));
J[3][1]=(x[2]*x[3]*(b[2]))/(b[1]) - (x[3]*(b[4])*(b[5]))/x[1]/x[1];
 
J[0][2]=x[3]*(b[3]) + x[1]/(b[1]) + x[1]*x[1]/((b[0])*(b[1])) + (x[1]*x[3]*(b[2]))/(b[1]) + 1;
J[1][2]=x[1]/(b[1]) + (2*x[1]*x[1])/((b[0])*(b[1])) + (x[1]*x[3]*(b[2]))/(b[1]);
J[2][2]=3*x[3]*(b[3]) + (3*x[1])/(b[1]) + (3*x[1]*x[1])/((b[0])*(b[1])) + (3*x[1]*x[3]*(b[2]))/(b[1]) + 3;
J[3][2]=x[3]*(b[3]) + (x[1]*x[3]*(b[2]))/(b[1]);

J[0][3]=x[2]*(b[3]) + (x[1]*x[2]*(b[2]))/(b[1]);
J[1][3]=((b[4])*(b[5]))/x[1] + (x[1]*x[2]*(b[2]))/(b[1]);
J[2][3]=3*x[2]*(b[3]) + ((b[4])*(b[5]))/x[1] + (3*x[1]*x[2]*(b[2]))/(b[1]);
J[3][3]=x[2]*(b[3]) + ((b[4])*(b[5]))/x[1] + (x[1]*x[2]*(b[2]))/(b[1]) + 1;
return;
}



void ETA_Aq::Maxwell_Stefan_Diffusion(vector<double> &dx2,vector<vector<double> > &allJi, vector<double> &Temp){
int ngrd;
ngrd=allConc[0].size();
vector<double> ct(ngrd,0);
vector<double> x(nAqSp,0);
vector<double> w(nAqSp,0);
vector<double> xm(nAqSp,0);
vector<double> d(nAqSp,0);
vector<double> Jn(nAqSp,0);
vector<double> Jn2(nAqSp-1,0);
vector<double> diam(nAqSp,0);
vector<double> xs(nAqSp,0);
double C_Davies=0.4;  // davies coefficient usually .1 .2 .3 or in our case .4
xs[0]=1;
vector<vector<double> > Dij(nAqSp,vector<double>(nAqSp,1e200));
vector<vector<double> > Bij(nAqSp,vector<double>(nAqSp,0));
vector<vector<double> > Di2(nAqSp,vector<double>(ngrd,0));
vector<vector<double> > dx (nAqSp,vector<double>(ngrd,0));  // Delta Mole fraction (not space) That is dx2
vector<double> dc(nAqSp,0);
double Farad= 9.65e4; //C /mol
double dphi = 0; //  Electorical potential gradient
double K_c  = 0 ; //Electrical Conductivity of a mixture
double sum1 = 0;
double dd;


vector<double> sx(nAqSp,0);


diam[0]=1.9e-10;   // water meters
diam[1]=1.8e-10;  // Smaller than water =/  oh-
diam[2]=9e-11;  // Very small h+
diam[3]=3.95e-10;  // Ca
diam[4]=3.94e-10;    // H2CO3
diam[5]=3.93e-10;    // HCO3-
diam[6]=3.92e-10;    //CO3--
diam[7]=6.5e-10;   //CaHCO3+
diam[8]=6.2e-10;   //CaCO3o
diam[9]=4.1e-10;   //CaOH+
diam[10]=3.716e-10;  //Na+
diam[11]=1.9e-10;  //Cl-    // tabulated as 189 but it will be larger than that, because of an extra electron

sx[0]=1;

double I = 0.0;                    // ionic strength
double dG=0;  // Chemical potential gradient
//double dxdg=0;  //  !!!!!    currently unused variable     ?? gradient of each species activity coefficient to a species j
int im;


getDiffusivities(Temp, Di2) ;

////////////////////////For loop should start here
//int delme3; //  !!!!!     unused variable
for (int k=0; k<ngrd ; k++){
for (int i=0; i<nAqSp; i++)
ct[k]+=allConc[i][k];
}
for (int ii=0; ii<ngrd+1; ii++){
    if( ii==0 | ii==ngrd)    //Boundaries  set fluxes to zero
        for (int k=0; k<nAqSp; k++)
            allJi[k][ii]=0;
    else{


im=ii-1;
dd = 2.0/( dx2[im] + dx2[ii] );

for (int i=0; i<nAqSp; i++){
w[i]=allConc[i][ii]*mwsp[i]/1000;
x[i]=allConc[i][ii]/ct[ii];
xm[i]=allConc[i][im]/ct[im];
dx[i][ii]=x[i] - xm[i];
dc[i]=allConc[i][ii]-allConc[i][im];
}
///////////////////////////////////////////////////////////////
////////////////////// NON - IDEAL Case ///////////////////////
///////////////////////////////////////////////////////////////
vector<double> Gamma(nAqSp,1.0);
/*for (int i=0; i<nAqSp; i++)
sum1+=x[i]*z_charges[i]*z_charges[i]; /// This value is 2x I
sum1+=-1.0/(ct[ii]*C_Davies);

double sum2=0.0;

for (int j=0; j<nAqSp-1; j++){
//for (int j=1; j<nAqSp; j++){
sum2+= x[i]*(z_charges[j]*z_charges[j]-z_charges[nAqSp-1]*z_charges[nAqSp-1])/sum1;
//sum2+= x[i]*(z_charges[j]*z_charges[j]-z_charges[0]*z_charges[0])/sum1;
//sum2+= x[i]*(z_charges[j]*z_charges[j])/sum1;
//cout << sum2 << "\n ";
//cin >> sum1;
}
Gamma[i]=sum2+1;
sum2=0;
}

sum1=0;


*/

////////////////////////////?Try two

//double I = 0.0;                    // ionic strength
//for(int i=0; i<nAqSp; i++)
//I += 0.5*Conc[i]*z_charges[i]*z_charges[i];
I=0.0;
for (int i=0; i<nAqSp; i++)
I+=allConc[i][ii]*z_charges[i]*z_charges[i]; /// This value is 2x I
I*=.5;

// Derivative of charges with respect to concentrations? // May need this for the solver some day =/
//cout <<"   I="<< I << "\n";
double eps = 78.46; 
double A = 1.8248E6 * pow(eps*Temp[ii],-1.5);
double sum3;
for (int i=0; i<nAqSp; i++){
for (int j=0; j<nAqSp-1; j++)
sum3+=-x[i]*log(10)*z_charges[i]*z_charges[i]*A*(1/((1/sqrt(I)+1)*(1/sqrt(I)+1))/(I*sqrt(I))*.5-C_Davies)*
(z_charges[j]*z_charges[j]-z_charges[nAqSp-1]*z_charges[nAqSp-1])*.5*ct[ii];

//Gamma[i]=sum3+1;

sum3=0;
}
/*
////
////////////////////////////////Try 3 three
vector<double> xL(nAqSp,0.0);
vector<double> xH(nAqSp,0.0);
vector<double> gL(nAqSp,0.0);
vector<double> g(nAqSp,0.0);
vector<double> gH(nAqSp,0.0);
vector<double> Gamma2(nAqSp,0.0);
double sum4=0;
double dx4=0.0001;
xL=x;
xH=x;
for (int i=0; i<nAqSp; i++){
for (int j=0; j<nAqSp-1; j++){
xH[j]+=dx4*x[j];
xH[nAqSp-1]-=dx4*x[j];
//xH[0]-=dx4*x[j];
xL[j]-=dx4*x[j];
xH[nAqSp-1]+=dx4*x[j];
//xH[0]+=dx4*x[j];

setActivityCoefficients3(Temp[ii], ct[ii], xH) ;
g=gamma;
setActivityCoefficients3(Temp[ii], ct[ii], xH) ;
gH=gamma;
setActivityCoefficients3(Temp[ii], ct[ii], xL) ;
gL=gamma;

sum4+=x[i]*(gH[j]-gL[j])/(xH[j]-xL[j])/g[j];
xL=x;
xH=x;
}

Gamma2[i]=sum4+1;
sum4=0;



}

int delme2;
*/
//for (int i=0; i<nAqSp; i++)
//Gamma[i]=1;
//    cout<<"species " << i<< "  "<<  Gamma[i] << "   " << Gamma2[i] << "    \n";

    //cin>>  delme2;
//////////////////////////////////////////////////////////
////////////// Get M-S Diffusivities /////////////////////
/////////////////////////////////////////////////////////
I=0;
for (int i=0; i<nAqSp; i++)
I+=x[i]*z_charges[i]*z_charges[i]; /// This value is 2x I
I*=.5;
//double D12=0; //  !!!!!     unused variable
//double D21=0; //  !!!!!     unused variable
double mu=0;
mu=getViscosity(Temp[ii]);

//cout <<" \n";
//cout <<" \n";

for (int i=1; i<nAqSp; i++)
for (int j=1; j<nAqSp; j++)
if (z_charges[i] !=0 && z_charges[j] !=0)
Dij[i][j]= (Di2[i][ii] +   Di2[j][ii])/2*pow(I,.55)/pow(abs(z_charges[i]*z_charges[j]*1.0),1);  //For electrolytes  // multiplied by 1.0 because abs doesn't have abs( int x)
else
/*
if (diam[i]/diam[j] >=3)
Dij[i][j]=8.31447*Temp[i]/(2*3.14159*6.022e23*mu*1.5*diam[i]);   /// Liquids

    else if (diam[i]/diam[j] <=.5){
Dij[i][j]=8.31447*Temp[i]/(2*3.14159*6.022e23*mu*0.5*diam[i]);   /// Liquids  

}
else*/
Dij[i][j]=8.31447*Temp[i]/(2*3.14159*6.022e23*mu*diam[i]*diam[i]/diam[j]);// liquids where diameters are about the same

//cout <<" \n";
//cout <<" \n";

for (int j=0; j<nAqSp; j++)
    Dij[0][j]=Di2[j][ii];
//for (int j=0; j<nAqSp; j++)
//    Dij[j][0]=Di2[j][ii];

//cout <<" asdflkasjflaskfjaslkfjaslkdfjaslkfj";
for (int i=1; i<nAqSp; i++)/////// Makes for sure that D_AB = D_BA
for (int j=0; j<nAqSp; j++)
Dij[i][j]=Dij[j][i];
/*
cout << "\n";
for (int j=0; j<nAqSp; j++){
for (int i=0; i<nAqSp; i++)
cout<< Dij[i][j] <<  "   ";
cout<<   "\n ";
}
*/
int delme;
//cin>> delme;
//cout << " ?????????";
double sum2;
sum2=0;
sum1=0;
for (int j=1; j<nAqSp; j++){  // because j=0 is water and we want to skip the solvent  j=0 would work too since water is neutral
         K_c       += allConc[j][ii]*  z_charges[j]*z_charges[j]* Di2[j][ii];  // Units Checked
         sum1      += z_charges[j]* Di2[j][ii] * (dc[j])*dd;            // Units Checked
         sum2      += z_charges[j]*allConc[j][ii];
}



dphi=-sum1/K_c;  // Voltage gradient, per Faradays constant




vector<vector<double> > B2(nAqSp-1,vector<double>(nAqSp,0));
vector<double> d2(nAqSp-1,0);



//cout<< " Dphi" << dphi<< "\n";

for (int i=1; i<nAqSp; i++)   // This loop is for the rows
for (int j=1; j<nAqSp; j++){   // This loop is for the columns
    Bij[i][j]=0;
if (i==j){
    for (int k=0; k<nAqSp; k++){  // This loop is for the sum
    if (k!=i)                
//    Bij[i][j]+=x[k]/Dij[j][k];  //Maxwell Stefan 
//    B2[i-1][j-1]+=xs[k]/Dij[j][k];  //Nernst-Planck 2 
  //  B2[i-1][j-1]+=x[k]/Dij[j][k];  //Maxwwell Stefa 2 
    B2[i-1][j-1]+=x[k]/Dij[j][k];  //Maxwwell Stefan Reformulated
}
  B2[i-1][j-1]+=x[0]/w[0]*w[i]/Dij[i][0];  //Maxwwell Stefan Reformulated
//Bij[i][j]=1/Dij[i][0]; // Nersnst-Planck
//B2[i-1][j-1]=1/Dij[i][0]; // Nersnst-Planck
//B2[i-1][j-1]=1/Di2[i][ii]; // Nersnst-Planck debug
}
else
//Bij[i][j]=-x[i]/Dij[i][j]; // Maxwell Stefan
//B2[i-1][j-1]=xs[i]/Dij[i][j]; // Nernst-Planck 2 
//B2[i-1][j-1]=x[i]/Dij[i][j]; // Maxwell Stefa 2 
B2[i-1][j-1]=-x[i]/Dij[i][j]+x[i]*x[0]/w[0]*w[j]/x[j]; // Maxwell Stefan Reformulated 
//Bij[i][j]=0;    // Nernst-Planck
//B2[i-1][j-1]=0;    // Nernst-Planck
}
//////////display loops
/*
for (int j=0; j<nAqSp-1; j++){
for (int i=0; i<nAqSp-1; i++)
cout<< B2[i][j] <<  "   ";

cout<<   "\n ";
}
*/
////////////display loops
//:comment out things here and up


//for (int k=0; k<nAqSp-1;k++)   // Sum of   Why is this to nAqSP - 1??
for (int i=0; i<nAqSp;i++){  // b vector in x=A\b
//dG+=  (-1+x[k]/gamma[k]*dxdg)*dx[k][ii]*dd;  // non-ideal  sum

//dG=  (-nAqSp+1+2);  // ideal    dG= x[i]/R/T * grad mu
//dG=1;


d[i]= -( Gamma[i]*dc[i]*dd + (dphi*z_charges[i]*allConc[i][ii]));  
//d[i]=-ct[ii]*(dG*dd + (dphi*z_charges[i]*x[i]));  
if (i!=0)
d2[i-1]= -( Gamma[i]*dc[i]*dd + (dphi*z_charges[i]*allConc[i][ii]) -w[i]*sum2*dphi );  
//d2[i-1]=-ct[ii]*( dx[i][ii]*dd + (dphi*z_charges[i]*x[i]));  

}
dG=0;
//QR(Bij, d, Jn);
QR(B2, d2, Jn2);

//d2[i-1]=-ct[ii]*( Gamma[i]*dx[i][ii]*dd + (dphi*z_charges[i]*x[i]));     // effectively the same as the next line, except Ct is factored out. and 1/Di moved to the other side
//for (int i=0; i<nAqSp;i++)  
//  allJi[i][ii]=(-ct[ii]*Dij[0][i]*dd*Gamma[i]*dx[i][ii]-allConc[i][ii]*z_charges[i]*Dij[0][i]*dphi)*mwsp[i];  // Easy S-B
double sum4;
double sum5;
double sum6;
sum4=0;
sum5=0;
sum6=0;
for (int i=1; i<nAqSp;i++){  // b vector in x=A\b
//allJi[i][ii]=Jn[i]*mwsp[i];//  Maxwell stefan


allJi[i][ii]=Jn2[i-1]*mwsp[i];  // Maxwell stefan and Nernst planck!!!
//sum4+= allJi[i][ii];
//sum5+=Jn2[i-1]*w[i]/x[i];


allJi[0][ii]+=Jn2[i-1]*w[i]/x[i];
}
allJi[0][ii]*=-x[0]/w[0];
//sum5+=allJi[0][ii]*w[0]/x[0];
allJi[0][ii]*=mwsp[0];

for (int i=0; i<nAqSp;i++){  // b vector in x=A\b
sum5+=allJi[i][ii];
//cout << allJi[i][ii] << "  " ;


}

//cout << sum5 << "  \n";
//cout << sum5 << "  \n";
//allJi[0][ii] -= sum4;
//sum1=0;

/////////////////////////////////////////////////////////////////
//////////////////////Solve for the current I2///////////////////
/////////////////////////////////////////////////////////////////
double I2;
double I3;
double elect;
double elect2;
I2= -Farad*sum1-Farad*dphi*K_c;

for (int i=0; i<nAqSp;i++) { // b vector in x=A\b
I3+=z_charges[i]*Jn[i];
elect+= z_charges[i]*dc[i]*dd;
elect2+= z_charges[i]*allConc[i][ii];
}

//cout << I2 <<"  "  <<   I3 << "  " << elect       <<" \n ";
//cout << -Farad*sum1 << "   " << -Farad*dphi*K_c<< "   "<< I2<< " \n";

//for (int i =0; i<nAqSp ; i++)
//cout << Jn[i] << "  ";
//cout << I3 << " \n ";

//for (int i =0; i<nAqSp ; i++)
//cout << z_charges[i] << "  ";
//cout << I3 << " \n ";


//for (int i =0; i<nAqSp ; i++)
//cout << z_charges[i]*(dc[i]*dd) << "  ";
//cout << elect << " \n ";
//cin >> delme;
if (abs(elect2)>1e-4 ){
for (int i =0; i<nAqSp ; i++)
cout << z_charges[i]*(allConc[i][ii]) << "  ";
cout << elect2 << " \n ";
cin >> delme;
}

elect2=0;
I2=0;
I3=0;
elect=0;
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

K_c=0;
sum1=0;
sum2=0;
sum3=0;
sum4=0;
sum5=0;
sum6=0;
}}
return;
}



       void ETA_Aq::computeFluxes(vector<double> &dd){ return;} 
       void ETA_Aq::computeSourceTerms() { return;} 
       void ETA_Aq::updateOdtLineVecs(bool updateRho ) { return;}
       void ETA_Aq::updateOdtLineVecs(int & i, bool updateRho ){ return;}

