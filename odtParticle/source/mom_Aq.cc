/**
 * @file mom.cc
 * Source file for class mom
 */
//#include<stdio.h>
#include "mom_Aq.h"
#include "processor.h"
#include <vector>
//#include"pdgen.h"
#include<cmath>
#include<iostream>
using namespace std;

extern processor proc;



// Constructor Function //

mom_Aq::mom_Aq( ETA_Aq *ETAp, odtline *odtlpi)
{


ETApx=ETAp; // this is so DQMOM can use QR decomposition, might be useful in otherways
odtlp=odtlpi;

//----------------------Problem Type---------------------///


    nmom= 4;   // Number of Moments Per polymorph
    npoly =  odtlp->nmom/nmom;  // Number of Polymorphs tracked
    np2=nmom/2;   // Number of Abscissa
    Lostwald=false;
    Ldqmom=true;
 //   Ldqmom=false;
  //Lvanoss=true ;
    Lvanoss=false ;




//----------------------Constants-------------------------//

    kb   = 1.3806503e-23 ;  //Boltzman J/K /# number of atoms
    R    = 8.31447  ;  // J/mol K
    MWCC = 100.09   ;  // g / gmol         Molecular weight of dry Calcium Carbonate
    unit =   1.0e-6 ;  //meter / Micrometers   
    unit2=unit*unit ;
    unit3=unit2*unit;

/*//--------------------------Table----------------------//

               Index  |  Species  
                 0    |  ACC
                 1    |  Vaterite
                 2    |  Aragonite
                 3    |  Calcite



//--------------------------------------------------------//


*/ 

//-----------------PolyMorph Propertie--------------------//

    gamma  = vector<double>(npoly,1.0) ;  
    gamma_o= vector<double>(npoly,1.0) ;  
      mp = vector<vector<double> > (npoly,vector<double>(npoly,-1.0));  // This parameter should be between .1 and .2
    if (npoly>=1){
    //    gamma[0]=.062*unit2;  // J/microm^2  ACC            EMPERICAL VALUE COMPARED WITH GEbaeur  // .062 old  .051
//        gamma[0]=.060*unit2;  // J/microm^2  ACC            EMPERICAL VALUE COMPARED WITH GEbaeur  // .062 old  .051
//        gamma[0]=.058*unit2;
  //  gamma[0]=.06015*unit2;   // Try 1
    gamma_o[0]=.0570*unit2;   // optimium 25 C hetero
    gamma_o[0]=.057*unit2;   // optimium 50 C homo  SEAN
    gamma_o[0]=.057*unit2;   // optimium 50 C homo Derek 
    gamma_o[0]=.05684*unit2;   // optimium 60 C homo Derek 
    gamma_o[0]=.05682*unit2;   // optimium 70 C homo Derek 
    gamma_o[0]=.05500*unit2;   // optimium 25 C homo Derek 
    gamma_o[0]=.05500*unit2;   // optimium 30 C homo Derek 
    gamma_o[0]=.05500*unit2;   // optimium 40 C homo Derek 
    gamma_o[0]=.05300*unit2;   // optimium 10 C homo Derek 
    gamma_o[0]=.05300*unit2;   // optimium 25 C homo Derek TAKE 2  QMOM
    gamma_o[0]=.05300*unit2;   // optimium 30 C homo Derek TAKE 2  QMOM
    gamma_o[0]=.05300*unit2;   // optimium 40 C homo Derek TAKE 2  QMOM
    gamma_o[0]=.05300*unit2;   // optimium 50 C homo Derek TAKE 2  QMOM
    gamma_o[0]=.05300*unit2;   // optimium 60 C homo Derek TAKE 2  QMOM
    gamma_o[0]=.05300*unit2;   // optimium 70 C homo Derek TAKE 2  QMOM
    gamma_o[0]=.05300*unit2;   // optimium 10 C homo Derek TAKE 2  QMOM
    gamma_o[0]=.05300*unit2;   // optimium 25 C homo Derek TAKE 2  QMOM
  //  gamma_o[0]=.0550*unit2;   // GEBAUR
  //  gamma_o[0]=.020*unit2;   // GEBAUR  pH of 9.0
   // gamma_o[0]=.00005*unit2;   // GEBAUR  pH of 9.0

 //         gamma[0]=.0520*unit2;  // Try 2
//        gamma[0]=.062*.95*unit2;
}
    if (npoly>=2){
   //     gamma[1]=.12*unit2; // J/mmicrom^2 .034  Vaterite     9.5e-3 to 108e-3  J/m^2// This is a factor of 3 too small for my program to handle =/
   //     gamma[1]=.1189*unit2; // J/mmicrom^2 .034  Vaterite     9.5e-3 to 108e-3  J/m^2// This is a factor of 3 too small for my program to handle =/
     //     gamma[1]=.1335*unit2; // J/mmicrom^2 .034  Vaterite     9.5e-3 to 108e-3  J/m^2// This is a factor of 3 too small for my program to handle =/
          gamma_o[1]=.139*unit2; // optimium 25 C Hetero
          gamma_o[1]=.1387*unit2; // optimium 50 C Homo  SEAN
          gamma_o[1]=.139*unit2; // optimium 50 C Homo  Derek 
          gamma_o[1]=.1430*unit2; // optimium 60 C Homo  Derek 
          gamma_o[1]=.1485*unit2; // optimium 70 C Homo  Derek 
          gamma_o[1]=.1290*unit2; // optimium 25 C Homo  Derek 
          gamma_o[1]=.1300*unit2; // optimium 30 C Homo  Derek 
          gamma_o[1]=.1340*unit2; // optimium 40 C Homo  Derek 
          gamma_o[1]=.1225*unit2; // optimium 10 C Homo  Derek 
          gamma_o[1]=.1274*unit2; // optimium 25 C Homo  Derek TAKE 2 QMOM
          gamma_o[1]=.1296*unit2; // optimium 30 C Homo  Derek TAKE 2QMOM
          gamma_o[1]=.1323*unit2; // optimium 40 C Homo  Derek TAKE 2QMOM
          gamma_o[1]=.1374*unit2; // optimium 50 C Homo  Derek TAKE 2QMOM
          gamma_o[1]=.1420*unit2; // optimium 60 C Homo  Derek TAKE 2QMOM
          gamma_o[1]=.1470*unit2; // optimium 70 C Homo  Derek TAKE 2QMOM
          gamma_o[1]=.1224*unit2; // optimium 10 C Homo  Derek TAKE 2QMOM
          gamma_o[1]=.12803*unit2; // optimium 25 C Homo  Derek ROUND 3 QMOM
       //   gamma_o[1]=.1296*unit2; // optimium 30 C Homo  Derek  ROUND 3 QMOM
       //   gamma_o[1]=.1329*unit2; // optimium 40 C Homo  Derek  ROUND 3 QMOM
 //         gamma_o[1]=.1370*unit2; // optimium 50 C Homo  Derek  ROUND 3 QMOM 
//         gamma_o[1]=.1426*unit2; // optimium 60 C Homo  Derek  ROUND 3 QMOM
         // gamma_o[1]=.1480*unit2; // optimium 70 C Homo  Derek  ROUND 3 QMOM

 //         gamma_o[1]=.12355*unit2; // optimium 10 C Homo  Derek TAKE 2QMOM

 //         gamma_o[1]=.1225*unit2; // optimium 25 C Homo  Derek DQMOM
 //       gamma[1]=.130*unit2; // try 2 
   //     gamma[1]=.123*unit2; // J/mmicrom^2 .034  Vaterite     9.5e-3 to 108e-3  J/m^2// This is a factor of 3 too small for my program to handle =/
//    gamma_o[1]=.14000*unit2;   // TESTING purposes only for DQMOM and QMOM
      gamma_o[1]=.12360*unit2; //Attempt 4 - 10 C 
      gamma_o[1]=.12800*unit2; //Attempt 4 - 25 C 
    //  gamma_o[1]=.12960*unit2; //Attempt 4 - 30 C 
    //  gamma_o[1]=.13290*unit2; //Attempt 4 - 40 C 
    //  gamma_o[1]=.1371*unit2; //Attempt 4 - 50 C 
    //  gamma_o[1]=.1425*unit2; //Attempt 4 - 60 C 
    //  gamma_o[1]=.1480*unit2; //Attempt 4 - 70 C 

   
}
    if (npoly>=3){
  //      gamma[2]=.145*unit2; // J/mmicrom^2   Aragonite    150e-3      J/m^2
  //    gamma[2]=.1635*unit2;
      gamma_o[2]=.1831*unit2;  // optimium 25 C hetero
      gamma_o[2]=.1616*unit2;  // optimium 50 C homo // SEAN
      gamma_o[2]=.1624*unit2;  // optimium 50 C homo // Derek 
      gamma_o[2]=.1652*unit2;  // optimium 60 C homo // Derek 
      gamma_o[2]=.1671*unit2;  // optimium 70 C homo // Derek 
      gamma_o[2]=.1590*unit2;  // optimium 30 C homo // Derek 
      gamma_o[2]=.1605*unit2;  // optimium 40 C homo // Derek 
      gamma_o[2]=.1565*unit2;  // optimium 10 C homo // Derek 
      gamma_o[2]=.1575*unit2;  // optimium 25 C homo // Derek TAKE 2QMOM
      gamma_o[2]=.1585*unit2;  // optimium 30 C homo // Derek TAKE 2QMOM
      gamma_o[2]=.1586*unit2;  // optimium 40 C homo // Derek TAKE 2QMOM
      gamma_o[2]=.1603*unit2;  // optimium 50 C homo // Derek TAKE 2 QMOM
      gamma_o[2]=.1621*unit2;  // optimium 60 C homo // Derek TAKE 2QMOM
      gamma_o[2]=.16485*unit2;  // optimium 70 C homo // Derek TAKE 2QMOM
      gamma_o[2]=.1560*unit2;  // optimium 10 C homo // Derek TAKE 2QMOM
      //-----------------------
      gamma_o[2]=.1580*unit2;  // optimium 25 C homo // Derek ROUND 3 QMOM
      gamma_o[2]=.1580*unit2;  // optimium 25 C homo // Derek ROUND 3 QMOM
      gamma_o[2]=.1595*unit2;  // optimium 30 C homo // Derek ROUND 3 QMOM
      gamma_o[2]=.1600*unit2;// optimium 40 C Homo  // Derek ROUND 3 QMOM
      gamma_o[2]=.16055*unit2;  // optimium 50 C homo //     ROUND 3 QMOM
      gamma_o[2]=.1626*unit2;  // optimium 60 C homo         ROUND 3 QMOM
      gamma_o[2]=.1653*unit2;  // optimium 70 C homo // De   Round 3 QMOM
      gamma_o[2]=.1560*unit2;  // optimium 10 C homo // Derek TAKE 2QMOM

      gamma_o[2]= .1555*unit2;  //Attempt 4 - 10 C
      gamma_o[2]= .1655*unit2;  //Attempt 4 - 25 C
   //  gamma_o[2]= .1700*unit2;  //Attempt 4 - 30 C
   //  gamma_o[2]= .1600*unit2;  //Attempt 4 - 40 C
   //  gamma_o[2]= .16075*unit2;  //Attempt 4 - 50 C
   //  gamma_o[2]= .1630*unit2;  //Attempt 4 - 60 C
   //  gamma_o[2]= .16525*unit2;  //Attempt 4 - 70 C
    }
    if (npoly>=4){
    //    gamma[3]=.280*unit2; // J/mmicrom^2 .034  Calcite      6.5e-3 to 280e-3 J/m^2  // try .17 for calcite
    //    gamma[3]=.1373*unit2; // J/mmicrom^2 .034  Calcite      6.5e-3 to 280e-3 J/m^2  // try .17 for calcite
 //      gamma[3]=.15625*unit2; // J/mmicrom^2 .034  Calcite      6.5e-3 to 280e-3 J/m^2  // try .17 for calcite
       gamma_o[3]=.1870*unit2; // optimium 25 C hetero 
       gamma_o[3]=.1605*unit2; // optimium 50 C homo_sean 
       gamma_o[3]=.1608*unit2; // optimium 50 C homo_Derek 
       gamma_o[3]=.1654*unit2; // optimium 60 C homo_Derek 
       gamma_o[3]=.1700*unit2; // optimium 70 C homo_Derek 
       gamma_o[3]=.15185*unit2; // optimium 25 C homo_Derek 
       
       gamma_o[3]=.15380*unit2; // optimium 30 C homo_Derek 
       gamma_o[3]=.15855*unit2; // optimium 40 C homo_Derek 
       gamma_o[3]=.144600*unit2; // optimium 10 C homo_Derek 
       gamma_o[3]=.15000*unit2; // optimium 25 C homo_Derek TAKE 2
       gamma_o[3]=.15250*unit2; // optimium 30 C homo_Derek TAKE 2
       gamma_o[3]=.15680*unit2; // optimium 40 C homo_Derek TAKE 2
       gamma_o[3]=.15870*unit2; // optimium 50 C homo_Derek  TAKE 2
       gamma_o[3]=.1620*unit2; // optimium 60 C homo_Derek TAKE 2
       gamma_o[3]=.1665*unit2; // optimium 70 C homo_Derek TAKE 2
       gamma_o[3]=.144520*unit2; // optimium 10 C homo_Derek TAKE 2
       gamma_o[3]=.1509505*unit2; // optimium 25 C homo_Derek DQMOM
     //  gamma_o[3]=.15325*unit2; // optimium 30 C homo_Derek TAKE 2
     //  gamma_o[3]=.15690*unit2; // optimium 40 C ROUND 3 QMOM
//       gamma_o[3]=.15980*unit2; // optimium 50 C ROUND 3 QMOM
 //     gamma_o[3]=.1620*unit2; // optimium 60 C ROUND 3 QMOM
    //   gamma_o[3]=.1657*unit2; // optimium 70 C ROUND 3 QMOM 
//       gamma_o[3]=.144990*unit2; // optimium 10 C homo_Derek TAKE 2

         gamma_o[3]=.145000*unit2; // attempt 4 - 10 C 
         gamma_o[3]=.151000*unit2; // attempt 4 - 25 C 
     //    gamma_o[3]=.153100*unit2; // attempt 4 - 30 C 
     //    gamma_o[3]=.156700*unit2; // attempt 4 - 40 C 
     //    gamma_o[3]=.159800*unit2; // attempt 4 - 50 C 
     //    gamma_o[3]=.162900*unit2; // attempt 4 - 60 C 
     //    gamma_o[3]=.1663600*unit2; // attempt 4 - 70 C 


    
        // Information from Donnet_2009

}


//----------------------------------_Van OSS-------------------------//
    if (npoly==0);
    else{
        ap= vector < double > (3*(npoly -1),0.0);
        bp= vector < double > (3*(npoly -1),0.0);

        gp_lw= vector < double > ((npoly -1),0.0);
        gp_p= vector < double > ((npoly -1),0.0);
        gp_m= vector < double > ((npoly -1),0.0);

        ap[0]= 0.213829257334127*1e-3*unit2;   bp[0] = 5.017830225947718*1e-3*unit2;
        ap[1]= 9.706136589531091*1e-3*unit2;   bp[1] =-38.834201910782532*1e-3*unit2;
        ap[2]= -0.284540556198167*1e-3*unit2;  bp[2] =27.683608547882827*1e-3*unit2;
        ap[3]= -0.167535133865372*1e-3*unit2;  bp[3] = 58.738807959130206*1e-3*unit2;
        ap[4]=  0.853012579626518*1e-3*unit2;  bp[4] = 40.840753499612710*1e-3*unit2;
        ap[5]=  0.007048670270405*1e-3*unit2;  bp[5] = 2.191086763253139*1e-3*unit2;
        ap[6]= -0.154422262531675*1e-3*unit2;  bp[6] =46.764198176644705*1e-3*unit2;
        ap[7]=  0.980765138455104*1e-3*unit2;  bp[7] =91.799098771147001*1e-3*unit2;
        ap[8]=  0.004198973656033*1e-3*unit2;  bp[8] = 2.990746396446802*1e-3*unit2;

        gp_lw[0]= 0.0000000000000199*unit2;  gp_p[0]= .2207556531477929*unit2;  gp_m[0]=.1076924853799937*unit2;// vaterite
        gp_lw[1]= 0.0003658317849058*unit2;  gp_p[1]= .1780614458701922*unit2;  gp_m[1]=.1685352653095696*unit2;// aragonite
        gp_lw[2]= 0.0715422444509332*unit2;  gp_p[2]= 0.0894071263931514*unit2;  gp_m[2]=0.0058486656334123*unit2;// Calcite

        gp_lw[0]= 0.0018980587343381*unit2;  gp_p[0]= .2100242425379460*unit2;  gp_m[0]= 0.0961600664048892*unit2;
        gp_lw[1]= 0.0055395803384098*unit2;  gp_p[1]= .1679578085945748*unit2;  gp_m[1]= .1436824924898262*unit2;
        gp_lw[2]= 0.0691907017890783*unit2;  gp_p[2]= .1009815838151020*unit2;  gp_m[2]= 0.0095428777571714*unit2;




    }

    rhos= vector<double>(npoly,1.0) ;  
 
    if (npoly>=1)
       rhos[0]=1.0/40.0*unit2; // gmol/microm^3    This is the approx. Density of ACC //From Sean // Fatz rho=1.9 g / cm^3
    if (npoly>=2)
       rhos[1]=2.66/MWCC*unit2;  //gmol/microm^3
    if (npoly>=3)
        rhos[2]=2.93/MWCC*unit2; //gmol/microm^3
    if (npoly>=4)
        rhos[3]=2.71/MWCC*unit2; //gmol/microm^3
    // Rhos from webmineral.com  and are fairly consistent with other lieteraute sources

    ro= vector<double>(npoly,80.0) ;  

    if (npoly>=1)
        ro [0]=80e-9/unit;    //micrometers
    if (npoly>=2)
        ro [1]=80e-9/unit;
    if (npoly>=3)
        ro [2]=80e-9/unit;
    if (npoly>=4)
        ro [3]=80e-9/unit;
     Lhet=false;

     //-------------------------------HETEROGENOEUS NUCLEATION PARAMETERS--------------------------///
    if( Lhet) {
        gammaAC = vector<vector< double > > (npoly,vector<double>(npoly,1.0));
        // rates Acc               vaterite     aragonite        calcite
        gammaAC[0][0]=0.0          ;       gammaAC[0][1]=1;            gammaAC[0][2]=1;    gammaAC[0][3]=1.0;  // Acc
        gammaAC[1][0]=gamma_o[1]*.000005  ;     gammaAC[1][1]=0.0;  gammaAC[1][2]=0;    gammaAC[1][3]=0;  // vaterite
        gammaAC[2][0]=1.0+gamma_o[2]*1.95   ;gammaAC[2][1]=1;            gammaAC[2][2]=1;    gammaAC[2][3]=1;  // Aragonite
        gammaAC[3][0]=gamma_o[3]*.1  ;     gammaAC[3][1]=gamma_o[3]*.01;  gammaAC[3][2]=1;    gammaAC[3][3]=gamma_o[3]*.9;  // Calcite
        // VECTOR   [  Current POlymorph     ][ Other Polymorph     ] 

        for (int kk=0; kk<npoly; kk++) // polymorph crystalizing
        for (int j=0; j<npoly; j++) // polymorph crystalizing
        mp[j][kk]= (gamma_o[j]-gammaAC[j][kk])/gamma_o[kk];
    }

//-------------------------Other Properties-------------------//

    A= vector<double>(npoly,0.0) ;  
    Ceq= vector<double>(npoly,1.0) ;  
    Temp= 300.0;  // intial value// should never be used.
    Jm   = 1e30;    //  particles / meter^3 /  sec   // Should be 1e36? @ 298  

    
    // Do  = 6.62867E-10/unit2;// microm^2/s  // CaCO3o  // Diffusivity of CaCO_3 in water at 298
    Diff0  = 6.628e-10/unit2;  // m^2/s

//------------------------------------------------------------//
   

    r1   = 10.0*1e-9/unit;    //nanometers should be reset


    Source   =  vector<double> (nmom*npoly ,0.0);
    mu_0   =  vector<double> (nmom ,0.0);
    Ksp=      vector<double> (npoly,0.0);
    double factor;
    factor =1.0e-5;   // Current fit was ran using this value, supercomputer used following value
   // factor =1.0;   // This value, initializes the code with a 3rd moment that is too large
   // factor =1.0e-20;

    // where factor = 1 this initializes a flat uniform particle distribution ranging from zero to .01 centimeters with .1*factor particles per meter cubed  
    if (npoly>0){                                                    
        for(int i=0; i<odtlp->ngrd; i++)
            for (int j=0; j<npoly; j++){
                odtlp->mom[0+j*nmom][i]=1e-3/odtlp->rho[i]                 * factor ; // hardcoded initial conditions, Derek, momx, momobj  // Ideal for QMOM
                odtlp->mom[1+j*nmom][i]=5e-6/odtlp->rho[i]*1e2/unit        * factor ;
                odtlp->mom[2+j*nmom][i]=3.33333e-8/odtlp->rho[i]*1e4/unit2 * factor ;
                odtlp->mom[3+j*nmom][i]=2.5e-10/odtlp->rho[i]* 1e6/unit3   * factor ;

  //              odtlp->mom[0+j*nmom][i]=0 ; // hardcoded initial conditions, Derek, momx, momobj  // Ideal for QMOM
  //              odtlp->mom[1+j*nmom][i]=0 ;
  //            odtlp->mom[2+j*nmom][i]=0 ;
  //              odtlp->mom[3+j*nmom][i]=0 ;
//--------------------------------For heteronucleation----------------------------//
           if (Lhet){  // changed initial conditions for heterogenous nucleation
                odtlp->mom[0+j*nmom][i]=1e-1/odtlp->rho[i]                 * factor ; // hardcoded initial conditions, Derek, momx, momobj  // Ideal for QMOM
                odtlp->mom[1+j*nmom][i]=5e-12/odtlp->rho[i]*1e2/unit        * factor ;  // these cause a segmentation fault in the adaptive ode function
                odtlp->mom[2+j*nmom][i]=3.33333e-8/odtlp->rho[i]*1e4/unit2 * factor ;
                odtlp->mom[3+j*nmom][i]=2.5e-10/odtlp->rho[i]* 1e6/unit3   * factor ;
            }
            if (nmom>4){                                                    
                odtlp->mom[4+j*nmom][i]=2e-12/odtlp->rho[i]*1e8/(unit3*unit)          *  factor ;             
                odtlp->mom[5+j*nmom][i]=1.66667e-14/odtlp->rho[i]*1e10/(unit3*unit2)  *  factor ;
            }
            if (nmom>6){
                odtlp->mom[6+j*nmom][i]=1.4285714e-16/odtlp->rho[i]*1e12/(unit3*unit3) *  factor ;             
                odtlp->mom[7+j*nmom][i]=1.25e-18/odtlp->rho[i]*1e14/(unit*unit3*unit3) *  factor ;

            }





            }

        for (int j=0; j<nmom; j++)
            mu_0[j]= odtlp->mom[j][0]*odtlp->rho[0]; // comparison vector   only element with index=3 is used


        if (Ldqmom){
            vector<double> mu(nmom,0.0);
            vector<double> wts(np2,0.0);
            vector<double> absc(np2,0.0);
            M3i= mu_0[3];//// index is zero in rho[0] because itital moments are uniform this would need to be put in the loop and made into array if  a non uniform profile with respect to x is used
            


            for(int i=0; i<odtlp->ngrd; i++){
                for (int j=0; j<npoly; j++){
        
                    for (int k=0; k<nmom; k++)
                        mu[k]=odtlp->mom[k+j*nmom][i]*odtlp->rho[i];

            //  if (mu [0+j*nmom] <= 0){
            //  wts[0]=0;
             // wts[1]=0;
              //absc[0]=0;
             // absc[1]=0;
             //   }
             // else
                        x.pdAlg(nmom, np2, mu, wts, absc);




                    for (int k=0; k<np2; k++){
                        odtlp->mom[k+j*nmom][i]=wts[k]/odtlp->rho[i]; // hardcoded initial conditions, Derek, momx, momobj
                        odtlp->mom[k+np2+j*nmom][i]=absc[k]/odtlp->rho[i]; // hardcoded initial conditions, Derek, momx, momobj


                        odtlp->mom[0+j*nmom][i]=1e6/odtlp->rho[i]; // hardcoded initial conditions, Derek, momx, momobj
                        odtlp->mom[1+j*nmom][i]=1e6/odtlp->rho[i]; // hardcoded initial conditions, Derek, momx, momobj
                        odtlp->mom[0+np2+j*nmom][i]=1e-9/odtlp->rho[i]/unit; // hardcoded initial conditions, Derek, momx, momobj
                        odtlp->mom[1+np2+j*nmom][i]=2e-9/odtlp->rho[i]/unit;  // hardcoded initial conditions, Derek, momx, momobj
                        if (nmom>4){
                        odtlp->mom[2+j*nmom][i]=1e6/odtlp->rho[i]; // hardcoded initial conditions, Derek, momx, momobj
                        odtlp->mom[2+np2+j*nmom][i]=5e-9/odtlp->rho[i]/unit;  // hardcoded initial conditions, Derek, momx, momobj
                        }
                         

                    }  // End loop for weights and Absicca
                } // End loop for moments
            }  // end loop for grid points
        }  // End DQMOM if statement  
    } // End if statement for Polymorphs  // without this statement function segfaults when moments are off because all moments are of zero size.
} // End constructor





void mom_Aq::EvalPropAtTemp(double H){
   
// In this function the K_sp*K_s / gamma_CaCO3 , however gamma_CaCO3 is assumed to be equal to one (molar gamma)
    Temp = ETApx->getTgivenH(H);

    double K_Reaction;
  //  K_Reaction= pow(10.0,-1228.732 - 0.299444  *Temp + 35512.75/Temp + 485.818 *log10(Temp)); // 
    K_Reaction= pow(10.0,-1228.732 - 0.299444  *Temp + 35512.75/Temp + 485.818 *log10(Temp));
  //  K_Reaction=1675.29 ;

    //double rho_o, rho_T, T0, Ratio_rho ; //  !!!!!     currently unused variable
    //--------------------------------Takes into expansion in the Crystal structure, not being used---------//
    /*
    T0=273.15+25.0;  // Reference Temperature!   
    rho_T = 1.8211e-6*Temp*Temp*Temp-1.9124e-3*Temp*Temp+6.4038e-1*Temp-1.3851e1 ;  // Dipr
    rho_o = 1.8211e-6*T0*T0*T0-1.9124e-3*T0*T0+6.4038e-1*T0-1.3851e1 ; // Dipr
    Ratio_rho= rho_T/rho_o;
    cout << R << "  RRRRRRRRRRRRRRR \n";
    for (int j=0; j<npoly; j++)
          A[j] = (16.0*3.14159265/3.0)*pow((gamma[j]/Temp),3.0)  * pow((1.0/Ratio_rho/rhos[j]),2.0)/kb;   
   Jm= 1e30 *pow( (Temp/298.0),1.5 ); // currently not used.  
   */
   //-------------------------------------------------------------------------------------------------------//
    if (npoly>=1)
        Ksp[0]=pow(10.0,-(6.1987 + 0.00053369 *(Temp-273.15) + 0.0001096 *(Temp-273.15)*(Temp-273.15)));   // Expressions from Nernk  // Gebaeur pH of 10.0
       // Ksp[0]=pow(10.0,-(6.1987 + 0.00053369 *(Temp-273.15) + 0.0001096 *(Temp-273.15)*(Temp-273.15)))*.025;   // Expressions from Nernk  // Gebaeur pH of 10.0
        //Ksp[0]=pow(10.0,-(6.1987 + 0.00053369 *(Temp-273.15) + 0.0001096 *(Temp-273.15)*(Temp-273.15)))*.1;   // Expressions from Nernk  // Gebaeur pH of  9.0
    if (npoly>=2)
        Ksp[1]=pow(10.0,-172.1295 - 0.077993* Temp + 3074.688 / Temp + 71.595* log10(Temp)); // gmol/L
    if (npoly>=3)
        Ksp[2]=pow(10.0,-171.9773 - 0.077993*Temp + 2903.293/Temp + 71.595* log10(Temp)); // gmol/L  // This is different than Nernk, however nernk is believed to have a typo in this eq.
    if (npoly>=4)
        Ksp[3]=pow(10.0,-171.9065 - 0.077993 *Temp + 2839.319 / Temp + 71.595* log10(Temp)); // gmol/L

    for (int i=0 ; i<npoly; i++)
        Ceq[i]=K_Reaction*Ksp[i]/CC_activity;   //Ceq=K_reaction * K_sp / CaCO3_activity_coefficient (molar)


    //   D  = Do/(298.0)*Temp;// microm^2/s  // CaCO3o  // Diffusivity of CaCO_3 in water at 298// BSL 17.4-8
    Diff  = Diff0*exp(-24380/8.31447*(1.0/Temp-1.0/298.15));   // From karlj



    for ( int j=0; j<npoly ; j++)
        gamma[j]=gamma_o[j];//*pow(R,2.0/3.0);
// ----------------------- using Van Oss theory--------------------------------------------// 

    if (Lvanoss){


        double gw_lw, gw_p, gw_m;
        //-------------Linear gw_p and gw_m--------------------------//
        gw_lw=( 22.8  + (22.8-21.0) / -38.0*(Temp-273.15))*1.0e-3*unit2;
        gw_p =( 25.5  + (25.5-32.4) / -18.0*(Temp-273.15-20.0))*1.0e-3*unit2;
        gw_m =( 25.5  + (25.5-18.5) / -18.0*(Temp-273.15-20.0))*1.0e-3*unit2;
        //-----------------------------------------------------------//A

        vector<double> al(3,0.0);
        vector<double> bl(3,0.0);
        al[0]=0.0004060217983651*unit2; bl[0]=   .1179976839237058*unit2;
        al[1]=0.0001815000000000*unit2; bl[1]=   .1522300000000001*unit2;
        al[2]=0.0003507029972752*unit2; bl[2]=   .1421328065395096*unit2;

        for (int j=0; j<npoly-1; j++){
            gamma[j+1]= pow(sqrt(gp_lw[j])+sqrt(gw_lw),2.0) + 2.0*(sqrt(gp_p[j]*gp_m[j]) +  sqrt(gw_p*gw_m) -  sqrt(gw_p*gp_m[j]) - sqrt(gp_p[j]*gw_m)) ;
            gamma[j+1]= al[j]*(Temp-273.15)+bl[j];
        }


    }


}





double mom_Aq::get_moment_source(vector<double> mu2, double Caq, double H,double CC_activity_dummy){
// This function operates on a per liter basis!
// This function can handle mol/L imput or kmol/m^3 input.
// It must return units in kg/s/m^3
    CC_activity=CC_activity_dummy;
   // cout << Ceq[0] << " \n ";
//--------------------------flag to ensure that if this function crashes its not its false----------------//
    if (isnan(Caq)){

        for (int i=0; i<nmom; i++)
            cout << mu2[i] << " ";
            cout << Caq << " NANANANNAN  FED to MOMENTS";
            int delme;
            cin >> delme;
            }
//--------------------------------------------------------------------------------------------------------//

    EvalPropAtTemp(H);
   // get_interfacial_energy(Caq);


    vector<double> mu(nmom,0.0);
   


    double dCaq=0.0;
    for (int i=0; i<npoly*nmom; i++)// set source terms to zero, they are global
        Source[i]=0.0;

         //    Caq=6.314556573856522e-03;  // DEERK DEBUG DELETE DESTROY
             
    for (int k=0; k<npoly; k++){           // Large Polymorph loop
        SS=Caq/Ceq[k];
//--------------_Quick fix to problems matrix inverting problem / dissolution--------//  DEREK DEBUG  // Also faster for 1-D cases
if (SS < 1.0) 
    continue;
//--------------------------------------------//

//-----------------For dissolution cases---------------//  // FOR SOME REASON mu2[3+k*nmom] is larger than it should be!  
    /*    if ( SS < 1.0 )
            if (Ldqmom ){
                double  DQsum=0.0;
                for (int i=0; i<np2;  i++)
                    DQsum+= mu2[i+k*nmom]*mu2[i+np2+k*nmom]*mu2[i+np2+k*nmom]*mu2[i+np2+k*nmom]; //sum( w * a ^3)
                if   (DQsum < M3i)
                    continue;
            }
            else{
                if (mu2[3+k*nmom] < mu_0[3]*2.0  ) 
                    continue;
                if (mu2[2+k*nmom] < mu_0[2]*.0001   )   // This is due to the negative moments generated by QMOM and so this value goes negative before Mu2[3] does
                    continue;
                if (mu2[2+k*nmom] < mu_0[1]*.0000001   )   // This is due to the negative moments generated by QMOM and so this value goes negative before Mu2[2] does
                    continue;
            }
            */
//-----------------------------------------------------//


            for (int j=0; j<nmom; j++)
                mu[j]=mu2[j+k*nmom];

            if (SS<1.0)

                r1=1;            // This is because r1 goes to infinity for small values of S, predicting an infinite growth term, NAN*0 = NAN increases stability, if step size is too large  (may not be needed)
            else
                r1=ro[k]/log(SS);// Negative r1 for SS < 1

//----------------Nucleation delta function-----------------------///
/*        if (SS > 1.0){
            B=Jm *exp(-A[k]/pow(log(SS),2.)); // Birth is ALWAYS positive, Turn it off for SS less than one =/
            for(int i=0; i<nmom; i++)
                S[i+k*nmom]=B*pow(r1, (double)  i);
        }
        */
//-------------------------  Ben Nucleation -------------------------/// 
            double DeltaCritGibbs, Critical_i, Zelda, Kf, N1 ;
            double pi=3.141592653589793;
            double nu;
            double r_c;
            double N_A=6.022141793e23; //Avagadro's number
            nu = 1./(rhos[k] *N_A);
            if (SS > 1.0){
            B=0;
            DeltaCritGibbs = (16*pi/3) * nu*nu * pow(gamma[k],3.0) / pow(kb * Temp * log(SS),2.0) 
            + kb * Temp * log(SS) - gamma[k] *pow( (36.0 * pi * nu*nu),(1.0/3.0)); 
            Critical_i = (32.0 * pi / 3.0) * nu*nu * pow(gamma[k],3.0) / pow((kb * Temp * log(SS)),3.0) ;

            Zelda = pow((DeltaCritGibbs / (3.0 * pi * kb * Temp * Critical_i*Critical_i)),(.5));

            Kf = Diff * pow((6.0 * pi*pi * nu),(1.0/3.0)) * pow(Critical_i,(2.0/3.0));   // Kf from equation 28 in the Kashchiev paper

 //           Kf =pow(243.0*pi*nu*nu/16.0,1.0/3.0)*Diff*pow(Critical_i,1.0/3.0)*pow(N1,1.0/3.0) ;     // Kf from equation 30 in the kashcheiv paper

            N1 = N_A * Ceq[k] * SS;   // Number / m^3
            B = Zelda * Kf * N1*N1 * exp(-DeltaCritGibbs / (kb * Temp))*unit3*1e6;  // 1e6 because Ceq is mol / L (1e3)^2

            r_c=2.0*nu*gamma[k]/kb/Temp/log(SS);
      //      B  =3.0/2.0*pow(gamma[k]/kb/Temp,.5)*nu*Diff*pow( Caq*N_A*1000*unit3  ,7.0/3.0 )*exp(-(16*pi/3) * nu*nu * pow(gamma[k],3.0) / pow(kb * Temp ,3.0)/log(SS)/log(SS));  // equation 30 kashchiev
        //     cout << B << "   \n";

           // B=Jm*exp(-A[k]/pow(log(SS),2.)); // Birth is ALWAYS positive, Turn it off for SS less than one =/ OLD nucleation
      //      cout <<setprecision(16)<< " dGc =" <<DeltaCritGibbs << "   \n";
         //   cout <<setprecision(6)<< " B =" <<B << "   \n";
       //     cout <<setprecision(16)<< " ZELDA =" <<Zelda << "   \n";
      //      cout <<setprecision(16)<< " N1 =" <<N1 << "   \n";
      //      cout <<setprecision(16)<< " Kf =" <<Kf << "   \n";
      //      cout <<setprecision(12)<< " nu =" <<nu << "   \n";
      //      cout <<setprecision(8)<< " gamma =" <<gamma[k] << "   \n";
      //      cout <<setprecision(8)<< " S =" <<SS << "   \n";
      //      cout <<setprecision(8)<< " Temp =" << Temp<< "   \n";
           // cout << " \n";

            for(int i=0; i<nmom; i++)
                Source[i+k*nmom]=B*pow(r_c, (double) i);
            }
//--------------HeteroGeneous nucleation--------------------------//
            if(Lhet){
            double w,fpp,f,X;//,Jp; //  Variable Jp is unused




                for (int i=0; i<npoly; i++) // polymorph being nucleated upon
                    for (int j=0; j<npoly; j++){ // polymorph crystalizing
                        if (mp[i][j]<-1.0)
                            mp[i][j]=-1.0;
                        if (mp[i][j] > 1.0)
                            mp[i][j]=1.0;}

                        for (int i=0; i<npoly; i++){
                            if (SS<1)
                            continue;
                    //      if (i==k) // Turns off nucleation reduction due to like crystals around
                    //         continue;
                        double Rs=mu2[1+i*nmom]/mu2[0+i*nmom];
                        X=Rs/r1;  //r1 is critical size, doesn't match seans 
              
                        X=X*.001;

                        w=pow((1.0+X*X-2.0*X*mp[k][i]),.5);

                        fpp=(1.0+(1.0-X*mp[k][i])/w)*.5;

                        f  = 1.0 +  pow((1.0-mp[k][i]*X)/w,3.0)
                        + X*X*X*(2.0-3.0*(X-mp[k][i])/w  + pow((X-mp[k][i])/w,3.0))
                        + 3.0*mp[k][i]*X*X*((X-mp[k][i])/w-1.0);
                        f*=.5;
//-------------------------Keep until you have lots of confidence in hetero--------//
                        if (f>1.1){
                            cout << X << "  " << mp[k][i] << "  " << f <<  "   "   << fpp << " \n";
                            f=1;
                        }   
                        if ( f<0){
                            f=abs(f);
                        }
//---------------------------------------------------------------------------------//


                        double DeltaCritGibbsH= DeltaCritGibbs;
                        DeltaCritGibbsH *=f;
                        Zelda = pow((DeltaCritGibbsH / (3 * pi * kb * Temp * Critical_i*Critical_i)),(.5));
                        double a=2.0*2.378e-10/unit;    // may be as high as 6.2e-10?
                        double B2;
                        B = Zelda * Kf * N1*N1 * exp(-DeltaCritGibbsH / (kb * Temp))*unit3*1e3*1e3;
                        B2 = 4.0*pi*a*Rs*Rs*mu2[0+i*nmom]*fpp*pow(f,.5)*B*unit3;



                        for(int ii=0; ii<nmom; ii++)
                            Source[ii+k*nmom]+=(B2)*pow(r1, (double)  ii);

                    }

            }

//----------------------------------------------------------------//

//---------------FAT delta function approach for Nucleation------------------------//
/*
if (SS>1.0){
      B=Jm*pow(2.71828,-A[k]/pow(log(SS),2)); // Birth is ALWAYS positive, Turn it off for SS less than one =/
    double dnuc;
    dnuc =.00000001*r1;
//    (for this value of dnuc, the function behaves the same!!! see if numerical stability increases!!!!!)


      S[0+k*nmom]=B ;     // Same as delta
      S[1+k*nmom]=B*r1 ; // Same as delta 
      S[2+k*nmom]=B*(r1*r1 +(1.0/3.0)*dnuc*dnuc) ;
      S[3+k*nmom]=B*(r1*r1*r1 + dnuc*dnuc*r1) ;
      if (nmom>4){
      S[4+k*nmom]=B*(r1*r1* r1*r1 +(2.0)*dnuc*dnuc*r1*r1 +0.2*dnuc*dnuc*dnuc*dnuc) ;
      S[5+k*nmom]=B*(r1*r1*r1*r1*r1 +10.0/3.0 * dnuc*dnuc * r1 *r1 *r1 + dnuc*dnuc*dnuc*dnuc*r1) ;
      }   
}
*/
//---------------------------------------------------------------//
        if(Ldqmom){

//---------------------------------------------DQMOM--------------------------------------------//
        vector<int> Dead_Absc(0);
        int dummy_n=0;
//----------------- Dissolution trigger -----------------------//
        if (SS < 1.0)
            for (int i=0; i<np2; i++)
                if (mu[i+np2]<1e-20  || mu[i]<1e-10 ){
                    dummy_n=1+Dead_Absc.size();
                    Dead_Absc.resize(dummy_n);
                    Dead_Absc[dummy_n-1]=i;
           //       if (mu[i+np2] < 1e-20)  cout << mu[i+np2]<< "  abscissa          MU \n";
           //       if (mu[i] < 1e-10)  cout << mu[i]<< "   weight          MU \n";
        }
//-------------------------------------------------------------//
        int nmom_tran=nmom-dummy_n*2;
        //int np_tran  = np2-dummy_n; //  !!!!!      currently unused variable
       //     cout << nmom_tran << " \n";

        vector<vector<double> >  A(nmom_tran,vector<double>( nmom_tran,0.0));
        vector<double> b(nmom_tran,0.0);
        vector<double> x(nmom_tran,0.0);

///*     
//if (nmom_tran==2)
//    for (int nm=0; nm<dummy_n; nm++)
//        cout << " Index " << Dead_Absc[nm] << "    ";

         //-------------- For Dissolution--------------//
            vector<double> dissolution_coefficient(np2,1.0);
           if (SS< 1.0) {
               double  r_d = .01;          // micrometers
               double delta       =   .05; // micrometers;
               for (int i=0; i<np2; i++)
              dissolution_coefficient[i] =  .5*(1.0 + tanh(2*(mu[np2+i] - r_d) / delta)  ) ;
            //  dissolution_coefficient[i] =  1.0 ;
 }

         vector<double> Growth(np2,0.0);
         vector<double> Dissolution(np2,0.0);
        //---------------------------------------------//
            for (int i=0; i<np2; i++){
            Growth[i]=(SS-Sb(mu[i + np2] ,k ) ) / mu[i+np2]  * Diff/rhos[k]*Ceq[k]*1e-15;
       //   Dissolution[i]=2*(SS-Sb(mu[row+ np2] ,k ) ) / mu[row+np2]/mu[row+np2]  * Diff/rhos[k]*Ceq[k]*1e-15;
            Dissolution[i]=2.0*Growth[i]/mu[i+np2] ;
            }

        int dummy2_n=0;
        for (int row=0; row<nmom_tran; row++){
            for (int column=0; column<np2; column++){
                if (dummy_n > dummy2_n)
                   if (Dead_Absc[dummy2_n]==column){
                      dummy2_n+=1;
              //        cout << " KILLING THE ABCISSA!!! \n";
                      continue;
                      }
          //        if (nmom_tran==2){
          //        cout << row<< " " << column << " \n";
          //        cout <<" ASSIGN these INDEXes " <<row <<"  "<<(column-dummy2_n)*2<< " " << 1+(column-dummy2_n)*2 << " \n";

            //      int delme;
             //     cin >> delme;
             //     }
                                  // abscissa^(a-1)                 
                A[row][(column-dummy2_n)*2] = pow(mu[column+np2],(double) row); // weights column
                                  // weight               // abscissa^(a-1)                 //a
                A[row][1+(column-dummy2_n)*2]=mu[column ]*pow(mu[column+np2], (double) row - 1.0)*(double) row; // Absicca Column


            Source[row+k*nmom]+= (double) row*pow(mu[column+np2],(double) row-1.0)*mu[column]*Growth[column]*dissolution_coefficient[column] + pow(mu[column+np2],(double) row)*mu[column]*Dissolution[column]*(1.0-dissolution_coefficient[column]);
            }      
            dummy2_n=0;
            b[row]=Source[row+k*nmom];   // Birth
        }

        vector<vector<double> >  Ax(nmom_tran,vector<double>( nmom_tran,0.0));
        Ax=A;

        ETApx->QR(A,b,x);



//--------------------Flag to Check to see if matrix is singular -------------------------//
if (nmom_tran>0)
        if (isnan(x[0])){
//---------------------------View A matrix, and b matrix-----------------------------------//
                cout <<" \n \n "<< setprecision(16) << SS<< " \n";
                cout << Dead_Absc.size() << " \n";
                cout << dummy2_n <<  "   " << nmom_tran << "       " << dummy_n<<"    " << nmom << "  \n";
                cout <<  mu[0] << "   " << mu[1] << "   " <<mu[2] << "   " << mu[3] <<" \n";
                cout << pow(mu[0+np2],(double) 0) << "  " << pow(mu[0+np2],(double) 1) << " \n";
                cout << pow(mu[1+np2],(double) 0) << "  " << pow(mu[1+np2],(double) 1) << " \n";
                cout << "\n";
cout<<mu[1]*pow(mu[1+np2], (double) 0 - 1.0)*(double) 0 << "   \n"<<mu[1 ]*pow(mu[1+np2], (double) 1 - 1.0)*(double) 1  << " \n";
cout <<(mu[0+np2]<1e-20  || mu[0]<1e-10 )<< "  " <<( mu[1+np2]<1e-20  || mu[1]<1e-10 )<< " \n";
cout << " A =";
                  if (nmom_tran==2){
                cout << " \n";
        for (int row=0; row<nmom_tran; row++){
            for (int column=0; column<nmom_tran; column++)
                cout <<    Ax[row][column] << "      ";
                cout << "                                         \n";
}
        cout << " \n";
        for (int row=0; row<nmom_tran; row++)
                cout <<  b[row] << " \n";
         //       int delme7;
         //       cin >> delme7;
        cout << " \n";
        for (int row=0; row<nmom_tran; row++)
                cout <<  x[row] << " \n";
        cout << " \n";
          }

        for (int row=0; row<nmom_tran; row++)
                cout <<  b[row] << " \n";
                cout << " \n";

        for (int row=0; row<nmom_tran; row++)
                cout << setprecision(16) <<  x[row] << " \n";
                int delme;
                cin >> delme;
//------------------------------------------------------------------------------------------//
            cout << "\n********* Singular matrix in MOM class NANANANANANANAN****************\n";
            cout << "   MOMENTS ARE:  "<<  mu[0] << "   " << mu[1] << "   " << mu[2] << "   " << mu[3] << " \n ";
            cout << "    nmom_tran =  "<<nmom_tran <<"  Dead_Absc.size = " << Dead_Absc.size()<<"    \n";
            int delme2;
            cin >> delme2;

        }
            
//----------------------------------------------------------------------------------------//

//---------------------------------source terms for weights and abscissa generated in x-------------------//
   //     double xx;
  //      xx=0.0;
 //       for(int j=0; j<nmom/2; j++)                          // 3 refers to the third moments
//            xx+=mu[j]*mu[j+np2]*(SS-Sb(mu[j+np2],k) )*Diff/rhos[k]*Ceq[k]*3.0*1e-15 +Source[3+k*nmom];  

            dCaq +=-Source[3+k*nmom]*4.0/3.0*3.14159265*rhos[k]*MWCC*.001;   // Units of kg/s/m^3  
        for (int i=0; i<nmom; i++)
            Source[i+k*nmom]=0;




        dummy2_n=0;
        for (int i=0; i<np2; i++)
                if (dummy_n > dummy2_n){
                   if (Dead_Absc[dummy2_n]==i){
                        dummy2_n+=1;
                        Source[i+k*nmom]=0.0;
                        Source[i+np2+k*nmom]=0.0;
                        }
                    }
                    else{
                        Source[i+k*nmom]     =x[(i-dummy2_n)*2];  // Weights_dot
                        Source[i+np2+k*nmom] =x[(i-dummy2_n)*np2+1];   //Abcissa_dot
                    }
     


            //      if (nmom_tran==2){
            //          cout << mu[0] << "   " << mu[1]<<" " <<Source[0] << "  " << Source[1]<<"   "<< Source[2]<<"  " << Source[3] << " \n";
            //          cout << dummy_n <<"  "<< Dead_Absc[dummy_n]<<"   " << " \n";
    //    cout << " \n";
      //            int delme;
       //           cin >> delme;
        //          }
//----------------- Dissolution trigger -----------------------//
//        for (int i=0; i<np2; i++)
//            if (mu[i+np2]<1e-15 && SS <1.0)
//                Source[i+np2+k*nmom]=0.0;
//-------------------------------------------------------------//

//cout << Source[0 +k*nmom] << "   " << Source[1+k*nmom] << "    " << Source[2+k*nmom] << "    " << Source[3+k*nmom]<< "  \n";
        }
        else{
///////////////////////////  JQMOM //////////////////////////
            vector<double> wts(np2,0.0);
            vector<double> absc(np2,1.0);


            if (mu [0] <= 1e-30);  // Check to see if moments are zero.  // This should be zero, but when it is the PD algorithm is called too early and this crashes.
           else if (abs( (mu [1] - mu [3]/ mu[2]*mu[0])/(mu [1])) <= 1e-12){  // Check to see if moments can be described by one environment
              wts[0]=mu[0];
              absc[0]=mu[1] / mu[0];
              cout << " Abbreviated PD algorithm called\n ";
              }
              else
                x.pdAlg(nmom, np2, mu, wts, absc);  
            
/////////////////////////////////////////////////////////////////////////////////////
            double xx; // I did this so that I didnt' need a modular nucleation and growth function, although may still be useful eventually, and to save on compuatation of not needing to multiply those 5 terms each time in loop;
            xx=0.0;    // Originally I did growth first, then I didn't need xx;
//-------------Diffusion limited growth ------------/
            for(int i=1; i<nmom; i++){
                for(int j=0; j<np2; j++)
                    xx+=wts[j]*pow(absc[j], (double) i-2.0 )*(SS-Sb(absc[j],k));   //  //double and int interact` // in comparison with the simpler diffusive growth

                Source[i+k*nmom]+=xx*Diff/rhos[k]*Ceq[k]*(double)i*1e-15;
                xx=0.0;
            }
//--------------------------------------------------//
//------------------Aggregation---------------------//
/*
double Beta=0 ;
double visc=1e-3 ;  //   kg / m / s
double efficiency=3e-5;
efficiency=1.0;
double Vel_Grad=10;   //   m / s / m   
            for(int z=0; z<nmom; z++)
                for(int j=0; j<np2; j++)
                    for(int i=0; i<np2; i++){
                           Beta=pow(absc[i] +absc[j],2.0 ) / absc[i]/absc[j] *kb*Temp*2.0/3.0/visc; // brownian collision  m^3/s
                           Beta=pow(absc[i] +absc[j],3.0 ) /  6.0 * Vel_Grad;
                           Source[z+k*nmom]+= wts[j]*wts[i]*Beta*(pow(absc[i]+absc[j],(double) z) -2.0*pow(absc[i], (double) z)); // ( micrometers^k / m^3 / s )
                           if (absc[i] < 0)
                               cout << " Abcissa is less than zero   " <<k <<"   " << absc[i] << "  " << i <<"   " << "\n";
                    }
                    */
//--------------------------------------------------//


            dCaq   +=-4.0/3.0*3.14159265*Source[3+k*nmom]*rhos[k]*MWCC/1000.;   // Units of kg/s/m^3  or g/s/L
//            cout << -4.0/3.0*3.14159265*S[3+k*nmom]*rhos[k]*MWCC/1000. << "   " << k  << "  " << Caq<< " \n " ;
/*
if (k==3){
            cout << -4.0/3.0*3.14159265*S[3+k*nmom]*rhos[k]*MWCC/1000. << "   " << k  << "  " << Caq<< " \n " ;
if (dCaq > 0 ){

            int delme7; 
 cin >> delme7  ;
}
}

*/

//-------------------checks if this function is returning nans that this function doesn't return nans---------------------------//

            if ( isnan(Source[3+k*nmom])  ){
                cout << "Error in  Polymorph  " << k << "  possibly due to PD algorithm\n"; 
                cout << SS << " \n";
                cout << "\n Derivative of second moment =  "<<Source[2+k*nmom] << "\n";
                cout << "\n Derivative of Third moment =  "<<Source[ 3+k*nmom] << "\n";

                for (int i=0; i<nmom; i++)
                    cout << "M_" <<i<< " = "<< mu[i] << " \n";

            }
//------------------------------------------------------------------------------------------------------------------------------//

        }  // End QMOM
//-----------------------------------------End  QMOM----------------------------------------------------------------------------//

    }   // End of Polymorph loop

//cout << " \n";
//for (int i=0; i<16 ; i++)
//    cout << Source[i] << "  ";
//    cout << dCaq << "   \n";
//    int B22;
//    cin >> B22 ;
    return dCaq;
}




double mom_Aq::Sb(double  radius, int k){    // Ostwald ripening function

  if (~Lostwald)     
      return 1.0;

   return exp(2.0*gamma[k]/(R*Temp*radius*rhos[k])) ;   // ostwald ripening

}








//-------------------old growth originated in get_moment_source  DQMOm--------------------//
//      double mkp1=0;
//     for (int i=0; i<np2; i++)
//      mkp1+=wts[i]*pow(absc[i],nmom);  // nmom = k+1
//  for(int i=0; i<nmom-1; i++)
//  S[i]=(K*Ceq[k]*(SS-1)*mu[i+1]*(double)i); // double and int mingle here
//  S[nmom-1]=(K*Ceq[k]*(SS-1)*mkp1*(nmom-1));
//--------------------------------------------------//








///////////////////////////////////////////////////////////////////////////////

/** 
 *
 * 
 */
// This function is not being used.
bool mom_Aq::check(vector<double>  mu){
vector<double> mu2(nmom);
vector<double> wts(np2);
vector<double> absc(np2);


for (int j=0; j<npoly; j++){

for (int i=0; i<nmom; i++)
mu2[i] = mu[i+j*nmom]*1000;

x.pdAlg(nmom, np2, mu2, wts, absc);

/*   for(int j=0; j<nmom; j++)
cout << mu[j] << "  ";
cout << "\n  ";
*/
//////// CHECK FOR PD ALGORITHM Useful if distribution is Illconditioned/////////////

   for(int jj=0; jj<np2; jj++)
       if (absc[jj]<0)
           return false;
}
       return true;
}




 //I would like to move the adaptive ODE class to the MOM class, however I can't because it needs varaibles defiend in the diffuser i could pass diffuser variables to this code =/
void mom_Aq::AdaptiveODE(double dtStep, int &RD, vector<double> &dxML ,vector<vector<double> > &flxEta){

//----------------------------------------------------------Warning flags ------------------------------------//
//    // Warning CALLED during gebaur!!!!   // this is a problem because the gebaur source terms are not implimented in this function
//    cout << " RECURSIVE CALLED DURING GEBAUR EXPERIMENT!!!!!";
//------------------------------------------------------------------------------------------------------------//
//cout <<" \n ***** Welcome to the Adaptive ODE solver, may heaven have mercy on your soul ******\n " ;
//int delme;
//cin >>delme;

///////////////////////////////////////////////////////////////////////////
    //bool PDError; //  !!!!!     currently unused variable
    vector<vector<double> > D_eta(odtlp->neta, vector<double> (odtlp->ngrd,0.0));
    vector<vector<double> > D_mom(odtlp->nmom, vector<double> (odtlp->ngrd,0.0));

    vector<double> DC(odtlp->ngrd,0.0);
    vector<double> momx(odtlp->nmom,0.0);  // place holder for moments and eta
    vector<vector<double> > momDot(odtlp->nmom,  vector<double> (odtlp->ngrd,0.0));
    vector<vector<double> > etaDot(odtlp->neta,vector<double> (odtlp->ngrd,0.0));
    int RC=2;  // The Recursion Fraction  2 or more
    RD+=1 ;
    for (int j=0; j<RC ;j++){

        bool Lx=true;  
        bool Lx2=true;  


        for (int i=0; i<odtlp->ngrd; i++){   //Take smaller time steps for all grid points


            for(int kk=0; kk<odtlp->nmom; kk++)   //old
                momx[kk]= odtlp->mom[kk][i]*odtlp->rho[i];  //old // 

            DC[i]=get_moment_source(momx,ETApx->CaCO3[i], odtlp->eta[odtlp->neta-1][i],ETApx->CaCO3_gamma[i]);

//----------------------------------------------------------Warning flags ------------------------------------//
        if (RD > 90000){
        for (int ii=0; ii<nmom; ii++)
        cout << i << "   " << ii <<"  "<<  momx[ii] << "  " << Source[ii]<< "  " << dtStep << "  " << dtStep*Source[ii] << "  " <<DC[i] <<" in ADAPTIVEi\n "  ;
        cout << " \n";
        }
//------------------------------------------------------------------------------------------------------------//



            for(int k=0; k<odtlp->nmom; k++){
                momDot[k][i] = Source[k]/odtlp->rho[i];  // ZERO DIFFUSION OF MOMENTS!
                D_mom[k][i]=odtlp->mom[k][i] + dtStep/(double) RC*momDot[k][i];
            }
//-------------------------Conservative condition ensuring that QMOM abscissa do not go negative, is really slow =/ -----------//
/*                            if (!Ldqmom ){
                                vector<double> momx2(odtlp->nmom);
                                for(int kk2=0; kk2<odtlp->nmom; kk2++) 
                                    momx2[kk2]  = D_mom[kk2][i]*rho ;
                            Lx= check(momx);
                        }
                        */
//----------------------------------------------------------------------------------------------------------------------//

            for(int k=0; k<odtlp->nmom; k++)
                if (D_mom[k][i] < 0.0){
                    Lx=false;
             //       cout << "New Moment " << k << " = " << D_mom[k][i]<<"   Old Mom =  " <<odtlp->mom[k][i]  <<" \n";
                }

            if (!Lx){   
             //   cout << "MOMENTS \n";
                AdaptiveODE(dtStep/(double) RC , RD, dxML, flxEta);
                i=odtlp->ngrd;
                break ;
            }
            else{
                double dd;
                dd = 1.0/(dxML[i]*odtlp->rho[i]);
                for(int k=0; k<odtlp->neta; k++) {
                    etaDot[k][i] = -dd*(flxEta[k][i+1] - flxEta[k][i]); 

                    if (k==0)                                           // These if statements are redundant since Yspel is zero for values other than 0,2,3, however watch out for enthalpy if you change this
                    // Carbon
                        etaDot[k][i]+=DC[i]/odtlp->rho[i]*ETApx->Yspel[8][k]    ;// Yspel converts mass of species to mass of element, species 8 is CaCO3
                    if (k==2)
                    //Oxygen
                        etaDot[k][i]+=DC[i]/odtlp->rho[i]*ETApx->Yspel[8][k]    ;// dd converst kg/s to 1/s?  kg/s * 1/(m * kg/m^3)
                    if (k==3)
                    //Calcium
                        etaDot[k][i]+=DC[i]/odtlp->rho[i]*ETApx->Yspel[8][k]   ;  // micro L / min *(min/s)*(L/microL) * .01 mol/ L *(1kmol/1000 mol) 
                } 

                for(int k=0; k<odtlp->neta; k++){
                        D_eta[k][i] = odtlp->eta[k][i] + dtStep/(double)RC*etaDot[k][i];
                    if (D_eta[k][i] < 0.0)
                    Lx2=false;
                }
                if (!Lx2){
                    //--------------------------------------------DELETE DEBUG------------//
//                    cout <<  "   ETA!!!!!!!!!!!!!!!! \n";
                    //--------------------------------------------------------------------//
                    AdaptiveODE(dtStep/(double) RC ,  RD, dxML, flxEta);
                    i=odtlp->ngrd;
                break; ///
                }
            }
      



        }
//---------------------------------------------Save good values----------//
////////////////////////////////////////////////////////
        if (Lx & Lx2){
            odtlp->eta = D_eta;
            odtlp->mom = D_mom;

//--------------------------------Warning flag--------------------------------------//
        for (int i=0; i<odtlp->ngrd; i++)   
            for(int k=0; k<odtlp->nmom; k++)
               if (odtlp->mom[k][i]<0)
                   cout << " NEGATIVE MOMENTS in ADAPTIVE ALGORITHM ";
//----------------------------------------------------------------------------------//

           // odtlp->etaObj.getEtaDiffusiveFluxes(odtlp->eta, odtP->rho_0, dxML, flxEta, true);
            ETApx->getEtaDiffusiveFluxes(odtlp->eta, odtlp->rho[0], dxML, flxEta, false);
            ETApx->equilibrateAll(odtlp->eta,odtlp->rho[0]);
        }   
//--------------------------------------------

    }
    return;
}








void mom_Aq::get_interfacial_energy(double Caq){  // By Mersmann

    double polymorph_activity;
// Need molar density 

for (int k = 0; k<npoly; k++){

 polymorph_activity= Ksp[k] / Ksp[3];
 
cout << " \n"; 
gamma[k]=.414*kb*Temp*pow(6.022e23*rhos[k],(2.0/3.0))*(log(polymorph_activity/CC_activity/Caq));
cout << log(polymorph_activity/CC_activity/Caq)<<"  " <<Caq<< " \n ";
}
int delme;

cout << gamma[0] << "   " << gamma[1] << "   "<< gamma[2] << "   "<< gamma[3] << " \n  ";
cout <<rhos[0] << "   " << rhos[1] << "   "<< rhos[2] << "   "<< rhos[3] << " \n  ";
cin >> delme;
return;
}






void  mom_Aq::transport_moments_and_eta( double dt, vector<double> dxML) {

    double dCaq=0;
    double dt_max=0.0;
    double temp=0;
    double factor=.1;
    vector<double> etaDot(odtlp->neta,0.0);
    vector<double> momx  (odtlp->nmom,0.0);
    vector<vector<double> > flxEta(odtlp->neta, vector<double>(odtlp->ngrdf,0.0));
    vector<double> Yelem(6,0);
    double dt_0=dt;
    double dd=0;
    int    ip;


          ETApx->getEtaDiffusiveFluxes(odtlp->eta, 1000.0, dxML, flxEta, true );  // This also populates the vector etaObj.CaCO3_gamma[j];
    for(int j=0; j < odtlp->ngrd; j++) {

        dt=dt_0;

        for (int k=0; k<odtlp->nmom; k++) // Derek
            momx[k]=odtlp->mom[k][j]*odtlp->rho[j];  // This removes the rho

        while(true){

            ip=j+1;
            dd = 1.0/(dxML[j]*odtlp->rho[j]);


            dCaq= get_moment_source(momx, ETApx->CaCO3[j],  odtlp->eta[odtlp->neta-1][j], ETApx->CaCO3_gamma[j]); // Compute Source term (no fluxes, for moments)

            for(int i=0; i<odtlp->neta; i++) {          // Compute d eta / dt  diffusion + Source
                etaDot[i] = -dd*(flxEta[i][ip] - flxEta[i][j]);  

                etaDot[i]+=dCaq/odtlp->rho[j]*ETApx->Yspel[8][i];// Yspel converts mass of species to mass of element, species 8 is CaCO3
            } 

                 dt_max= dt;
            for (int i=0; i< odtlp->nmom; i++){             //Determine the maximium allowable step size   // adapt time-wise on mom
                temp  = abs( momx[i] /  Source[i]  *factor ) ;
                dt_max =  min(temp,dt_max) ;
         //       cout << dt_max << "  " << temp << " \n";
            } 
     //       for (int i=0; i< odtlp->neta-1; i++){             //Determine the maximium allowable step size  // adapt time-wise on eta
     //           temp   = abs( odtlp->eta[i][j] /  etaDot[i] * factor ) ;
     //           dt_max =  min(temp,dt_max) ;
     //          }

            for (int i=0; i< odtlp->nmom; i++)      // transport moments: explicit Euler
                momx[i] = momx[i] + dt_max*Source[i];




            for (int i=0; i< odtlp->neta; i++)      // transport eta: explicit Euler
                odtlp->eta[i][j] = odtlp->eta[i][j] + dt_max*etaDot[i];


            if (dt_max == dt){
                for (int i=0; i<odtlp->nmom; i++) // Derek
                    odtlp->mom[i][j]= momx[i]/odtlp->rho[j];  // Returns the rho

                break;
            }
            dt-=dt_max;

            for(int i=0; i<odtlp->neta-1; i++)
            Yelem[i] = odtlp->eta[i][j];


            ETApx->equilibrate(Yelem, 1000.0,ETApx->TA,false);  // This also populates the vector etaObj.CaCO3_gamma[j];
            ETApx->CaCO3[j]=ETApx->Conc[8];
            ETApx->CaCO3_gamma[j]=ETApx->gamma[8];
  //      for (int i=0; i<16; i++)
   //     cout << Source[i] << "  ";
    //    cout << " \n";
    //    cout << dt_max << " \n";
        }


    }
    return;
}

         void mom_Aq::computeFluxes(std::vector<double> &dd) {return;}
         void mom_Aq::setVolumeFraction(std::vector<double> &volFrac){return;} 

         void mom_Aq::computeSourceTerms(){                 // This isn't computing the source, in fact, it is transporting ETA and MOM using an adaptive first order euler method
             transport_moments_and_eta( diff->Class_dtStep , diff->dxML) ;
          // diff->rhsSrc[diff->iptEta + k][i]
             return; 
         }  





