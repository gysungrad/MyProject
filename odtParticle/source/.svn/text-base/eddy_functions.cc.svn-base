
#include "eddy_functions.h"
#include <iostream>


using namespace std;
void eddy_functions::apply_single_map(anyline &anyl, double leftEdge, double rightEdge){

   /****************************************************************************************************************
   *Adapt the constants cFrac1 and cFrac2 to get another version of the triplet map                                *
   *cfrac1 is the size relative to the map size for the first and third segment of the map, cfrac for the second.  *
   ****************************************************************************************************************/

    int    ngrd  = anyl.ngrd;
    int    ngrd2 = ngrd*2;
    int    ngrd3 = ngrd*3;
    double cFrac1   = 1./3.;
    double cFrac2   =  1./3.;
    int i,j,k;

    ///////////// make space for 2nd, 3rd segments
    vector<double> oldposf=anyl.posf;
    for(i=0;i<oldposf.size();i++){
     oldposf[i]-=leftEdge;
    }
    anyl.pos.resize(ngrd3);
    anyl.posf.resize(ngrd3+1);
    anyl.rho.resize(ngrd3);
    anyl.molec.resize(ngrd3);
    anyl.lambda.resize(ngrd3);
    anyl.phase.resize(ngrd3);

    for(k=0; k<anyl.nprops; k++) 
        (*anyl.props[k]).resize(ngrd3); 
//     for(k=0; k<anyl.nauxVar; k++) 
//         (**(anyl.auxVar+k)).resize(ngrd3);  

    ///////////// fill second segment
    //D.M.comments
    //take the values from the still existing first segment and
    //fill them (last entry into first entry, second last entry into second entry etc)

    for(i=ngrd, j=ngrd-1; i<ngrd2; i++, j--) {
        anyl.rho[i]  = anyl.rho[j];
        anyl.molec[i] = anyl.molec[j];
        anyl.lambda[i] = anyl.lambda[j];
        anyl.phase[i] = anyl.phase[j];
        for( k=0; k<anyl.nprops; k++ )
            (*anyl.props[k])[i] = (*anyl.props[k])[j];
// 	for( k=0; k<anyl.nauxVar; k++ )
//             (**(anyl.auxVar+k))[i] = (**(anyl.auxVar+k))[j];
    }

    ///////////// fill third segment
    //D.M. segment
    //copy first segment into last segment.
    //chance of speed improvement?
    for(i=ngrd2, j=0; i<ngrd3; i++, j++) {
        anyl.rho[i]  = anyl.rho[j];
        anyl.molec[i] = anyl.molec[j];
        anyl.lambda[i] = anyl.lambda[j];
        anyl.phase[i] = anyl.phase[j];
        for( k=0; k<anyl.nprops; k++ )
            (*anyl.props[k])[i] = (*anyl.props[k])[j];
// 	for( k=0; k<anyl.nauxVar; k++ )
//             (**(anyl.auxVar+k))[i] = (**(anyl.auxVar+k))[j];
    }

    //////////// write new cell and face positions

    //--------- cell face positions

    double newSegSize1 = (rightEdge-leftEdge)*cFrac1;
    double newSegSize2 = (rightEdge-leftEdge)*cFrac2;

    anyl.posf[ngrd3] = anyl.posf[ngrd];       // last face doesn't change

    for(i=1; i<ngrd; i++)
        anyl.posf[i] = leftEdge + (anyl.posf[i]-leftEdge)*cFrac1;
    anyl.posf[ngrd] = leftEdge + newSegSize1;

//     double secondsegmentbegin=anyl.posf[ngrd];
//     for(i=ngrd+1,j=1;i<ngrd2; i++,j++)
//         anyl.posf[i] = oldposf[j]*cFrac2+secondsegmentbegin;
//     anyl.posf[ngrd2] = leftEdge + newSegSize2+newSegSize1;
// 
//     double thirdsegmentbegin=anyl.posf[ngrd2];
//     for(i=ngrd2+1,j=1; i<ngrd3; i++,j++)
//         anyl.posf[i] = oldposf[j]*cFrac1+thirdsegmentbegin;
    for(int i=ngrd+1, j=ngrd-1; i<ngrd2; i++, j--) 
        anyl.posf[i] = anyl.posf[ngrd] + (anyl.posf[ngrd]-anyl.posf[j]);
        anyl.posf[ngrd2] = anyl.posf[ngrd] + newSegSize2;                

    for(int i=ngrd2+1, j=ngrd2-1; i<ngrd3; i++, j--)
        anyl.posf[i] = anyl.posf[ngrd2] + (anyl.posf[ngrd2]-anyl.posf[j]);

    //--------- cell center positions

    for(i=0; i<ngrd3; i++)
        anyl.pos[i] = 0.5*(anyl.posf[i]+anyl.posf[i+1]);



    anyl.ngrd  = ngrd3;
    anyl.ngrdf = anyl.ngrd+1;



    }


///////////////////////////////////////////////////////////////////////////////
