/**
 * @file adaptMesh.cc
 * Header file for class adaptMesh
 */

#include "adaptMesh.h"
#include "odtParam.h"
#include "processor.h"
#include <fstream>
#include <cmath>            
#include <cstdlib>            
#include <algorithm>            


using namespace std;

extern processor proc;


///////////////////////////////////////////////////////////////////////////////
/*! adaptMesh constructor function
 *
 * @param anyll  \input set anyline pointer with.
 * @param odtPP  \input set odtParam pointer with.
 * @param phiP   \input set vector pointer with.
 * @param Ladapt \input enables the splitLargeCells() funktion
 */

adaptMesh::adaptMesh(anyline *anyll, odtParam *odtPP, vector<vector<double> *> phiP, bool Ladpt) {
#ifdef COMPSCI
    init(anyll,odtPP,Ladpt);
#else
    anyl = anyll;
    odtP = odtPP;

    ngrd  = anyl->ngrd;
    ngrdf = anyl->ngrdf;

    iLower = 0;
    iUpper = ngrd-1;
    posLower = anyl->posf[iLower];
    posUpper = anyl->posf[iUpper+1];

    xnf = anyl->posf;
    setDxArray();

    // added by Falko; call splitLargeCells only during initialisation
    // For a bad choice of gDens and max grid spacing this call of splitting
    // large cells for an eddy can result in non-increasing odtline faces.
    if(Ladpt)
        splitLargeCells();
    
    bdy = &(anyl->phase);
#endif
    phi = phiP;
}



///////////////////////////////////////////////////////////////////////////////
/** User interface to the mesh adapter.
 *
 * @param iLowerDummy \input lower cell to adapt from (approx).
 * @param iUpperDummy \input upper cell to adapt to (approx).
 */

void adaptMesh::adaptGrid(int iLowerDummy, int iUpperDummy) {
    
    for (int i=0; i < anyl->ngrd; i++){
        if(fabs(anyl->posf.at(i) - anyl->posf.at(i+1)) < 1.0e-15){
            *proc.ostrm << endl << "   erasing cell i=" << i;
            anyl->eraseCell(i);
            i-=3; // i-=2; // -2 should be enough
        }
        else if(anyl->posf.at(i) > anyl->posf.at(i+1)){
            *proc.ostrm << endl << "ERROR: Non-increasing odtLine!";
            anyl->outputProperties("ERROR_anyline.dat");
            exit(0);
        }
    }
    if(!odtP->LmultiPhase){
        // case is for a single phase
        if(iLowerDummy == iUpperDummy)
            return;
        
        if(iLowerDummy < iUpperDummy) 
            adaptGrid_details(iLowerDummy, iUpperDummy);
        
        //---------- for periodic where adaption region wraps, adapt twice: once on each section
        
        else {
            
            posUpper = anyl->posf[iUpperDummy+1];
            
            adaptGrid_details(iLowerDummy, anyl->ngrd-1);
            
            iUpperDummy = anyl->linePositionToIndex(posUpper, false);
            
            adaptGrid_details(0, iUpperDummy);
        }
    }
    else{ // case is for multiple immiscable phases
        if(iLowerDummy == iUpperDummy)
            return;
        if(iLowerDummy < iUpperDummy){
            
            double posUpperDummy = anyl->posf[iUpperDummy+1];
            double posMiddDummy;
            double phase = anyl->phase[iLowerDummy];
            
            for (int iDummy = iLowerDummy+1; iDummy <= iUpperDummy; iDummy++){
                if (phase != anyl->phase[iDummy]){
                    
                    // change to new Phase
                    phase = anyl->phase[iDummy];
                    
                    // save position of phase change
                    posMiddDummy = anyl->posf[iDummy];
                    
                    // adapt grid within the old phase
                    adaptGrid_details(iLowerDummy, iDummy-1);
                    
                    // recover indexes
                    iUpperDummy = anyl->linePositionToIndex(posUpperDummy, false);
                    iLowerDummy = anyl->linePositionToIndex(posMiddDummy, true);
                    iDummy = iLowerDummy;
                }
            }
            
            // adaption of the last phase
            adaptGrid_details(iLowerDummy, iUpperDummy);
        }
        
        //---------- for periodic where adaption region wraps, adapt twice: once on each section
        else{
            posUpper = anyl->posf[iUpperDummy+1];
            
            // call adaptGrid twice, once for the lower part and once for the
            // upper part. adaptGrid is needet to take different phases into
            // account.
            adaptGrid(iLowerDummy, anyl->ngrd-1);
            
            iUpperDummy = anyl->linePositionToIndex(posUpper, false);
            
            adaptGrid(0, iUpperDummy);
        }
    }

}


///////////////////////////////////////////////////////////////////////////////
/** 
 * Makes several passes over the grid applying different criteria in succession.
 *
 * @param iLowerDummy \input lower cell to adapt from (approx).
 * @param iUpperDummy \input upper cell to adapt to (approx).
 */
void adaptMesh::adaptGrid_details(int iLowerDummy, int iUpperDummy) {


    // set xf and yf in the region to adapt

    iLower   = iLowerDummy;
    iUpper   = iUpperDummy;
    posLower = anyl->posf[iLower];
    posUpper = anyl->posf[iUpper+1];

    xf = vector<double>(anyl->posf.begin()+iLower,    // face positions
            anyl->posf.begin()+iUpper+2); 

    yf = vector<vector<double> >(phi.size(), vector<double>(xf.size(), 0.0));

    for(int iProf = 0; iProf<(int)phi.size(); iProf++) {

        if(iLower==0)                            
            yf[iProf][0] = (*phi[iProf])[0];
        else 
            yf[iProf][0] = (*phi[iProf])[iLower-1] + (anyl->posf[iLower]-anyl->pos[iLower-1])/
                (anyl->pos[iLower]-anyl->pos[iLower-1])*
                ((*phi[iProf])[iLower]-(*phi[iProf])[iLower-1]);
        if(iUpper==anyl->ngrd-1)
            yf[iProf][xf.size()-1] = (*phi[iProf])[iUpper];
        else
            yf[iProf][xf.size()-1] = (*phi[iProf])[iUpper] + (anyl->posf[iUpper+1]-anyl->pos[iUpper])/
                (anyl->pos[iUpper+1]-anyl->pos[iUpper])*
                ((*phi[iProf])[iUpper+1]-(*phi[iProf])[iUpper]);
        for(int i=1,j=iLower+i; i<(int)xf.size()-1; i++, j++) {          
            yf[iProf][i] = (*phi[iProf])[j-1] + (anyl->posf[j]-anyl->pos[j-1])/
                (anyl->pos[j]-anyl->pos[j-1])*
                ((*phi[iProf])[j]-(*phi[iProf])[j-1]);
        }
    }
    //---------- Scale xf and yf

    //double xfMax = *max_element(xf.begin(), xf.end());
    //double xfMin = *min_element(xf.begin(), xf.end());
    double xfMax = *max_element(anyl->posf.begin(), anyl->posf.end());
    double xfMin = *min_element(anyl->posf.begin(), anyl->posf.end());
    for(int i=0; i<(int)xf.size(); i++) 
        xf[i] = (xf[i]-xfMin)/(xfMax-xfMin);

    for(int iProf = 0; iProf<(int)phi.size(); iProf++) {
        //double yfMax = *max_element(yf[iProf].begin(), yf[iProf].end());
        //double yfMin = *min_element(yf[iProf].begin(), yf[iProf].end());
        double yfMax = *max_element(phi[iProf]->begin(), phi[iProf]->end());
        double yfMin = *min_element(phi[iProf]->begin(), phi[iProf]->end());
        for(int i=0; i<(int)xf.size(); i++) 
            yf[iProf][i] = (yf[iProf][i]-yfMin)/(yfMax-yfMin);
    }


    //---------- Compute the new grid distance function

    vector<double> s_dist(xf.size(),0.0);                 // distance function for x,y
    vector<double> sn_dist;                                // new distance function for new grid X,Y

    double s_distCum = calcDistance(xf, yf, s_dist);   

    int nNewGrd = odtP->gDens*s_distCum+1;      
    sn_dist.resize(nNewGrd+1, 0.0);
    xnf.resize(nNewGrd+1, 0.0);

    double dS = s_distCum / nNewGrd;
    sn_dist[0] = 0.0;
    for(int i=1, im=0; i<nNewGrd+1; i++, im++)
        sn_dist[i] = sn_dist[im] + dS;

    //---------- Get the new grid using the old grid, old dist func, new dist func

    interpVec(s_dist, xf, sn_dist, xnf);                  // interp from xf --> xnf

    for(int i=0; i<(int)xnf.size(); i++)                      // Unscale xnf
        xnf[i] = xfMin + xnf[i]*(xfMax-xfMin);

    xnf[0] = posLower;
    xnf[xnf.size()-1] = posUpper;
    for(int i=0; i<(int)xnf.size(); i++)                      // doldb ?
        if(xnf[i] < posLower || xnf[i] > posUpper) {
            cout << endl << "Error in mesher: X is out of bounds, i= " << i << endl;
            cout << "        posLower, posUpper, xnf[i] ,xnf[0], xnf[end] = " << 
                posLower << ", " << posUpper << ", " << xnf[i]<< ", " << xnf[0]<< ", " << xnf[xnf.size()-1] << endl;
            exit(0);
        }

    for(int i=1; i<(int)xnf.size(); i++)                     // doldb ?
        if(xnf[i] <= xnf[i-1]) {
            anyl->outputProperties("Line_before_error.dat");
            cout << endl << endl;
            cout << "i = " << i << "  xnf.size = " << xnf.size() << "  xnf = "
                 << xnf[i-1] <<" "<< xnf[i] << endl;
            cout << endl << "Error in mesher: xnf is nonincreasing " << endl;
            exit(0);
        }

    //--------- Add in the rest of the old grid so that we can adapt desired region of WHOLE line

    if(iLower > 0) 
        xnf.insert(xnf.begin(), anyl->posf.begin(), anyl->posf.begin()+iLower);
    if(iUpper < anyl->ngrd-1) 
        xnf.insert(xnf.end(),   anyl->posf.begin()+iUpper+2, anyl->posf.end());

    /////////////////////////////////////////////////////////////////////////////////////////////

    //---------- Adapt the new grid 

    ngrdf = xnf.size();
    ngrd  = ngrdf-1;
    setDxArray();
    set_iLower_iUpper();
    
    //----------- split large cells
    for (int i=0; i<ngrd; i++) {
        if (dx[i] > odtP->dxmax){
            dx[i] *= 0.5;
            dx.insert(dx.begin()+i,dx[i]);
            ngrd++;
            ngrdf++;
            i--;
        }
    }
    xnf.resize(ngrdf);
    xnf[0] = anyl->posf[0];               // recover the grid
    for(int i=1; i<ngrdf; i++)
        xnf[i] = xnf[i-1]+dx[i-1];


    //----------- Small Cells
    //if (!odtP->LmultiPhase)
    //    mergeSmallCells();   // modifies dx, ngrd, ngrdf (but xnf is now inconsistent)
    //else
        mergeSmallCellsMP(); // Function does the same as the old one, but has a phase check

    xnf.resize(ngrdf);
    xnf[0] = anyl->posf[0];               // recover the grid
    for(int i=1; i<ngrdf; i++)
        xnf[i] = xnf[i-1]+dx[i-1];

    set_iLower_iUpper();        // recover the bounds

    //----------- 2.5 rule

    impose2point5rule();        // modifies dx, ngrd, ngrdf (but xnf is now inconsisten)

    xnf.resize(ngrdf);
    xnf[0] = anyl->posf[0];               // recover the grid
    for(int i=1; i<ngrdf; i++)
        xnf[i] = xnf[i-1]+dx[i-1];

    //---------- Apply new grid to odt line:

    //---------- split anyl where new grid has faces that don't match anyl.
    // anyl     |          |          |          |          |          |
    // xnf      |     ::   |          |    : :                   :     |
    // anyl new |     ::   |          |    : :   |          |    :     |

    vector<double> whereSplit;
    int iCell; //, iCellD;  //  !!!!!  currently unused variable

    double rtol = 1.0E-10*odtP->domainLength/anyl->ngrd;
    for(int j=1; j<(int)xnf.size()-1; j++) {         // loop over new grid cell faces
        iCell = anyl->linePositionToIndex(xnf[j], true);
        //if(anyl->posf[iCell] != xnf[j]) {
        if(abs(anyl->posf[iCell] - xnf[j]) > rtol) {
            //if(anyl->posf[iCell+1] != xnf[j]) {
            if(abs(anyl->posf[iCell+1] - xnf[j]) > rtol) { // far away from next face
                // this second if is needed due to the fact that somtimes there is
                // a face in xnf that's exact equal to one face in anyl->posf. In
                // this case a cell with dx=0 is generated which may cause trouble.
                whereSplit.push_back(xnf[j] - anyl->posf[iCell]);
                j++;
                for( ; j<(int)xnf.size()-1; j++) {
                    //iCellD = anyl->linePositionToIndex(xnf[j], true);
                    //if(iCellD != iCell) { 
                    //    j--;
                    //    break;
                    //}
                    //else
                    //    whereSplit.push_back(xnf[j] - anyl->posf[iCell]);
                    
                    if( (abs(anyl->posf[iCell+1] - xnf[j]) > rtol) 
                        && (anyl->posf[iCell+1] > xnf[j]) ) {

                        whereSplit.push_back( xnf[j] - anyl->posf[iCell] );
                    }
                    else{
                        j--;
                        break;
                    }
                }
                anyl->splitCell(iCell, whereSplit.size(), whereSplit, true);
                whereSplit.clear();
            }
        }
    }

    //---------- merge anyl where faces present that are not in new grid
    // anyl     |     ||   |          |    | |   :          :    |     |
    // xnf      |     ||   |          |    | |                   |     |
    // anyl new |     ||   |          |    | |                   |     |
    // now anyl matches xnf

    for(int i=1; i<anyl->ngrdf-1; i++) {         // loop over interior faces of anyl
        //if(anyl->posf[i] != xnf[i]) {
        if( abs(anyl->posf[i]-xnf[i]) > rtol) {
            anyl->merge2cells(i-1);
            i--;
        }
    }

}


/////////////////////////////////////////////////////////////////////////////
/** Finds iLower and iUpper based on the new cell face positions
 *
 */
void adaptMesh::set_iLower_iUpper() {

    iLower = 0;
    iUpper = ngrd-1;

    if(posLower==xnf[ngrd] || posUpper==xnf[0]) {
        cout << posLower << " " << posUpper << " " << xnf[0] << " " << xnf[ngrd] << endl;
        cout << endl << "ERROR in adapMesh posLower/Upper on wrong boundary" << endl;
        exit(0);
    }
    for(int i=0; i<ngrd; i++) {
        if(posLower < xnf[i+1]) {
            iLower = i;
            break;
        }
    }
    for(int i=ngrd-1; i>=0; i--) {
        if(posUpper > xnf[i]) {
            iUpper = i;
            break;
        }
    }

}


///////////////////////////////////////////////////////////////////////////////
/** Cells that are too small are merged with their neighbors.
 *  Cells here are marked simply as whether they are small or not.
 *  The subsequent merging is with either neighbor.  Elsewhere marked
 *  cells are to be merged with the right neighbor.
 *  Small cells are usually generated by triplet map events that
 *  compress the grid (or by accelerating flows in spatial domains that compress the grid).
 *  Mark the small cells by filling mark (member) with indicies in grid
 *  Merge the small cells in size order.  That is, initial size order, not
 *  accounting for intermediate changes.  This eliminates directional bias too.
 *  So rearrange the small cell indicies into sorted order, scMark (local).
 *  Then loop over each small cell till big enough, merging with its smaller neighbor.
 *  The smaller neighbor may or may not be a small cell.  When merge cells,
 *  need to decrement larger indicies (cells) in the scMark array, and delete cells
 *  from scMark if a small cell is merged with another.
 *  No special treatment for periodic boundaries.  Periodic is the same as nonperiodic.
 *  The periodic boundary position is maintained, and no merging occurs across 
 *  a periodic boundary (which would move the boundary from zero).                               
 */
void adaptMesh::mergeSmallCells() {

    int i, j;

    //---------- mark the small cells

    mark.clear();
//    if(!odtP->Lperiodic) 
//        for(i=iLower; i<=iUpper; i++)    
//            if(dx[i] < odtP->dxmin)
//                mark.push_back(i);
//    else 
        for(i=0; i<=ngrd-1; i++)            // periodic: the combustion line can have
            if(dx[i] < odtP->dxmin)         // small cells outside the adapt range.
                mark.push_back(i);          // just check the whole domain

    //---------- sort the small cells to remove any directional bias
    
    vector<int> ind(mark.size());           // index map for the sort
    for(i=0; i<(int)ind.size(); i++)        // sort "ind" using dx as sort condition
        ind[i]=i;
    sort(ind.begin(), ind.end(), *this);    // *this invokes the functor "operator()"
    vector<int> scMark(mark.size());        // reorder the mark array
    for(i=0; i<(int)ind.size(); i++)        // could just use mark[ind[i]],
        scMark[i] = mark[ind[i]];           // but reorder for simplicity 


    //---------- merge small cells
    
    int isc;                                // small cell index
    int iSmallerNB;                         // smaller of 2 neighbors
    int iend = scMark.size();               // may change as merge cells
    bool LcellDone;                         // flag to repeat small cells till big

    //---------------------------------

    for(i=0; i<iend; i++) {                 // loop over small cells
        
        isc = scMark[i];
        isc = scMark[i];
        LcellDone = false;     

        //---------------------------------

        while(!LcellDone) {                 // keep going till small cell is done

            if(isc == 0)
                iSmallerNB = 1;
            else if(isc == ngrd-1)
                iSmallerNB = ngrd-2;
            else {
                iSmallerNB = (dx[isc-1] <= dx[isc+1]) ? isc-1 : isc+1;

                // For a location with two neighbors, if the phase of the chosen               
                // neighbor is different, instead choose the other neighbor.                   

                //if (fabs((*bdy)[isc]-(*bdy)[iSmallerNB]) > 0.01)                               
                //    iSmallerNB = 2 * isc - iSmallerNB;                                         
            }                                                                              

            if(dx[iSmallerNB] + dx[isc] >= odtP->dxmin) 
                LcellDone = true;           // if new cell big enough then done 

            //---------------------------------

            if(iSmallerNB > isc) {

                // merge cells

                ngrdf--;
                ngrd--;
                dx[isc] = dx[isc] + dx[isc+1];
                dx.erase( dx.begin() + isc+1 );

                // delete cell iSmallerNB from scMark array (if present)
                // and decrement all cells > isc in scMark

                for(j=i+1; j<iend; j++) 
                    if(scMark[j] > isc)
                        scMark[j]--;
                for(j=i+1; j<iend; j++) 
                    if(scMark[j]==isc){
                        scMark.erase(scMark.begin()+j);
                        iend--;
                        break;
                    }

            }         // if(iSmallerNB > isc)

            //---------------------------------

            else {    

                // merge cells

                ngrdf--;
                ngrd--;
                dx[isc-1] = dx[isc-1] + dx[isc];
                dx.erase( dx.begin() + isc );
                 
                // decrement isc and scMark[i] which is isc
                // delete cell isc-1 in scMark if present
                // and decrement 

                scMark[i]--;
                isc--;

                for(j=i+1; j<iend; j++) 
                    if(scMark[j] > isc)
                        scMark[j]--;
                for(j=i+1; j<iend; j++) 
                    if(scMark[j]==isc){
                        scMark.erase(scMark.begin()+j);
                        iend--;
                        break;
                    }

            }     // else (i.e., iSmallerNB < isc)
            //---------------------------------
        }         // end while
        //---------------------------------
    }             // end loop over small cells

    //---------- update the eddy regions

}                 // end function



//////////////////////////////////////////////////////////////////////////////////
/** This function does the same as mergeSmallCell() excapt that it checks, in
 *  addition, the phase of the merged cell. If there is no opportunity for
 *  merging the cell, it is skipped.
 *  
 *  ToDo:
 *  If there is a too small cell and both neigbor cells have different phases,
 *  the cell has to be handeled. For this case, I have no idea how to handle.
 */
void adaptMesh::mergeSmallCellsMP() {
    
    double x  = anyl->posf[0];
    double xp = 0; // midd point of previous cell
    double xc = 0; // midd point of current cell
    double xn = 0; // midd point of next cell
    double phase_p = 0;
    double phase_c = 0;
    double phase_n = 0;
    
    
    for (int i=0; i<ngrd; i++){
        if (dx[i] >= odtP->dxmin)
            x += dx[i];
        else{
            // the current cell is too small
            
            // if it is the first cell
            if (i == 0){
                xc = x + dx[i]/2;
                xn = x + dx[i]+dx[i+1]/2;
                phase_c = anyl->phase[ anyl->linePositionToIndex( xc, true ) ];
                phase_n = anyl->phase[ anyl->linePositionToIndex( xn, true ) ];
                
                if (phase_c == phase_n){
                    ngrdf--;
                    ngrd--;
                    dx[i] = dx[i] + dx[i+1];
                    dx.erase( dx.begin() + i+1 );
                    i--;
                } else{
                    x += dx[i];
                }
                continue;
            }
            
            // if it is the last cell
            if (i == ngrd-1){
                xp = x - dx[i-1]/2;
                xc = x + dx[i]/2;
                phase_p = anyl->phase[ anyl->linePositionToIndex( xp, true ) ];
                phase_c = anyl->phase[ anyl->linePositionToIndex( xc, true ) ];
                
                if (phase_p == phase_c){
                    ngrdf--;
                    ngrd--;
                    x -= dx[i-1]; // subtracting the dx of the previous cell, added in the last loop step
                    dx[i-1] = dx[i-1] + dx[i];
                    dx.erase( dx.begin() + i );
                    i--; i--;
                } else{
                    x += dx[i];
                }
                continue;
            }
            
            // it is a cell in the midd of the domain
            xp = x - dx[i-1]/2;         // cell center of previous cell
            xc = x + dx[i]/2;           // cell center of current cell
            xn = x + dx[i] + dx[i+1]/2; // cell center of next cell
            phase_p = anyl->phase[ anyl->linePositionToIndex( xp, true ) ];
            phase_c = anyl->phase[ anyl->linePositionToIndex( xc, true ) ];
            phase_n = anyl->phase[ anyl->linePositionToIndex( xn, true ) ];
            
            // the case, that the previous cell has another phase
            if (phase_p != phase_c && phase_c == phase_n){
                // the next cell has to be taken for merging cells
                ngrdf--;
                ngrd--;
                dx[i] = dx[i] + dx[i+1];
                dx.erase( dx.begin() + i+1 );
                i--;
                continue;
            }
            
            // the case, that the next cell has another phase
            if (phase_p == phase_c && phase_c != phase_n){
                // the previous cell has to be taken for merging cells
                ngrdf--;
                ngrd--;
                x -= dx[i-1]; // subtracting the dx of the previous cell, added in the last loop step
                dx[i-1] = dx[i-1] + dx[i];
                dx.erase( dx.begin() + i );
                i--; i--;
                continue;
            }
            
            // the case, that both neighbor cells have different phases
            if (phase_p != phase_c && phase_c != phase_n){
                x += dx[i];
                continue;
            }
            
            // the last case, that all three cells have the same phase
            if (phase_p == phase_c && phase_c == phase_n){
                // the smaller one of both neighbor cells has to be taken
                if (dx[i-1] < dx[i+1]){
                    // the previous neigbor has to be taken
                    ngrdf--;
                    ngrd--;
                    x -= dx[i-1]; // subtracting the dx of the previous cell, added in the last loop step
                    dx[i-1] = dx[i-1] + dx[i];
                    dx.erase( dx.begin() + i );
                    i--; i--;
                } else{
                    // the next neigbor has to be taken
                    ngrdf--;
                    ngrd--;
                    dx[i] = dx[i] + dx[i+1];
                    dx.erase( dx.begin() + i+1 );
                    i--;
                }
                continue;
            } else{
                *proc.ostrm << "\n\ni = " << i;
                *proc.ostrm <<   "\nphase_p = " << phase_p;
                *proc.ostrm <<   "\nphase_c = " << phase_c;
                *proc.ostrm <<   "\nphase_n = " << phase_n;
                *proc.ostrm <<   "\n(phase_p != phase_c && phase_c == phase_n) = " << (phase_p != phase_c && phase_c == phase_n);
                *proc.ostrm <<   "\n(phase_p != phase_c && phase_c == phase_n) = " << (phase_p == phase_c && phase_c != phase_n);
                *proc.ostrm <<   "\n(phase_p != phase_c && phase_c == phase_n) = " << (phase_p != phase_c && phase_c != phase_n);
                *proc.ostrm <<   "\n(phase_p != phase_c && phase_c == phase_n) = " << (phase_p == phase_c && phase_c == phase_n);
                *proc.ostrm << "\nERROR: In the function adaptMesh::mergeSmallCells_FS()"
                << endl << "This error might not be there caused by logic one of the"
                << endl << "previous cases have to occur."
                << endl << "There must be a NaN problem.";
                exit(0);
            }
        }
    }
}


//////////////////////////////////////////////////////////////////////////////////
/**
 * Imposes the 2.5 rule for Mesh adaption to the grid.
 */
void adaptMesh::impose2point5rule(){

    int    i;
    double d1;

    //---------- set the mark array with offending cells

    mark.clear();
    int iHi = (iUpper < ngrd-1) ? iUpper   : iUpper-1;
    int iLo = (iLower > 0)      ? iLower-1 : iLower;    
    //int iHi = ngrd-2;
    //int iLo = 0;
    for(i=iLo; i<=iHi; i++) {                      
        d1 = dx[i]/dx[i+1];
        if(d1 > 2.5 || d1 < 0.4)
            mark.push_back(i);
    }
    i = ngrd-1;
    if(odtP->Lperiodic && iUpper==i) {
        d1 = dx[i]/dx[0];
        if(d1 > 2.5 || d1 < 0.4)
            mark.push_back(i);
    }

    //-------- loop over each marked cell and fix the 2.5 offenders

    for(i=0; i<(int)mark.size(); i++)        
        fix2point5offender(mark[i], i);   

}

///////////////////////////////////////////////////////////////////////////////
/** Marks cells that offend their neighbors with respect to the 2.5 rule
 *
 *           Routine is recursive.  Returns if the test is passed, or if we are at
 *                the ends of the domain.                                                \n 
 *           If it starts with no "small cells" routine, it won't create any small cells.\n
 *           The order of operations doesn't matter.                                     \n
 *           Find the larger of mPos and its next neighbor (mPos+1) and split that
 *                cell on the side of the smaller cell.  Split it nsplit times in half.  \n
 *                So, <tt><pre>
 *                      |               ||  ==>  |   *   :    :   || </pre></tt>
 *                when splitting twice.\n
 *           Then go to the opposite side of the split cell (marked with \c * above) and
 *                compare to its neighbor on the other side.  In this example, we would
 *                "mark" the cell before the \c * and call the same routine on it.          \n
 *           Hence, we can traverse through the domain.  Generally you keep traversing 
 *                till satisfy the rule or hit the edge.  To keep going in one direction
 *                the cells need to keep increasing (so we won't go too far since motion
 *                implies geometric growth).  You also will stop traversing when you run into
 *                a smaller cell.                                                        \n
 *           Consider:\n <tt><pre>
 *                      ||        |            ||            ||                                  |
 *                       0    1          2     3       4     5                  6
 *           mark offenders: 0 2 3 4 5
 *           Start at 0:
 *           0 vs 1 --> split 1 once
 *                      ||    :   |            ||            ||                                  |
 *                       0  1   2        3     4       5     6                  7
 *           2 vs 3 --> split 3 once
 *                      ||    :   |      :     ||            ||                                  |
 *                       0  1   2     3     4  5       6     7                  8
 *           4 vs 5 --> split 4 once
 *                      ||    :   |      :  :  ||            ||                                  |
 *                       0  1   2     3   4  5 6       7     8                  9
 *           3 vs 4 --> pass --> done
 *           now do next in mark array (which was updated as we went to ( X 3 6 7 8)
 *           Start at 3:
 *           3 vs 4 --> pass (trivially) --> done
 *           Start at 6:
 *           6 vs 7 --> split 7 once
 *                      ||    :   |      :  :  ||      :     ||                                  |
 *                       0  1   2     3   4  5 6   7      8  9                  10
 *           8 vs 9 --> pass (say) --> done:  mark array ( X X X 7 9 )
 *           Start at 7:
 *           7 vs 8 --> pass
 *           Start at 9:
 *           9 vs 10 --> split 4 times
 *                      ||    :   |      :  :  ||      :     || :  :    :       :                |
 *                       0  1   2     3   4  5 6   7      8  9 10 11 12    13           14
 *           14 at the end --> DONE
 *           </pre></tt>
 *
 * @param mPos    \input  marks a cell that offends its next neighbor with respect to the 2.5 rule.
 * @param iglobal \input is the index of which loop it is on in the mark arrary,
 *       for the mark array update. (Could just leave this out and do all of them
 *       but only need to update subsequent ones).
 *
 */
void adaptMesh::fix2point5offender(int mPos, const int &iglobal) {


    //--------- Simply return if you pass the test

    int ip = (mPos == ngrd-1) ? 0 : mPos+1;
    double ratio = (dx[mPos] < dx[ip]) ? dx[ip]/dx[mPos] : dx[mPos]/dx[ip];
    if(ratio < 2.5)
        return;

    //--------- Split the larger of the two cells

    int     i;
    int     isplt = (dx[mPos] < dx[ip]) ? ip : mPos;          // split larger of 2 cells 
    int     nsplt = static_cast<int>(log2(ratio / 2.5)) + 1 ; // how many splits to do

    vector<double> icp(nsplt);                 // pos of new intrnl cell fcs (rel to lft edge)

    bool splitCellOnRight = (isplt > mPos);
    if(mPos == ngrd-1 && isplt == 0)           // for periodic
        splitCellOnRight = true;

    if(splitCellOnRight) {                     // | |  :  :    :        |  --> nsplt=3
        icp[nsplt-1] = dx[isplt]*0.5;
        for(i=nsplt-2; i>=0; i--)
            icp[i] = icp[i+1]*0.5;

        //splitCell(isplt, nsplt, icp);               

        for(i=nsplt-1; i>=1; i--)              // update dx arr, reuse icp arr
            icp[i] -= icp[i-1];                // i.e. we inserted cells, so insert dx
        dx[isplt] *= 0.5;                      // icp was positions, now are dx's
        dx.insert(dx.begin()+isplt, icp.begin(), icp.end());
    }
    else {                                     // |        :    :  :  | |   --> nsplt=3
        icp[0] = dx[isplt]*0.5;                
        for(i=1; i<nsplt; i++)                 // set icp as dx's : 1/2, 1/4, 1/8 ...
            icp[i] = icp[i-1]*0.5;
        for(i=1; i<nsplt; i++)                 // now offset to 0.5, 0.75, 0.876 ...
            icp[i] += icp[i-1];

        //splitCell(isplt, nsplt, icp);               
        
        for(i=0; i<nsplt-1; i++)               // update dx arr, reuse icp arr
            icp[i] = icp[i+1] - icp[i];
        if(nsplt >=2) 
            icp[nsplt-1] = icp[nsplt-2];
        dx[isplt] *= 0.5;
        dx.insert(dx.begin()+isplt+1, icp.begin(), icp.end());
    }

    ngrd  += nsplt;                            // update adaptMesh's ngrd, ngrdf (odtl's are done)
    ngrdf += nsplt;

    //--------- Update the mark array

    for(i=iglobal+1; i<(int)mark.size(); i++)
        if(mark[i] > isplt)
            mark[i] += nsplt;

    //--------- Now compare the other half of the split cell with neighbor

    int inext;
    if(splitCellOnRight) {
        inext = isplt + nsplt;
        if(inext == ngrd-1 && !odtP->Lperiodic)
            return;
    }
    else {
        inext = mPos-1;
        if(inext == -1) {
            if(!odtP->Lperiodic)
                return;
            else
                inext = ngrd-1;
        }
    }

    fix2point5offender(inext, iglobal);                 // recursive call 

}


/////////////////////////////////////////////////////////////////////////////// 
/** Split any cells larger than the limit.
*  Not really needed, except maybe once at the very start.  The odt solution and
*  mesh adaptation will not create large cells to be split.  The initial condition
*  should be set up with enough cells to not violate the condition.  That is,
*  the user should be consistent in the specification of domain size, ncells, and
*  large cell size.
*  The user should ensure that the large limit is at least twice the small limit
*/
void adaptMesh::splitLargeCells() {
    vector<double> icp(1);
    for (int i=0; i < ngrd; i++) {
        if (anyl->posf[i+1]-anyl->posf[i] > odtP->dxmax){
            // Cell has to be split
            icp[0] = anyl->pos[i]-anyl->posf[i];
            anyl->splitCell(i, 1, icp, true);
            dx[i] *= 0.5;
            dx.insert(dx.begin()+i,dx[i]);
            ngrd++;
            ngrdf++;
            i--;
        }
    }
    xnf = anyl->posf;
}


///////////////////////////////////////////////////////////////////////////////
/** Resizes and sets the dx array to reflect the correct distances.
 *
 */
void adaptMesh::setDxArray() {

    dx.resize(ngrd);
    for(int i=0, ip=1; i<ngrd; i++, ip++)
        dx[i] = xnf[ip] - xnf[i];

}


////////////////////////////////////////////////////////////////////////////////
/** Given a position, find an index (use with pos)
 * 
 * @param x      \input vector of cell positions
 * @param val    \input location which is converted to index
 * @param istart \input start looking at this index
 */
int adaptMesh::findPos(const vector<double> &x, const double val, int &istart) {

    if(val <= x[0])
        return 0;
    if(val >= x[x.size()-1])
        return x.size()-2;
    for(int i=istart; i<(int)x.size(); i++)
        if(x[i] >= val) {
            return i-1;
        }
    return 0;
}


////////////////////////////////////////////////////////////////////////////////
/** Linear interpolation of a single point, given two vectors \c x and <code>y</code>.
 * 
 * @param x      \input  vector of cell positions
 * @param y      \input  vector of cell values
 * @param xval   \input  position to interpolate
 * @param yval   \output value at interpolated position
 * @param istart \output start looking at this index (for findPos)
 */
void adaptMesh::interp1pt(const vector<double> &x, const vector<double> &y, 
        const double &xval, double &yval, int &istart) {

    int i, ip;

    i = findPos(x, xval, istart);
    ip = i+1;
    if(i>0)
        istart = i-1;

    yval = y[i] + (xval-x[i])*(y[ip]-y[i])/(x[ip]-x[i]);

}


////////////////////////////////////////////////////////////////////////////////
/** Linear interpolation of a vector of points (call interp1pt each time)
 * 
 * @param x      \input  vector of cell positions
 * @param y      \input  vector of cell values
 * @param xn     \input  vector of positions to interp to
 * @param yn     \output vector of interpolated values
 */
void adaptMesh::interpVec(const vector<double> &x, const vector<double> &y,
               const vector<double> &xn, vector<double> &yn) {

    if(x.size() != y.size() || xn.size() != yn.size()) {
        cerr << "\nERROR IN INTERPVEC" << endl;
        exit(0);
    }
    int istart = 0;
    for(int i=0; i<(int)xn.size(); i++) 
        interp1pt(x,y,xn[i],yn[i], istart);
}


////////////////////////////////////////////////////////////////////////////////
/** Computes and returns total distance (length) along a curve specified by vectors x and y
 * 
 * @param x      \input  vector of cell positions
 * @param y      \input  vector of vector cell values for profiles (phi's) to compare
 * @param sDist  \input  vector that stores the running distance along the curve
 *
 * @returns the total distance(length) along the specified curve.
 */
double adaptMesh::calcDistance(const vector<double> &x, const vector<vector<double> > &y,
                               vector<double> &sDist) {

    double dx2, dy2;    // dx^2 and dy^2
    double dmb;

    sDist[0] = 0.0;
    for (int i=1; i<(int)x.size(); i++) {

        dx2 = x[i]-x[i-1];
        dx2 *= dx2;

        dy2 = 0.0; 
        for(int iProf=0; iProf<(int)y.size(); iProf++) {
            dmb = y[iProf][i] - y[iProf][i-1];
            dmb *= dmb;
            if(dmb > dy2) dy2 = dmb;
        }

        sDist[i]= sDist[i-1]+sqrt(dx2 +dy2);
    }   

    return sDist[sDist.size()-1]; 
}


////////////////////////////////////////////////////////////////////////////////






