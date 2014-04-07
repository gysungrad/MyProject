/**
 * @file anyline.cc
 * Source file for class anyline
 */

#include "anyline.h"
#include "processor.h"
#include <iomanip>
#include <cmath>          // fabs
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <algorithm>

using namespace std;

extern processor proc;

///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor function.
 * Class is the parent of odtline, which is usually what one creates.
 * Note, props is not set here, its set by its derived classes.  So, upBprops
 * and loBprops are set there too, as is propNames
 *
 * @param npts \input initial number of cells in the line.  Has default. (Why does this say "input:"?)
 * @param Ld \input initial line length.  Has default.
 */
anyline::anyline(int npts, double Ld) { 


    ngrd    = npts;
    ngrdf   = ngrd+1;
    Ldomain = Ld;
    pos     = vector<double>(ngrd,  0.0);
    posf    = vector<double>(ngrdf, 0.0);
    rho     = vector<double>(ngrd,  1.186);    // kg/m3 air at stp
    molec   = vector<double>(ngrd,  1.86E-5);  // kg/m*s air at stp
    lambda  = vector<double>(ngrd,  2.62E-2);  // W/m*K air at stp
    phase   = vector<double>(ngrd,  0.0);    

    nprops = 0;
    uptr   = 0;
    nspc   = 0;
    odtP   = 0;

    double dp = Ldomain/ngrd;
    for(int i=1; i<ngrdf; i++)
        posf[i] = posf[i-1] + dp;
    pos[0] = dp/2;
    for(int i=1; i<ngrd; i++)
        pos[i] = pos[i-1] + dp;

    upBrho    = 1.E20;
    loBrho    = 0.0;
    upBmolec  = 1.E20;
    loBmolec  = 0.0;
    upBlambda = 1.E20;
    loBlambda = 0.0;
    upBuptr   = 1.0E20;
    loBuptr   = -1.0E20;
 
    LpropsHasYi = false;
    iPropsYi    = -1;

}

///////////////////////////////////////////////////////////////////////////////
/**Find index of cell for given position (residing in cell).
 * Start search assuming a uniform grid,
 * then search forward or back till hit the cell index.
 * If position is on cell face j, then if LowSide true, return j, else j-1.                                                                         \n
 * For start of eddy region, set LowSide to true                                                                                                    \n
 * For end of eddy region, set LowSide to false                                                                                                     \n
 * (This is so triplet maps don't overlap cells)
 *                                                                                                                                                  <pre><code>
 * e.g., usual:   | { | | | | } |    5 pts, eddy pos between cell faces
 *       okay:    {   | | | |   }    5 pts, eddy pos on cell faces (1 or both)
 *       bad:     |   { | | }   |    5 pts, eddy pos on internal faces (1 or both)
 *                                                                                                                                                  </code></pre>
 * @see eddy::tripMap
 * 
 * @param position \input position to find corresponding index.
 * @param LowSide \input flag true, then return j if position is on cell face j, else j-1.
 * @return index of position.
 */

int anyline::linePositionToIndex(double position, bool LowSide) {

    if(fabs(position-posf[0]) < 1.0E-15)
        return 0;
    if(fabs(position-posf[ngrd]) < 1.0E-15)
        return ngrd-1;

    // if(position == posf[0])
    //     return 0;
    // if(position == posf[ngrd])  
    //     return ngrd-1;

    if(position < posf[0])         // for periodic (from eddies only)
        position += Ldomain;
    if(position > posf[ngrd])
        position -= Ldomain;

    if(position < posf[0] || position > posf[ngrd]) {
       cerr << "\n ERROR anyline::linePositionToIndex position < posf[0] or > posf[ngrd] \n"
               "at the seed #---> "<< odtP->seed + proc.myid + 12800 << " and @ processor's id---> " << proc.myid
               <<" Value of position is---> "<<position << " and values of posf[0] and posf[ngrd] are "
               <<posf[0]<< " and "<<posf[ngrd] <<" respectively "<< endl; 
       exit(0);
    }

    int i;
    int ipos = static_cast<int>((position-posf[0])/Ldomain*ngrd);  

    if(posf[ipos+1] > position) {      // case 1: grd skewed more pts on right half
        for(i=ipos+1; i>=0; i--)  {
            if(posf[i] <= position) {
                if(position == posf[i]) {
                    if(LowSide)
                        return i;
                    else
                        return i-1;
                }
                else
                    return i;
            }
        }
    }

    else  {                           // case 2: grd skewed more pts on left half
        for(i=ipos+1; i<=ngrdf; i++) {
            if(posf[i] >= position) {
                if(position == posf[i]) { 
                    if(LowSide)
                        return i;
                    else
                        return i-1;
                }
                else
                    return i-1;
            }
        }
    }
    
    cout << "\n\n******** ERROR IN linePositionToIndex " 
         << position << '\t' << posf[0] << '\t' << posf[ngrd] << '\t' << endl << endl;

    return -1;  
    
}

///////////////////////////////////////////////////////////////////////////////

/** Delete old eddy region of original line, insert new triplet mapped, kerneled eddy in its place.
 *
 *  The eddy was done as follows:
 *                                                                                    \vc{
 *  orig grid:       | | | | | { | | } | | | | |
 *  eddy grid:       | { | | } |
 *  trip map it:     | { |||||||| } |
 *                                                                                    }
 *  When kernels were applied they were given to the whole cell
 *                                                                                    <pre><code>
 *      | { |.
 *                                                                                    </code></pre>
 *  We want to split the original grid from
 *                                                                                    <pre><code>
 *      | { |...  to  | | |...,
 *                                                                                    </code></pre>
 *  Then delete the second portion (only the eddy part),
 *  Then change the eddy from
 *                                                                                    <pre><code>
 *      | { |||... to { |||..., (e.g., just the eddy part)
 *                                                                                    </code></pre>
 *  Then insert the old eddy.
 *  SEE ADDITIONAL COMMENTS BELOW FOR periodic eddies
 *  
 *  @param eddyLine  \input line to work on.
 *  @param iStart    \inout starting index of the eddy.
 *  @param iEnd      \inout ending index of the eddy. (What does "input/output" mean?)
 *  @param leftEdge  \inout left edge of the eddy.
 *  @param rightEdge \inout right edge of the eddy.
 *  @param Lwrap     \input flag indicating whether eddy wraps domain (for periodic).


 *
 */
void anyline::insertEddy(anyline &eddyLine, int &iStart, int &iEnd, 
                         double &leftEdge, double &rightEdge, bool Lwrap) {


  
    bool flag1=0, flag2=0;                     // used below for iStart/iEnd
    vector<double> interPos(1);                // used for cell split


    //---------- adjust eddyLine boundaries from | { |||| } | to { |||| }

    eddyLine.posf[0]              = leftEdge;
    eddyLine.posf[eddyLine.ngrd]  = rightEdge;
    eddyLine.pos[0]               = 0.5*(eddyLine.posf[0]+eddyLine.posf[1]);
    eddyLine.pos[eddyLine.ngrd-1] = 0.5*(eddyLine.posf[eddyLine.ngrd] +
                                         eddyLine.posf[eddyLine.ngrd-1]);

 
    //-----------------------------------------------------------------------------
    // No Wrap: as above

    if(!Lwrap) { 

        //---------- split the cells that the eddy bndry intersects

        if(fabs(posf[iEnd+1] - rightEdge) > 1.0e-15) {
            interPos[0] = rightEdge-posf[iEnd];
            splitCell(iEnd, 1, interPos, false);       // don't interpolate! no intrp when grab 
            flag2 = 1;                                   // eddy --> if interp here then violate conservation when insert (ARK found)
        }
        if(fabs(posf[iStart] - leftEdge) > 1.0e-15) {
            interPos[0] = leftEdge-posf[iStart];
            splitCell(iStart, 1, interPos, false);    // again, don't interp
            flag1 = 1;
            iStart++;
            iEnd++;
        }

        ngrd  = ngrd - (iEnd - iStart + 1) + eddyLine.ngrd;
        ngrdf = ngrd + 1;

        //---------- delete the old eddy region, insert the new

        vector<double>::iterator it1, it2;     

        it1 = pos.erase( pos.begin()+iStart, pos.begin()+iEnd+1 );
        pos.insert(it1, eddyLine.pos.begin(),  eddyLine.pos.end());

        it1 = posf.erase( posf.begin()+iStart, posf.begin()+iEnd+2 );
        posf.insert(it1, eddyLine.posf.begin(),  eddyLine.posf.end());

        it1 = rho.erase( rho.begin()+iStart, rho.begin()+iEnd+1 );
        rho.insert(it1, eddyLine.rho.begin(),  eddyLine.rho.end());

        it1 = molec.erase( molec.begin()+iStart, molec.begin()+iEnd+1 );
        molec.insert(it1, eddyLine.molec.begin(), eddyLine.molec.end()); 

        it1 = lambda.erase( lambda.begin()+iStart, lambda.begin()+iEnd+1 );
        lambda.insert(it1, eddyLine.lambda.begin(), eddyLine.lambda.end()); 

        it1 = phase.erase( phase.begin()+iStart, phase.begin()+iEnd+1 );
        phase.insert(it1, eddyLine.phase.begin(), eddyLine.phase.end()); 

        for(int k=0; k<nprops; k++) {
            it1 = (*props[k]).erase( (*props[k]).begin() + iStart, (*props[k]).begin() + iEnd+1 );
            (*props[k]).insert( it1, (*eddyLine.props[k]).begin(), (*eddyLine.props[k]).end() );
        }

        //---------- update iStart, iEnd

        iEnd = iStart + eddyLine.ngrd - 1;           // end indx of eddy in new mainline
        if(flag1) iStart--;                          // not eddy region indicies, but indicies
        if(flag2) iEnd++;                            //   AFFECTED by the eddy (for later meshing)
        // (remember, we split the cells)

    }

    //-----------------------------------------------------------------------------
    /** WRAP:                                                                                \n
    * Orig grid:
    *                                                                                        \vc{
    *           |   |   |   |   |   |   |   |   |
    *                                                                                        }
    * Eddy regn:
    *                                                                                        <pre><tt>
    *           |   | } |   |   |   |   | { |   |            10 cells  (4 in eddy region)
    *                                                                                        </tt></pre>
    * Eddy:
    *                                                                                        <pre><tt>
    *           { | |   | | }                                5 cells
    *                                                                                        </tt></pre>
    * Split the Grid at
    *                                                                                        <pre><tt>
    *           {, }
    *                                                                                        </tt></pre>
    * Split the eddy at domain edge \c +
    *                                                                                        <pre><tt>
    *           { | | + | | }                                6 cells
    *                                                                                        </tt></pre>
    * Delete old eddy region from the grid and insert the new eddy region.                   \n
    * New grid:
    *                                                                                        <pre><tt>
    * + | |e} |   |   |   |   | {s| | +                      12 cells
    *                                                                                        </tt></pre>
    * iStart is at
    *                                                                                        <pre><tt>
    * {s|     on the right
    *                                                                                        </tt></pre>
    * iEnd   is at
    *                                                                                        <pre><tt>
    * |e}     on the left
    *                                                                                        </tt></pre>
    */
    
    else {  

        //---------- split the cells that the eddy bndry intersects

        vector<double> interPos(1); 
        if(posf[iStart] != leftEdge) {
            interPos[0] = leftEdge-posf[iStart];
            splitCell(iStart, 1, interPos, false);       // again, don't interpolate
            flag1 = 1;
            iStart++;
        }
        if(posf[iEnd+1] != rightEdge-Ldomain) {
            interPos[0] = rightEdge-Ldomain-posf[iEnd];
            splitCell(iEnd, 1, interPos, false);         // again, don't interpolate
            flag2 = 1;
            iStart++;
        }

        //---------- split the eddy at the domain boundary (right boundary)

        int ies = eddyLine.linePositionToIndex(posf[ngrd], true);
        if(eddyLine.posf[ies] != posf[ngrd]) {
            interPos[0] = posf[ngrd] - eddyLine.posf[ies];
            eddyLine.splitCell(ies, 1, interPos, false);
            ies++;    // eddy cell on the right of the split
        }

        //---------- shift the right portion of the eddy back

        for(int i=ies; i<eddyLine.ngrd; i++)
            eddyLine.pos[i] -= Ldomain;
        for(int i=ies; i<eddyLine.ngrdf; i++)
            eddyLine.posf[i] -= Ldomain;
        eddyLine.posf[ies] = 0.0;              // redunant, but make it exactly 0.0

        //---------- delete the old eddy region, insert the new

        ngrd  = ngrd - (ngrd-iStart  +  iEnd+1) + eddyLine.ngrd;
        ngrdf = ngrd+1;

        pos.erase(  pos.begin()+iStart, pos.end() );
        pos.erase(  pos.begin(),        pos.begin()+iEnd+1 );
        pos.insert( pos.end(),   eddyLine.pos.begin(),      eddyLine.pos.begin()+ies );
        pos.insert( pos.begin(), eddyLine.pos.begin()+ies,  eddyLine.pos.end() );

        posf.erase(  posf.begin()+iStart, posf.end() );
        posf.erase(  posf.begin(),        posf.begin()+iEnd+2 );
        posf.insert( posf.end(),   eddyLine.posf.begin(),      eddyLine.posf.begin()+ies+1 );
        posf.insert( posf.begin(), eddyLine.posf.begin()+ies,  eddyLine.posf.end() );

        posf[ngrd] = odtP->domainLength;           // (should equal Ldomain)

        rho.erase(  rho.begin()+iStart, rho.end() );
        rho.erase(  rho.begin(),        rho.begin()+iEnd+1 );
        rho.insert( rho.end(),   eddyLine.rho.begin(),      eddyLine.rho.begin()+ies );
        rho.insert( rho.begin(), eddyLine.rho.begin()+ies,  eddyLine.rho.end() );

        molec.erase(  molec.begin()+iStart, molec.end() );
        molec.erase(  molec.begin(),        molec.begin()+iEnd+1 );
        molec.insert( molec.end(),   eddyLine.molec.begin(),      eddyLine.molec.begin()+ies );
        molec.insert( molec.begin(), eddyLine.molec.begin()+ies,  eddyLine.molec.end() );

        lambda.erase(  lambda.begin()+iStart, lambda.end() );
        lambda.erase(  lambda.begin(),        lambda.begin()+iEnd+1 );
        lambda.insert( lambda.end(),   eddyLine.lambda.begin(),      eddyLine.lambda.begin()+ies );
        lambda.insert( lambda.begin(), eddyLine.lambda.begin()+ies,  eddyLine.lambda.end() );

        phase.erase(  phase.begin()+iStart, phase.end() );
        phase.erase(  phase.begin(),        phase.begin()+iEnd+1 );
        phase.insert( phase.end(),   eddyLine.phase.begin(),      eddyLine.phase.begin()+ies );
        phase.insert( phase.begin(), eddyLine.phase.begin()+ies,  eddyLine.phase.end() );

        for(int k=0; k<nprops; k++) {
            (*props[k]).erase(  (*props[k]).begin()+iStart, (*props[k]).end() );
            (*props[k]).erase(  (*props[k]).begin(),        (*props[k]).begin()+iEnd+1 );
            (*props[k]).insert( (*props[k]).end(),   (*eddyLine.props[k]).begin(),      (*eddyLine.props[k]).begin()+ies );
            (*props[k]).insert( (*props[k]).begin(), (*eddyLine.props[k]).begin()+ies,  (*eddyLine.props[k]).end() );
        }

        //---------- update iStart, iEnd

        iEnd   = eddyLine.ngrd - ies - 1;
        iStart = ngrd - (eddyLine.ngrd - ies) - 1;

        if(flag1==1) iStart--;                 // split cell --> set AFFECTED eddy region back a cell
        if(flag2==1) iEnd++;                   // split cell --> set AFFECTED eddy region forward a cell
        if(iStart < iEnd) {iEnd--; iStart++; } // undo it if your eddy is so big you cross

    }


}


///////////////////////////////////////////////////////////////////////////////
/** This funktion erases one cell of the line 
 *
 *  @param  ieras   erases cell ieras.
 */
void anyline::eraseCell(int ieras){
    
    pos.erase(  pos.begin() +ieras);
    posf.erase( posf.begin()+ieras);
    rho.erase(  rho.begin() +ieras);
    molec.erase( molec.begin()+ieras);
    lambda.erase( lambda.begin()+ieras);
    phase.erase( phase.begin()+ieras);
    for(int k=0; k<nprops; k++) 
        (*props[k]).erase( (*props[k]).begin()+ieras );
    
    ngrdf--;
    ngrd--;

}


///////////////////////////////////////////////////////////////////////////////
/**Merge two cells conservatively, preserving mass and momentum, etc.
 * Constant profiles are assumed in each cell when merging
 * imrg is the lower index of the two cells to merge
 *
 * @param imrg \input merge cells imrg and imrg+1.
 */
void anyline::merge2cells(int imrg) {


    int k;
    int ip  = imrg+1;
    int ipp = ip+1;

    double rho3;
    double rhodx1, rhodx2, rhodx3;
    double rhoU3, u3;
    double rhoUdx1, rhoUdx2, rhoUdx3;

    vector<double> props3(nprops);
    double molec3, lambda3;

    double dx1 = posf[ip] - posf[imrg];
    double dx2 = posf[ipp] - posf[ip];
    double dx3 = dx1 + dx2;
    
    if (phase.at(imrg) != phase.at(ip)){
        *proc.ostrm << "\n\nip   = " << ip;
        *proc.ostrm <<   "\nimrg = " << imrg;
        *proc.ostrm << endl << endl << "WARNING:" << endl
                    << "Two cells are marked to be merged but have different phases.";
        this->outputProperties("LineBeforeCrash.dat");
        exit(0);
        return;
    }
    
    //---------- evaluate new cell contents
    if(!odtP->Lspatial) {
        rhodx1 = rho[imrg]*dx1;
        rhodx2 = rho[ip]  *dx2;
        rho3   = (rhodx1 + rhodx2)/dx3; 
        rhodx3 = rho3*dx3;

        for(k=0; k<nprops; k++) {
             props3[k] = (rhodx1*(*props[k])[imrg]
                         + rhodx2*(*props[k])[ip])/rhodx3;
         }

    } 

    else {
        rhoUdx1 = rho[imrg]*(*uptr)[imrg]*dx1;
        rhoUdx2 = rho[ip]*(*uptr)[ip]*dx2;
        rhoU3   = (rhoUdx1 + rhoUdx2)/dx3;
        rhoUdx3 = rhoU3*dx3;
        for(k=0; k<nprops; k++) {
             props3[k] = (rhoUdx1*(*props[k])[imrg]
                         + rhoUdx2*(*props[k])[ip])/rhoUdx3;
         }
        
        u3    = (rhoUdx1*(*uptr)[imrg]+rhoUdx2*(*uptr)[ip])/rhoUdx3;
        rho3  = rhoUdx3/(u3*dx3);
    }

    molec3  = dx3/(dx1/molec[imrg] + dx2/molec[ip]); 
    lambda3 = dx3/(dx1/lambda[imrg] + dx2/lambda[ip]);

    //---------- fill new cell contents

    pos[imrg]    = 0.5*(posf[imrg] + posf[ipp]);
    // don't need to do posf, its just gone
    rho[imrg]    = rho3; 
    molec[imrg]  = molec3;
    lambda[imrg] = lambda3;
    //no change in phase value on merge 
    //(so don't merge cells with different phase values)
    for(k=0; k<nprops; k++) 
        (*props[k])[imrg] = props3[k];

    //---------- delete cell

    pos.erase(  pos.begin() +ip);
    posf.erase( posf.begin()+ip);
    rho.erase(  rho.begin() +ip);
    molec.erase( molec.begin()+ip);
    lambda.erase(lambda.begin()+ip);
    phase.erase( phase.begin()+ip);
    for(k=0; k<nprops; k++) 
        (*props[k]).erase( (*props[k]).begin()+ip );

    ngrdf--;
    ngrd--;

}

///////////////////////////////////////////////////////////////////////////////

/**Merge several cells conservatively.
 * Constant profiles are assumed in each cell when merging.
 * ifirst is the start cell and ncells is how many cells to do.
 *
 * @param ifirst \input merge cell index ifirst with next cells.
 * @param ncells \input total cells to merge.
 */
void anyline::mergeSeveralCells(int ifirst, int ncells) {
    int            i, j, k; 
    int            ip     = ifirst+1;
    int            ilast  = ifirst + ncells - 1; 
    int            ilastp = ilast+1;
    double         sumdx  = 0.0;
    double         rhodx  = 0.0;
    double         rhoUdx = 0.0;
    double         molecC = 0.0;
    double         lambdaC= 0.0;
    double         u3     = 0.0;

    vector<double> rhodxP(nprops,0.0);
    vector<double> rhoUdxP(nprops,0.0);
    vector<double> dx(ncells);
    vector<double> rhodxA(ncells);
    vector<double> rhoUdxA(ncells);

    for(i=0, j=ifirst; i<ncells; i++, j++) {
        dx[i] = posf[j+1] - posf[j];
        sumdx += dx[i];
    }


    if(!odtP->Lspatial) {
        for(i=0, j=ifirst; i<ncells; i++, j++) 
           rhodxA[i] = rho[j]*dx[i];
        for(i=ifirst, j=0; i<=ilast; i++, j++){
           rhodx   += rhodxA[j];
           molecC  += dx[j]/molec[i];
           lambdaC += dx[j]/lambda[i];
           for(k=0; k<nprops; k++)
                rhodxP[k] += rhodxA[j]*(*props[k])[i];
        }
           
       pos[ifirst]  = 0.5*(posf[ifirst]+posf[ilastp]);
       // don't need to change posf
       rho[ifirst]  = rhodx/sumdx;
       molec[ifirst] = sumdx/molecC;
       lambda[ifirst]= sumdx/lambdaC;
       //no change in phase value on merge 
       //(so don't merge cells with different phase values)
       for(int k=0; k<nprops; k++)
            (*props[k])[ifirst] = rhodxP[k]/rhodx;

    }

    else {
        //make u*rho conserved instead of just rho
        for(i=0, j=ifirst; i<ncells; i++, j++) {
             rhoUdxA[i] = rho[j]*(*uptr)[j]*dx[i];
             u3 += rhoUdxA[i]*(*uptr)[j]; //normalized bellow
        }
   
        for(i=ifirst, j=0; i<=ilast; i++, j++) {
             rhoUdx  += rhoUdxA[j];
             molecC  += dx[j]/molec[i];
             lambdaC += dx[j]/lambda[i];
             for(k=0; k<nprops; k++)
                rhoUdxP[k] += rhoUdxA[j]*(*props[k])[i];
        }
             
        pos[ifirst]  = 0.5*(posf[ifirst]+posf[ilastp]);
        // don't need to change posf
        u3 /= rhoUdx;
        rho[ifirst]  = rhoUdx/(u3*sumdx);
        molec[ifirst] = sumdx/molecC;
        lambda[ifirst]= sumdx/lambdaC;
        //no change in phase value on merge 
        //(so don't merge cells with different phase values)
        for(int k=0; k<nprops; k++)
            (*props[k])[ifirst] = rhoUdxP[k]/rhoUdx;
    }

    pos.erase(   pos.begin() +ip,  pos.begin() +ilastp);
    posf.erase(  posf.begin()+ip,  posf.begin()+ilastp);
    rho.erase(   rho.begin() +ip,  rho.begin() +ilastp);
    molec.erase( molec.begin()+ip, molec.begin()+ilastp);
    lambda.erase(lambda.begin()+ip,lambda.begin()+ilastp);
    phase.erase( phase.begin()+ip, phase.begin()+ilastp);
    for(int k=0; k<nprops; k++) 
        (*props[k]).erase( (*props[k]).begin()+ip, (*props[k]).begin()+ilastp );

    ngrdf -= ncells-1;
    ngrd  -= ncells-1;

}

///////////////////////////////////////////////////////////////////////////////

/**Output the anyline properties.  
 *
 * @param fname \input file name to write.
 */
void anyline::outputProperties(std::string fname) {

   ofstream ofile(fname.c_str()); 
   if(!ofile) 
       cout << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
   
   ofile << "# grid points = " << ngrd;
   ofile << "\n# Domain Size = " << Ldomain;
   if(fabs(posf[ngrd] -posf[0] - Ldomain) > 1.0E-6)
       ofile << "\n# last posf-first posf != Ldomain, last posf = " << posf[ngrd];
   ofile << endl;
   ofile << "#  1_pos            "
         << "2_posf            "
         << "3_rho             "
         << "4_molec           "
         << "5_lambda          "
         << "6_phase           ";
   for(int k=0; k<nprops; k++) 
       ofile << k+7 << "_" << propNames[k] << "            ";

   ofile << scientific;
   ofile << setprecision(10);
   for(int i=0; i<ngrd; i++) {
       ofile << endl;
       ofile << setw(19) << pos[i] 
             << setw(19) << posf[i]
             << setw(19) << rho[i]
             << setw(19) << molec[i]
             << setw(19) << lambda[i]
             << setw(19) << phase[i];
       for(int k=0; k<nprops; k++) 
           ofile << setw(19) << (*props[k])[i];
   }

   ofile.close();
}

///////////////////////////////////////////////////////////////////////////////

/**Read anyline properties from an anyline output file.
 * Members like propNames, and the bounds variables are not read
 *
 * @param fname \input file name to write.
 */
void anyline::readProperties(std::string fname) {

   *proc.ostrm << endl << "Reading anyline property file " << fname << endl;

   ifstream ifile(fname.c_str()); 
   if(!ifile) {
       cout << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
       exit(0);
   }
   
   string       s1;
   stringstream ss1;

   getline(ifile, s1);                            // read line "# grid point = 100"
   ss1.str(s1);
   ss1 >> s1 >> s1 >> s1 >> s1 >> ngrd;
   ss1.clear();
   getline(ifile, s1);                            // read line "# Domain Size = 2"
   ss1.str(s1);
   ss1 >> s1 >> s1 >> s1 >> s1 >> Ldomain;
   getline(ifile, s1);

   ngrdf = ngrd + 1;

   pos.resize(ngrd);
   posf.resize(ngrdf);
   rho.resize(ngrd);
   molec.resize(ngrd);
   lambda.resize(ngrd);
   phase.resize(ngrd);

   for(int k=0; k<nprops; k++) 
     (*props[k]).resize(ngrd);

   for(int i=0; i<ngrd; i++) {
       ifile >> pos[i]
             >> posf[i]
             >> rho[i] 
             >> molec[i]
             >> lambda[i]
             >> phase[i];
       for(int k=0; k<nprops; k++) 
           ifile >> (*props[k])[i];
   }

   posf[ngrd] = posf[0] + Ldomain;

   ifile.close();

}


///////////////////////////////////////////////////////////////////////////////

/**  Split cell isplt into nsplt+1 pieces (cell is split nsplt times).
 *
 * @param isplt          \input cell index to split.
 * @param nsplt          \input number of cells to split.
 * @param interCellPosns \input positions of added faces.
 * @param Linterp        \input true to interpolate split variables with conservation.
 *
 * If Linterp is false, all contents are the same except the positions, 
 * which are in interCellPosns.
 * IntercellPosns are relative to left edge of cell
 * If Linterp is true, then conservatively interpolate the split values
 * to be between the original neighbors of isplt.
 * A parabola (or a line, or a constant) is fit and used to interpolate the variables.
 * Conditions:  if isplt is a min/max in the given variable, then don't
 * allow the variable to go out of bounds.
 * For example: \vc{| 1 |        2        | 1 |    ---->  | 1 | 1.5   | 2.5 |   1.5    | 1 |}
 * That is, we split the middle cell twice, then its profile goes up in the middle and
 * down on the sides to be conservative.  This is allowed, but stay in the realizable 
 * bounds.  If the original profile is monotonic, then it has to stay monotonic.
 * For profiles like phi = 1 1 2.  To split the third cell, have to retain the tophat.
 * For profiles like 1 2 10, fitting a parabola can result in a new min or max in the
 * split cell.  This is not allowed.
 * 
 * For species with interpolation, rho is interpolated independent of rho*Yi,
 * but then sum of rho*Yi is not equal to rho after the interpolation.
 * Instead, if LpropsHasYi is true, then compute rho after interpolation from
 * rho*Yi
 * This flag assumes that the props array is enthalpy, the Yi for all species
 * @sideeffect changes rho, molec, phase, pos, posf vector sizes\n
 *             changes ngrd, ngrdf values
 */
void anyline::splitCell(int isplt, int nsplt, std::vector<double> &interCellPosns, bool Linterp) {

    //if(odtP->Lspatial) Linterp = false; //eim no interpolation

    //-------------------------------------------------------------------------

    if(!Linterp) {

      
        rho.insert( rho.begin() +isplt, nsplt, rho[isplt]);
        molec.insert(molec.begin()+isplt, nsplt, molec[isplt]);
        lambda.insert(lambda.begin()+isplt, nsplt, lambda[isplt]);
        phase.insert(phase.begin()+isplt, nsplt, phase[isplt]);

        for(int k=0; k<nprops; k++) 
            (*props[k]).insert( (*props[k]).begin()+isplt, nsplt, (*props[k])[isplt] );

        pos.insert( pos.begin() +isplt, nsplt, 0.0);
        posf.insert(posf.begin()+isplt, nsplt, 0.0);

        posf[isplt] = posf[isplt+nsplt];

        int i,j;
        for(i=isplt+1, j=0; i<=isplt+nsplt; i++, j++)
            posf[i] = interCellPosns[j]+posf[isplt];
        for(i=isplt; i<=isplt+nsplt; i++)
            pos[i] = 0.5*(posf[i] + posf[i+1]);

        ngrd  += nsplt;
        ngrdf += nsplt;

        return;
    }

    ///------------------------------------------------------------------------
    ///
    /// For Linterp = true, the spatial and temporal formulations interpolation
    /// is done seperately in their own functions
    
    if(!odtP->Lspatial)
        splitCellTemporal(isplt,nsplt,interCellPosns);
    else
        splitCellSpatial(isplt,nsplt,interCellPosns);

    
    //----------------------------------------------------------------------------------

    //----------------- do the positions 

    pos.insert( pos.begin() +isplt, nsplt, 0.0);
    posf.insert(posf.begin()+isplt, nsplt, 0.0);

    posf[isplt] = posf[isplt+nsplt];

    int i,j;
    for(i=isplt+1, j=0; i<=isplt+nsplt; i++, j++)
        posf[i] = interCellPosns[j]+posf[isplt];
    for(i=isplt; i<=isplt+nsplt; i++)
        pos[i] = 0.5*(posf[i] + posf[i+1]);

    ngrd  += nsplt;
    ngrdf += nsplt;

}

//////////////////////////////////////////////////////////////////////////////////////////
/**  Called by anyline::splitCell. Interpolation used to split cells in temporal cases
 *  @param isplt \input index of cell to split
 *  @param nsplt \input number of times to split the cell
 *  @param interCellPosns \input locations from left edge to split
*/

void anyline::splitCellTemporal(int isplt, int nsplt, std::vector<double> &interCellPosns){
                        
    
    //-------------------------------------------------------------------------
    /**
    * Fit a parabola through the cell to be split and its
    * two neighbors.  \fun{y = a*x^2 + b*x + c}.
    * Three degrees of freedom: a, b, c.  Three constraints: pass through the
    * neighboring y's and the integral of y in the cell to be split (s) is
    * \fun{y_s*dx}.  Then, once a,b,c are known, split the cell at positions interCellPosns
    * and compute the y in those cells by equating the integral of the parabola
    * in that smaller region to \fun{y*dx} of the region to compute y.
    * Integration is for the conserved variables, \fun{\rho}, \fun{\rho*u}, \fun{\rho*v},
    * \fun{\rho*w}, \fun{\rho*Yi} etc.
    * If the cell of interest has the same value as its two neighbors, then do
    * as for Linterp false.
    * If the cell to split has the same value as one of its neighbors, do a line
    * instead of a parabola.  Likewise if the cell is on a boundary.
    * Note, molec is just interpolated as if a regular variable, no harmonic mean
    * buisiness, unlike for the merging.
    */
    
    double         a, b, c;
    //double         dx = posf[isplt+1] - posf[isplt]; //  !!!!!  unused variable
    int            ip = isplt+1;
    int            im = isplt-1;
    vector<double> var(nsplt+1,0.0);
    double         y1, y2, y3;
    int            iCheck;

    //-------------------------------------------------------------------------


    if(ip == ngrd)
        ip = 0;
    if(im == -1)
        im = ngrd-1;

    // mass being conserved
    double oldRho1 = rho[im];
    double oldRho2 = rho[isplt];
    double oldRho3 = rho[ip];

    //---------- molec

    y1 = molec[im];
    y2 = molec[isplt];
    y3 = molec[ip];

    iCheck = getABC(a,b,c, isplt, y1, y2, y3); 
    interpTheVar(var, a,b,c, y1,y2,y3, isplt, iCheck, interCellPosns);
    if(iCheck==1) {                         // don't interp if over/undershoot
        double maxV = *max_element(var.begin(), var.end());
        double minV = *min_element(var.begin(), var.end());
        if(maxV > upBmolec || minV < loBmolec) 
            for(int i=0; i<(int)var.size(); i++)
                var[i] = y2;
    }
    molec.erase(molec.begin()+isplt);
    molec.insert(molec.begin()+isplt, var.begin(), var.end());

    //---------- lambda

    y1 = lambda[im];
    y2 = lambda[isplt];
    y3 = lambda[ip];

    iCheck = getABC(a,b,c, isplt, y1, y2, y3); 
    interpTheVar(var, a,b,c, y1,y2,y3, isplt, iCheck, interCellPosns);
    if(iCheck==1) {                         // don't interp if over/undershoot
        double maxV = *max_element(var.begin(), var.end());
        double minV = *min_element(var.begin(), var.end());
        if(maxV > upBlambda || minV < loBlambda) 
            for(int i=0; i<(int)var.size(); i++)
                var[i] = y2;
    }
    lambda.erase(lambda.begin()+isplt);
    lambda.insert(lambda.begin()+isplt, var.begin(), var.end());

    //---------- phase

    phase.insert(phase.begin()+isplt, nsplt, phase[isplt]);

    ///----------------------------------------------------------------------------------
    /**
    * Different treatment for odtline and scaline here. \n
    * odtline is here, scaline is below                 \n
    * scaline implies props is enth, Yspecies           \n
    */

    if(!LpropsHasYi) {  

        //---------- rho

        y1 = rho[im];
        y2 = rho[isplt];
        y3 = rho[ip];

        iCheck = getABC(a,b,c, isplt, y1, y2, y3);
        interpTheVar(var, a,b,c, y1,y2,y3, isplt, iCheck, interCellPosns);
        if(iCheck==1) {                         // don't interp if over/undershoot
            double maxV = *max_element(var.begin(), var.end());
            double minV = *min_element(var.begin(), var.end());
            if(maxV > upBrho || minV < loBrho)
                for(int i=0; i<(int)var.size(); i++)
                    var[i] = y2;
        }
        rho.erase(rho.begin()+isplt);
        rho.insert(rho.begin()+isplt, var.begin(), var.end());

        //---------- props

        for(int k=0; k<nprops; k++) {

            y2 = (*props[k])[isplt]*oldRho2;
            if(isplt > 0)
                y1 = (*props[k])[im]*oldRho1;
            if(isplt==0 && odtP->Lperiodic) 
                y1 = ((*props[k])[im]-odtP->pJump[k])*oldRho1;
            if(isplt < ngrd-1)
                y3 = (*props[k])[ip]*oldRho3;
            if(isplt==ngrd-1 && odtP->Lperiodic) 
                y3 = ((*props[k])[ip]+odtP->pJump[k])*oldRho3;

            iCheck = getABC(a,b,c, isplt, y1, y2, y3);
            interpTheVar(var, a,b,c, y1,y2,y3, isplt, iCheck, interCellPosns);
            for(int i=0; i<(int)var.size(); i++)
                var[i] /= rho[isplt+i];
            if(iCheck==1) {                         // don't interp if over/undershoot
                double maxV = *max_element(var.begin(), var.end());
                double minV = *min_element(var.begin(), var.end());
                if(maxV > upBprops[k] || minV < loBprops[k]) {
                    y2 /= oldRho2; 
                    for(int i=0; i<(int)var.size(); i++)
                        var[i] = y2;
                } 
            } 
            (*props[k]).erase((*props[k]).begin()+isplt);
            (*props[k]).insert((*props[k]).begin()+isplt, var.begin(), var.end());
        }
    }
    
    /**----------------------------------------------------------------------------------
    *
    * Scaline here.
    * For interpolation, each variable \fun{\rho*props} is done independently, so
    * density will not be consistent with \fun{\sum \rho*Y}, so get density from this
    * sum rather than by interpolating density.  Things are a bit reordered
    * because of this.  It may be worth reworking this (but only for clarity)
    */

    else {

        vector<double>          newRho(nsplt+1, 0.0);
        vector<vector<double> > varYs(nspc, vector<double>(nsplt+1));

        //---------- props (species)

        for(int k=iPropsYi; k<nprops; k++) {        // note k=1, not 0 (0 is enth, done below)

            y2 = (*props[k])[isplt]*oldRho2;
            if(isplt > 0)
                y1 = (*props[k])[im]*oldRho1;
            if(isplt==0 && odtP->Lperiodic) 
                y1 = ((*props[k])[im]-odtP->pJump[k])*oldRho1;
            if(isplt < ngrd-1)
                y3 = (*props[k])[ip]*oldRho3;
            if(isplt==ngrd-1 && odtP->Lperiodic) 
                y3 = ((*props[k])[ip]+odtP->pJump[k])*oldRho3;

            iCheck = getABC(a,b,c, isplt, y1, y2, y3);
            interpTheVar(var, a,b,c, y1,y2,y3, isplt, iCheck, interCellPosns);
            if(iCheck==1) {                         // don't interp if over/undershoot
                double maxV = *max_element(var.begin(), var.end());
                double minV = *min_element(var.begin(), var.end());
                // note, here bounds are on rho*Y and go from 0 to inf
                if(maxV > upBprops[k] || minV < loBprops[k]) 
                    for(int i=0; i<(int)var.size(); i++)
                        var[i] = y2;
            }
            
            // computing density from sum(rho*Y) rather than interpolate rho
            // can't interpolate rho else its not consistent with interpolated rhoY
            for(int i=0; i<(int)newRho.size(); i++)
                newRho[i] += var[i];

            // store the interpolated rhoYi till rho is computed, then do the Y's
            varYs[k-iPropsYi] = var;
            
        }
        // now that we have the newRho, get the Y's
        for(int k=0; k<nspc; k++)
            for(int i=0; i<(int)var.size(); i++)
                varYs[k][i] /= newRho[i];

        // now reinsert the arrays
        for(int k=iPropsYi; k<nprops; k++) {
            (*props[k]).erase((*props[k]).begin()+isplt);
            (*props[k]).insert((*props[k]).begin()+isplt, varYs[k-iPropsYi].begin(), varYs[k-iPropsYi].end());
        }

        //---------- rho

        rho.erase(rho.begin()+isplt);
        rho.insert(rho.begin()+isplt, newRho.begin(), newRho.end());

        //---------- props (e.g., enthalpy)
        // enthalpy need rho which is now done, so get em here

        for(int k=0; k<iPropsYi; k++) {
        
            y2 = (*props[k])[isplt]*oldRho2;
            if(isplt > 0)
                y1 = (*props[k])[im]*oldRho1;
            if(isplt==0 && odtP->Lperiodic) 
                y1 = ((*props[k])[im]-odtP->pJump[k])*oldRho1;
            if(isplt < ngrd-1)
                y3 = (*props[k])[ip]*oldRho3;
            if(isplt==ngrd-1 && odtP->Lperiodic) 
                y3 = ((*props[k])[ip]+odtP->pJump[k])*oldRho3;

            iCheck = getABC(a,b,c, isplt, y1, y2, y3);
            interpTheVar(var, a,b,c, y1,y2,y3, isplt, iCheck, interCellPosns);
            for(int i=0; i<(int)var.size(); i++)
                var[i] /= rho[isplt+i];
            if(iCheck==1) {                         // don't interp if over/undershoot
                double maxV = *max_element(var.begin(), var.end());
                double minV = *min_element(var.begin(), var.end());
                if(maxV > upBprops[k] || minV < loBprops[k]) {
                    y2 /= oldRho2; 
                    for(int i=0; i<(int)var.size(); i++)
                        var[i] = y2;
                }   
            }   
            (*props[k]).erase((*props[k]).begin()+isplt);
            (*props[k]).insert((*props[k]).begin()+isplt, var.begin(), var.end());
        }
    }
    return;
    
}

//////////////////////////////////////////////////////////////////////////////////////////

/** Called by #splitCell. Interpolation used to split cells in spatial cases
 *  @param isplt \input index of cell to split
 *  @param nsplt \input number of times to split the cell
 *  @param interCellPosns \input locations from left edge to split
*/

void anyline::splitCellSpatial(int isplt, int nsplt, std::vector<double> &interCellPosns) {
                        

    if(odtP->Lperiodic) {
        cout << "\n***************** ERROR Periodic does not work in spatial " << endl;
        exit(0);
    }


    //-------------------------------------------------------------------------
    /**
    * Fit a parabola through the cell to be split and its
    * two neighbors.  \fun{y = a*x^2 + b*x + c}.
    * Three degrees of freedom: a, b, c.  Three constraints: pass through the
    * neighboring y's and the integral of y in the cell to be split (s) is
    * \fun{y_s*dx}.  Then, once a,b,c are known, split the cell at positions interCellPosns
    * and compute the y in those cells by equating the integral of the parabola
    * in that smaller region to \fun{y*dx} of the region to compute y.
    * Integration is for the conserved variables, \fun{\rho*u} \fun{\rho*u*u}, \fun{\rho*u*v}, \fun{\rho*u*Y_i} etc.
    * If the cell of interest has the same value all its two neighbors, then do
    * as for Linterp false.
    * If the cell to split has the same value as one of its neighbors, do a line
    * instead of a parabola.  Likewise if the cell is on a boundary.
    * Note, molec is just interpolated as if a regular variable, no harmonic mean
    * buisiness, unlike for the merging.
    */

    double         a, b, c;
    //double         dx = posf[isplt+1] - posf[isplt]; //  !!!!!  unused variable
    int            ip = isplt+1;
    int            im = isplt-1;
    vector<double> var(nsplt+1,0.0);
    vector<double> rhoU(nsplt+1,0.0);
    vector<double> newRho(nsplt+1,0.0);
    double         y1, y2, y3;
    int            iCheck;

    //-------------------------------------------------------------------------


    if(ip == ngrd)
        ip = 0;
    if(im == -1)
        im = ngrd-1;
    
    // mass flux being conserved rather then mass
    double oldRhoU1 = rho[im]*(*uptr)[im];
    double oldRhoU2 = rho[isplt]*(*uptr)[isplt];
    double oldRhoU3 = rho[ip]*(*uptr)[ip];

    //---------- molec

    y1 = molec[im];
    y2 = molec[isplt];
    y3 = molec[ip];

    iCheck = getABC(a,b,c, isplt, y1, y2, y3); 
    interpTheVar(var, a,b,c, y1,y2,y3, isplt, iCheck, interCellPosns);
    if(iCheck==1) {                         // don't interp if over/undershoot
        double maxV = *max_element(var.begin(), var.end());
        double minV = *min_element(var.begin(), var.end());
        if(maxV > upBmolec || minV < loBmolec) 
            for(int i=0; i<(int)var.size(); i++)
                var[i] = y2;
    }
    molec.erase(molec.begin()+isplt);
    molec.insert(molec.begin()+isplt, var.begin(), var.end());

    //---------- lambda

    y1 = lambda[im];
    y2 = lambda[isplt];
    y3 = lambda[ip];

    iCheck = getABC(a,b,c, isplt, y1, y2, y3); 
    interpTheVar(var, a,b,c, y1,y2,y3, isplt, iCheck, interCellPosns);
    if(iCheck==1) {                         // don't interp if over/undershoot
        double maxV = *max_element(var.begin(), var.end());
        double minV = *min_element(var.begin(), var.end());
        if(maxV > upBlambda || minV < loBlambda) 
            for(int i=0; i<(int)var.size(); i++)
                var[i] = y2;
    }
    lambda.erase(lambda.begin()+isplt);
    lambda.insert(lambda.begin()+isplt, var.begin(), var.end());

    //---------- phase

    phase.insert(phase.begin()+isplt, nsplt, phase[isplt]);

    
    /**----------------------------------------------------------------------------------
    *
    * Different treatment for odtline and scaline here
    * odtline is here, scaline is below
    * scaline implies props is enth, Yspecies
    */

    if(!LpropsHasYi) {  

        //---------- rho*uvel(mass flux)

        y1 = rho[im]*(*uptr)[im];
        y2 = rho[isplt]*(*uptr)[isplt];
        y3 = rho[ip]*(*uptr)[ip];

        iCheck = getABC(a,b,c, isplt, y1, y2, y3);
        interpTheVar(var, a,b,c, y1,y2,y3, isplt, iCheck, interCellPosns);
        
        rhoU = var;

        if(iCheck==1) {                         // don't interp if over/undershoot
            double maxV = *max_element(var.begin(), var.end());
            double minV = *min_element(var.begin(), var.end());
            if(maxV > upBuptr || minV < loBuptr)
                for(int i=0; i<(int)var.size(); i++)
                    var[i] = y2;
        }
        //-------uvel to find rho

        // does not do periodic
        y1 = oldRhoU1*(*uptr)[im];
        y2 = oldRhoU2*(*uptr)[isplt];
        y3 = oldRhoU3*(*uptr)[ip];

        iCheck = getABC(a,b,c, isplt, y1, y2, y3);
        interpTheVar(var, a,b,c, y1,y2,y3, isplt, iCheck, interCellPosns);
        for(int i=0; i<(int)var.size(); i++)
            var[i] /= rhoU[i];
       
        if(iCheck==1) {                         // don't interp if over/undershoot
            double maxV = *max_element(var.begin(), var.end());
            double minV = *min_element(var.begin(), var.end());
            if(maxV > upBuptr || minV < loBuptr) {
                y2 /= oldRhoU2; 
                for(int i=0; i<(int)var.size(); i++)
                    var[i] = y2;
            } 
        }
        
        //---------- rho
        
        for(int i=0; i<(int)var.size(); i++)
            newRho[i] = rhoU[i]/var[i];


        rho.erase(rho.begin()+isplt);
        rho.insert(rho.begin()+isplt, newRho.begin(), newRho.end());

        //---------- props

        for(int k=0; k<nprops; k++) {

            y2 = (*props[k])[isplt]*oldRhoU2;
            if(isplt > 0)
                y1 = (*props[k])[im]*oldRhoU1;
            if(isplt==0 && odtP->Lperiodic) 
                y1 = ((*props[k])[im]-odtP->pJump[k])*oldRhoU1;
            if(isplt < ngrd-1)
                y3 = (*props[k])[ip]*oldRhoU3;
            if(isplt==ngrd-1 && odtP->Lperiodic) 
                y3 = ((*props[k])[ip]+odtP->pJump[k])*oldRhoU3;

            iCheck = getABC(a,b,c, isplt, y1, y2, y3);
            interpTheVar(var, a,b,c, y1,y2,y3, isplt, iCheck, interCellPosns);
            for(int i=0; i<(int)var.size(); i++)
                var[i] /= rhoU[i];
            if(iCheck==1) {                         // don't interp if over/undershoot
                double maxV = *max_element(var.begin(), var.end());
                double minV = *min_element(var.begin(), var.end());
                if(maxV > upBprops[k] || minV < loBprops[k]) {
                    y2 /= oldRhoU2; 
                    for(int i=0; i<(int)var.size(); i++)
                        var[i] = y2;
                } 
            } 
            (*props[k]).erase((*props[k]).begin()+isplt);
            (*props[k]).insert((*props[k]).begin()+isplt, var.begin(), var.end());
        }
    }
    
    /**----------------------------------------------------------------------------------
    *
    * Scaline here.
    * For interpolation, each variable \fun{\rho_U*props} is done independently, so
    * \fun{\rho_U} will not be consistent with \fun{\sum \rho_U*Y}, so get \fun{\rho_U} from this
    * sum rather than by interpolating \fun{\rho_U}.  Things are a bit reordered
    * because of this.  It may be worth reworking this (but only for clarity)
    */

    else {

        vector<vector<double> > varYs(nspc, vector<double>(nsplt+1));

        //---------- props (species)

        for(int k=iPropsYi; k<nprops; k++) {        // note k=1, not 0 (0 is enth, done below)

            y2 = (*props[k])[isplt]*oldRhoU2;
            if(isplt > 0)
                y1 = (*props[k])[im]*oldRhoU1;
            if(isplt==0 && odtP->Lperiodic) 
                y1 = ((*props[k])[im]-odtP->pJump[k])*oldRhoU1;
            if(isplt < ngrd-1)
                y3 = (*props[k])[ip]*oldRhoU3;
            if(isplt==ngrd-1 && odtP->Lperiodic) 
                y3 = ((*props[k])[ip]+odtP->pJump[k])*oldRhoU3;

            iCheck = getABC(a,b,c, isplt, y1, y2, y3);
            interpTheVar(var, a,b,c, y1,y2,y3, isplt, iCheck, interCellPosns);
            if(iCheck==1) {                         // don't interp if over/undershoot
                double maxV = *max_element(var.begin(), var.end());
                double minV = *min_element(var.begin(), var.end());
                
                // note, here bounds are on rho*uvel*Y and go from inf to inf,
                // because of uvel
               // if(maxV > upBuptr || minV < loBuptr) eimdb
                if(maxV > upBprops[k] || minV < loBprops[k]) //eimdb
                    for(int i=0; i<(int)var.size(); i++)
                        var[i] = y2;
            }
            
            // computing rho*uvel from sum(rho*uvel*Y) rather than interpolate rho
            // can't interpolate rho*uvel else its not consistent with interpolated rhoY
            for(int i=0; i<(int)rhoU.size(); i++)
                rhoU[i] += var[i];

            // store the interpolated rhoUYi till rho is computed, then do the Y's
            varYs[k-iPropsYi] = var;
            
        }
        // now that we have the rhoU, get the Y's
        for(int k=0; k<nprops-iPropsYi; k++)
            for(int i=0; i<(int)var.size(); i++)
                varYs[k][i] /= rhoU[i];

        // now reinsert the arrays
        for(int k=iPropsYi; k<nprops; k++) {
            (*props[k]).erase((*props[k]).begin()+isplt);
            (*props[k]).insert((*props[k]).begin()+isplt, varYs[k-iPropsYi].begin(), varYs[k-iPropsYi].end());
        }

        //-------uvel to find rho

        // does not do periodic
        y1 = oldRhoU1*(*uptr)[im];
        y2 = oldRhoU2*(*uptr)[isplt];
        y3 = oldRhoU3*(*uptr)[ip];

        iCheck = getABC(a,b,c, isplt, y1, y2, y3);
        interpTheVar(var, a,b,c, y1,y2,y3, isplt, iCheck, interCellPosns);
        for(int i=0; i<(int)var.size(); i++)
            var[i] /= rhoU[i];
       
        if(iCheck==1) {                         // don't interp if over/undershoot
            double maxV = *max_element(var.begin(), var.end());
            double minV = *min_element(var.begin(), var.end());
            if(maxV > upBuptr || minV < loBuptr) {
                y2 /= oldRhoU2; 
                for(int i=0; i<(int)var.size(); i++)
                    var[i] = y2;
            } 
        }
        
        //---------- rho
        
        for(int i=0; i<(int)var.size(); i++)
            newRho[i] = rhoU[i]/var[i];

        rho.erase(rho.begin()+isplt);
        rho.insert(rho.begin()+isplt, newRho.begin(), newRho.end());

        //---------- props (e.g., u v w enthalpy)
        // u, v, w, enthalpy need rho which is now done, so get enth here

        for(int k=0; k<iPropsYi; k++) {
        
            y2 = (*props[k])[isplt]*oldRhoU2;
            if(isplt > 0)
                y1 = (*props[k])[im]*oldRhoU1;
            if(isplt==0 && odtP->Lperiodic) 
                y1 = ((*props[k])[im]-odtP->pJump[k])*oldRhoU1;
            if(isplt < ngrd-1)
                y3 = (*props[k])[ip]*oldRhoU3;
            if(isplt==ngrd-1 && odtP->Lperiodic) 
                y3 = ((*props[k])[ip]+odtP->pJump[k])*oldRhoU3;

            iCheck = getABC(a,b,c, isplt, y1, y2, y3);
            interpTheVar(var, a,b,c, y1,y2,y3, isplt, iCheck, interCellPosns);
            for(int i=0; i<(int)var.size(); i++)
                var[i] /= rhoU[i];
            if(iCheck==1) {                         // don't interp if over/undershoot
                double maxV = *max_element(var.begin(), var.end());
                double minV = *min_element(var.begin(), var.end());
                if(maxV > upBprops[k] || minV < loBprops[k]) {
                    y2 /= oldRhoU2; 
                    for(int i=0; i<(int)var.size(); i++)
                        var[i] = y2;
                }   
            }   
            (*props[k]).erase((*props[k]).begin()+isplt);
            (*props[k]).insert((*props[k]).begin()+isplt, var.begin(), var.end());
        }
    }

    return;
}      


//////////////////////////////////////////////////////////////////////////////////////////////////

/** Called by splitCell, and also by itself (recursive).
 *  Given parabola coefficients, interpolate the cell to split into var values in 
 *  vector var.
 *  
 *  @param var    \output interpolated values of variable.
 *  @param a      \input parabola coefficient \fun{ a } as in \fun{ ax^2+bx+c }
 *  @param b      \input parabola coefficient \fun{ b } as in \fun{ ax^2+bx+c }
 *  @param c      \input parabola coefficient \fun{ c } as in \fun{ ax^2+bx+c }
 *  @param y1     \input variable value at point x1.
 *  @param y2     \input variable value at point x1.
 *  @param y3     \input variable value at point x1.
 *  @param isplt  \input index of cell being split.
 *  @param iCheck \input (but also reset for recusion): integer to flag properties of the problem.
 *  @param px     \input cell positions relative to left face.
 */
void anyline::interpTheVar(std::vector<double> &var, double &a, double &b, double &c, 
                           double y1, double y2, double y3,
                           int isplt, int iCheck, std::vector<double> px) {


    double xw = posf[isplt];
    double xe = posf[isplt+1];

    for(int i=0; i<(int)px.size(); i++)  // px are locations relative to left face so shift 
        px[i] += xw;                // to absolute position

    int    j = (int)px.size()-1;
    double d1, d2, d3;

    //---------- apply the interpolation: The integral over the profile 
    // in the interval (new cell) should equal the new value times the interval size
    
    if(a == 0 && b == 0)           // trival case: no interp for a constant prof. 
        for(int i=0; i<(int)var.size(); i++)
            var[i] = c;
    else if(a == 0) {              // easy case of a linear profile
        var[0]            = b*(px[0]*px[0]-xw*xw)/(2.0*(px[0]-xw)) + c;
        var[var.size()-1] = b*(xe*xe-px[j]*px[j])/(2.0*(xe-px[j])) + c;
        for(int i=1, im=0; i<(int)var.size()-1; i++, im++) 
            var[i] = b*(px[i]*px[i]-px[im]*px[im])/((px[i]-px[im])*2.0) + c;
    }
    else {                         // parabolic profile (general, so could just use this case)
        d1 = px[0]*px[0];
        d2 = xw*xw;
        d3 = px[0]-xw;
        var[0]  = a*(d1*px[0] - d2*xw)/(d3*3.0) + b*(d1-d2)/(d3*2.0) + c;

        d1 = xe*xe;
        d2 = px[j]*px[j];
        d3 = xe-px[j];
        var[var.size()-1] = a*(d1*xe - d2*px[j])/(d3*3.0) + b*(d1-d2)/(d3*2.0) + c;

        for(int i=1, im=0; i<(int)var.size()-1; i++, im++) {
            d1 = px[i]*px[i];
            d2 = px[im]*px[im];
            d3 = px[i] - px[im];
            var[i] = a*(d1*px[i] - d2*px[im])/(d3*3.0) + b*(d1-d2)/(d3*2.0) + c;
        }
    }

    //---------- the interpolated profile may not be good, depending on iCheck

    if(iCheck == 0)           // no need to check the interpolation
        return;

    if(iCheck == 2) {         // orig was monotonic, redo if interp makes a new min/max

        bool Lflag = false; 
        j = var.size()-1;

        if( (var[0]-y1)*(var[0]-var[1]) > 0.0 || (var[j]-y3)*(var[j]-var[j-1]) > 0.0 ) 
            Lflag = true;
        if(!Lflag) 
            for(int i=1, im=0, ip=2; i<(int)var.size()-1; i++, im++, ip++)
                if( (var[i]-var[im])*(var[i]-var[ip]) > 0.0) {
                    Lflag = true;
                    break;
                }
        if(Lflag) {
            iCheck = getABC(a,b,c,isplt, y1,y2,y3, true);
            for(int i=0; i<(int)px.size(); i++)
                px[i] -= xw;
            interpTheVar(var, a,b,c, y1,y2,y3, isplt, iCheck, px);
        }

    }

    return;

}

///////////////////////////////////////////////////////////////////////////////

/**Called by splitCell.  Routine finds an interpolating parabola to determine
 * variable values for the newly created cells.
 * grid:
 *                                                                                                                                                  <tt><pre>
 *            |    *     |     *    |     *    |
 *                      xw          xe
 *                 x1          x2        x3
 *                                                                                                                                                  </pre></tt>
 * parabolic interpolant: \fun{ y = ax^2 + bx + c }.
 * three degrees of freedom a,b,c, with 3 constaints: pass through points 1,3 and
 * \fun{ \int x_2 = y_2(x_e-x_w) }.
 * If a=0, then have a linear prof.  If a,b = 0 then have the usual no interp const profile
 * Return value is 0 if the abc values are fine as is, 1 if the profile has a min/max,
 * 2 if the profile is monotonic but the parabola has a min/max.
 *
 *  @param a \input parabola coefficient \fun{ a } as in \fun{ ax^2+bx+c }
 *  @param b \input parabola coefficient \fun{ b } as in \fun{ ax^2+bx+c }
 *  @param c \input parabola coefficient \fun{ c } as in \fun{ ax^2+bx+c }
 *  @param isplt \input index of cell being split.
 *  @param y1 \input variable value at point x1.
 *  @param y2 \input variable value at point x1.
 *  @param y3 \input variable value at point x1.
 *  @param Llinear \input flag to reduce to a conservative linear interpolation.
 *  @return integer flag of the state of the profile (0, 1, 2)
 */
int anyline::getABC(double &a, double &b, double &c, int isplt, 
                     double y1, double y2, double y3, bool Llinear) {

    int ip = isplt+1;
    int im = isplt-1;
    double x2 = pos[isplt];
    double xw = posf[isplt];
    double xe = posf[ip];
    double x1;
    double x3;

    if(im==-1){
        x3 = pos[ip];
        if(odtP->Lperiodic)
            x1 = pos[ngrd-1] - Ldomain;
        else {
            a = 0;
            b = 0;
            c = y2;
            return 0;
        }
        im = ngrd-1;
    }
    else if(ip==ngrd){
        x1 = pos[im];
        if(odtP->Lperiodic)
            x3 = pos[0] + Ldomain;
        else{
            a = 0;
            b = 0;
            c = y2;
            return 0;
        }
        ip = 0;
    }
    else {
        x1 = pos[im];
        x3 = pos[ip];
    }

    //---------------- Don't interpolate across a phase change,
    // e.g., 0 1 0 or 1 1 0 or 0 1 1

    if( fabs(phase[im]-phase[isplt]) > 0.01 &&
      fabs(phase[isplt]-phase[ip]) > 0.01 ) {
        a = 0;
        b = 0;
        c = y2;
        return 0;
    }
    else if( fabs(phase[im]-phase[isplt]) > 0.01 ) {
        a = 0;
        b = (y3-y2)/(x3-x2);
        c = y2-b*x2;
        return 0;
    }
    else if( fabs(phase[ip]-phase[isplt]) > 0.01 ) {
        a = 0;
        b = (y2-y1)/(x2-x1);
        c = y2-b*x2;
        return 0;
    }

    //---------------- check for linear bypass: the profile is monotnic but the 
    // parabola has an extremum that generates a new extremum.  Hence approximate
    // the cell profile with a line through point 2 and point 1 or 3 of min absolute slope

    if(Llinear) {  
        a = 0;
        b = (y3-y2)/(x3-x2);
        double dmb = (y2-y1)/(x2-x1);
        if(fabs(dmb) < fabs(b))
            b = dmb;
        c = y2-b*x2;
        return 0;
    }

    //---------------- central region

    if(y1 == y2 || y2 == y3) {    // no interp: prof. is like |___|___|---|  or |---|---|___|
        a = 0;
        b = 0;
        c = y2;
        return 0;
    }


    double e = (xe*xe*xe-xw*xw*xw)/(3.0*(xe-xw));
    double f = (xe*xe-xw*xw)/(2.0*(xe-xw));
    double dA  = x1*x1*(x3-f) - x1*(x3*x3-e) + x3*(x3*f-e);
    if(abs(dA) < 1E-14) {
        a = 0;
        b = 0;
        c = y2;
        return 0;
    }
    else {
        double dA1 = y1*(x3-f)          - x1*(y3-y2)         + y3*f-x3*y2;
        double dA2 = x1*x1*(y3-y2)      - y1*(x3*x3-e)       + x3*x3*y2-e*y3;
        double dA3 = x1*x1*(x3*y2-f*y3) - x1*(x3*x3*y2-e*y3) + y1*(x3*x3*f-e*x3);
        a = dA1/dA;
        b = dA2/dA;
        c = dA3/dA;
    }

    if( (y1 < y2 && y2 > y3) || (y1 > y2 && y2 < y3) ) {
        a = 0;  //doldb
        b = 0; //doldb
        c = y2; //doldb
        return 1;       
    }
    else {
        if(a==0)
            return 0;
        double xcrit = -b/(2*a);
        if(xcrit >= x3 || xcrit <= x1)
            return 0;
        else
            return 2;
    }

}

///////////////////////////////////////////////////////////////////////////////

/** Split two adjacent cells (i, i+1) and merge the center two cells to 
 * form three cells.
 *                                                                                                                                                   <pre><tt>
 * |   |   |   ==>   |  : | :  |   ==>   |  |  |  |
 *                                                                                                                                                   </tt></pre>
 * Called from adapt mesh to handle high changes in a variable phi
 *
 * @param i \input split this cell and i+1.
 * @param f1 \input lower face of new inner cell.
 * @param f2 \input upper face of new inner cell.
 */
void anyline::split2mergeMiddle(int i, double f1, double f2) {


    int ip  = i + 1;
    
    double d1 = posf[ip] - f1;
    double d2 = f2 - posf[ip];   
    double d3 = d1+d2;
    double w1, w2;
    double rho3, rhoU3, u3;

    posf[ip] = f2;
    posf.insert( posf.begin() + ip, f1 );

    pos[i]  = 0.5 * (posf[i] + f1);
    pos[ip] = 0.5 * (f2 + posf[ip+1+1]);          // note ip-->ip+1 after face insert
    pos.insert( pos.begin()+ip, 0.5*(f1+f2) );

    if(!odtP->Lspatial) {
        w1   = rho[i]  * d1/d3;
        w2   = rho[ip] * d2/d3;
        rho3 = w1+w2;

        w1 /= rho3;
        w2 /= rho3;
    }

    else {
        
        w1   = rho[i]  *(*uptr)[i]  * d1/d3;
        w2   = rho[ip] *(*uptr)[ip] * d2/d3;
        rhoU3 = w1+w2;

        u3    = (w1 *(*uptr)[i] + w2*(*uptr)[ip]) / rhoU3;

        w1 /= rhoU3;
        w2 /= rhoU3;

        rho3  = rhoU3 / (u3);

    }

    rho.insert(  rho.begin() +ip, rho3 ); 
    molec.insert( molec.begin()+ip, d3/(d1/molec[i] + d2/molec[ip]) );
    lambda.insert( lambda.begin()+ip, d3/(d1/lambda[i] + d2/lambda[ip]) );
    //Don't expect to be in here if phase is not constant across i, ip
    //Just us phase[i] value in new cell
    phase.insert( phase.begin()+ip, phase[i] );

    for(int k=0; k<nprops; k++) 
        (*props[k]).insert( (*props[k]).begin()+ip, (*props[k])[i]*w1+(*props[k])[ip]*w2 );
    
    
    ngrd++;
    ngrdf++;

}

///////////////////////////////////////////////////////////////////////////////

/**Output is useful for plotting when you want to assume cell contents are
 * peicewise constants.  This routine prints this profile.  Good for
 * looking at mesh adaption, etc.
 *
 * @param fname \input write to this file name.
 */
void anyline::outputProperties2(std::string fname) {

   ofstream ofile(fname.c_str()); 
   if(!ofile) 
       cout << "\n\n***************** ERROR OPENING FILE " << fname << endl << endl;
   
   ofile << "# grid points = " << ngrd;
   ofile << "\n# Domain Size = " << Ldomain;
   if(fabs(posf[ngrd] -posf[0] - Ldomain) > 1.0E-6)
       ofile << "\n# last posf-first posf != Ldomain, last posf = " << posf[ngrd];
   ofile << endl;
   ofile << "#  1_posf            "
         << "2_rho             "
         << "3_molec           "
         << "4_lambda          "
         << "5_phase           ";
   for(int k=0; k<nprops; k++) 
       ofile << k+6 << "_" << propNames[k] << "            ";

   ofile << scientific;
   ofile << setprecision(10);
  
   int i=0;
   ofile << endl;
   ofile << setw(19) << posf[i] 
         << setw(19) << rho[i]
         << setw(19) << molec[i]
         << setw(19) << lambda[i]
         << setw(19) << phase[i];
   for(int k=0; k<nprops; k++) 
       ofile << setw(19) << (*props[k])[i];

   for(i=1; i<ngrd; i++) {
       ofile << endl;
       ofile << setw(19) << posf[i] 
             << setw(19) << rho[i-1]
             << setw(19) << molec[i-1]
             << setw(19) << lambda[i-1]
             << setw(19) << phase[i-1];
       for(int k=0; k<nprops; k++) 
           ofile << setw(19) << (*props[k])[i-1];
       ofile << endl;
       ofile << setw(19) << posf[i] 
             << setw(19) << rho[i]
             << setw(19) << molec[i]
             << setw(19) << lambda[i]
             << setw(19) << phase[i];
       for(int k=0; k<nprops; k++) 
           ofile << setw(19) << (*props[k])[i];
   }

   i=ngrd;
   ofile << endl;
   ofile << setw(19) << posf[i] 
       << setw(19) << rho[i-1]
       << setw(19) << molec[i-1]
       << setw(19) << lambda[i-1]
       << setw(19) << phase[i-1];
   for(int k=0; k<nprops; k++) 
       ofile << setw(19) << (*props[k])[i-1];

   ofile.close();
}

///////////////////////////////////////////////////////////////////////////////

/** Copy properties from one index to another in vectors rho, molec.
 *
 *  @param from \input copy data in this index into index "to".
 *  @param to   \input copy data from index "from" into this index.
 */
void anyline::copyProps(int from, int to) {


  rho[to] = rho[from];
  molec[to] = molec[from];
  lambda[to] = lambda[from];

}

///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//subdomain decompsition
void anyline::addL(double inL){
    int i;
 for(i=0;i<ngrd;i++){
  pos[i]+=inL;
  posf[i]+=inL;
 }
posf[ngrdf-1]+=inL;

}







