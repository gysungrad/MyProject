/**
 * @file adaptMesh.h
 * Header file for class adaptMesh
 */

#ifndef ADAPTMESH_H
#define ADAPTMESH_H

#include "anyline.h"
#include "odtParam.h"
#include <vector>

///////////////////////////////////////////////////////////////////////////////

/** 
 * @class adaptMesh
 * 
 * Mesh adaptor class.  This class will adapt odtlines or scalines etc.  
 *  The specific line to adapt is set through a pointer.
 *  A single variable profile is used for adaption of a given line (e.g. velocity
 *  or mixture fraction).  The mesh adapter has several parts (see adaptODTgrid).
 *  This is meant to be called after eddy events or significant diffusion.  The
 *  adaption inserts or removes cells depending on cell size, gradients, curvature,
 *  and neighboring cell sizes.
 *  
 *  @author David O. Lignell
 */

#ifndef OLDMESHER

using namespace std;
//NEW MESHER
class adaptMesh {

    ////////////////////// DATA MEMBERS /////////////////////

    public:

        anyline             *anyl;     ///< pointer to odt line to adapt
        odtParam            *odtP;     ///< pointer to the odt parameters

        vector<vector<double> *> phi;  ///< vector of pointers to vector<double> profiles to adapt
        vector<double> *bdy;           ///< pointer to phase variable

        vector<double> dx;             ///< vector of cell sizes
        vector<int>    mark;           ///< dummy small cell index array for sorting
        vector<double> theta;          ///< vector of angles for curvature

        vector<double> xf;             ///< vector of cell face positions
        vector<vector<double> > yf;    ///< vector of cell values
        vector<double> xnf;            ///< vector of new cell face positions
        vector<double> X;              ///< vector of cell center positions

        int                 ngrd;      ///< same as odtl->ngrd
        int                 ngrdf;     ///< same as odtl->ngrdf

        int                 iLower;    ///< region of grid to adapt (the cell w/ left eddy edge)
        int                 iUpper;    ///< region of grid to adapt (the cell w/ right eddy edge)
        double              posLower;  ///< physical region of eddy (lower bound)
        double              posUpper;  ///< physical region of eddy (upper bound)


        ////////////////////// MEMBER FUNCTIONS  /////////////////////

        void adaptGrid(int iLower, int iUpper);       // Public interface
        
        inline bool operator()(const int &a, const int &b) const {
            return dx[mark[a]] < dx[mark[b]];            
        }                                                // Functor used for sorting

    private:

        void adaptGrid_details(int iLower, int iUpper);  // called by adaptGrid
        void mergeSmallCells();                          // 1 of 5 adapt operations
        void mergeSmallCellsMP();                        // 1 of 5 adapt operations
        void impose2point5rule();                        // 4 of 5 adapt operations
        void splitLargeCells();                          // 5 of 5 adapt operations

        void fix2point5offender(int mPos, const int &iglobal);
        void setDxArray();

        int  findPos(const vector<double> &x, const double val, int &istart);
        void interp1pt(const vector<double> &x, const vector<double> &y,
                       const double &xval, double &yval, int &istart);
        void interpVec(const vector<double> &x, const vector<double> &y,
                       const vector<double> &xn, vector<double> &yn);
        double calcDistance(const vector<double> &x,
                            const vector<vector<double> > &y,
                            vector<double> &sDist);
        void set_iLower_iUpper();



        ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:

        adaptMesh(anyline *anyll, odtParam *odtPP, vector<vector<double> *> phiP, bool Ladpt);
#ifdef COMPSCI
        inline void init(anyline *anyll, odtParam *odtPP, bool Ladpt);
#endif
        adaptMesh(){};
        ~adaptMesh(){ anyl=0; odtP=0; bdy=0;}

};
#ifdef COMPSCI
inline void adaptMesh::init(anyline *anyll, odtParam *odtPP, bool Ladpt) {
    
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
    
}
#endif

#else

using namespace std;
//OLD MESHER
class adaptMesh {

    ////////////////////// DATA MEMBERS /////////////////////

    public:

        anyline             *anyl;      ///< pointer to odt line to adapt
        odtParam            *odtP;      ///< pointer to the odt parameters

        vector<double> *phi;       ///< pointer to variable to adapt on (make sure this is in anyline::props so it gets adapted)
        vector<double> *bdy;       ///< pointer to phase variable

        vector<double> dx;         ///< vector of cell sizes
        vector<int>    mark;       ///< dummy small cell index array for sorting
        vector<double> theta;      ///< vector of angles for curvature

        int                 ngrd;       ///< same as odtl->ngrd
        int                 ngrdf;      ///< same as odtl->ngrdf

        double              largeGrad;  ///< largest allowed difference in grad var phi
        double              smallGrad;  ///< smallest allowed difference in grad var phi
        double              largeCurv;  ///< largest allowed difference in theta
        double              smallCurv;  ///< smallest allowed difference in theta

        int                 iLower;     ///< region of grid to adapt (the cell w/ left eddy edge)
        int                 iUpper;     ///< region of grid to adapt (the cell w/ right eddy edge)
        double              posLower;   ///< physical region of eddy (lower bound)
        double              posUpper;   ///< physical region of eddy (upper bound)


        ////////////////////// MEMBER FUNCTIONS  /////////////////////

        void adaptGrid(int iLower, int iUpper);       ///< Public interface
        
        inline bool operator()(const int &a, const int &b) const {
            return dx[mark[a]] < dx[mark[b]];            
        }                                                ///< Functor used for sorting

    private:

        void mergeSmallCells();                          ///< 1 of 5 adapt operations
        void mergeSmallGradsCurvs();                     ///< 2 of 5 adapt operations
        void splitHighGradsCurvs();                      ///< 3 of 5 adapt operations
        void impose2point5rule();                        ///< 4 of 5 adapt operations
        void splitLargeCells();                          ///< 5 of 5 adapt operations

        void split2mergeMiddle(int i, double f1, double f2);
        void merge2cells(int imrg);
        void mergeSeveralCells(int imrg, int ncells);
        void setLargeSmallGradCurv();
        void mergeTheMergeGroup(int mPos, int gSize, bool LdoAll);
        void fix2point5offender(int mPos, const int &iglobal);
        void setDxArray();
        bool testSmallGradAndCurv(const int &i);
        bool testCellLargeGradCurv(int i, vector<double> &newFacesL,
                                          vector<double> &newFacesH);
        void setThetaArray();                            ///< Used for curvature calculation

        void interpAdtpVar2p5(int is, int ie);

        ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:

        adaptMesh(anyline *anyll, odtParam *odtPP, vector<vector<double> *> phiP);
        adaptMesh(){}

};

#endif

#endif

