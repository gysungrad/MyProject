/**
 * @file anyline.h
 * Header file for class anyline
 */

#ifndef ANYLINE_H
#define ANYLINE_H

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "odtParam.h"

///////////////////////////////////////////////////////////////////////////////

/**
 * @class anyline
 * 
 * A base class for constructing various types of ODT lines. Especially
 * ODT lines consisting of the basic properties, like velocity, as 
 * well as scalar lines.  Generally, create odtlines or scalines.
 *  
 *  @author David O. Lignell
 */

class anyline {

  ////////////////////// DATA MEMBERS /////////////////////

    public:
    
        int                                 ngrd;         ///< number of grid points
        int                                 ngrdf;        ///< number of face points (ngrd+1)
        double                              Ldomain;      ///< domain length
        std::vector<double>                 pos;          ///< vector of cell center positions
        std::vector<double>                 posf;         ///< vector of face positions
        std::vector<double>                 rho;          ///< vector of densities (kg/m3)
        std::vector<double>                 molec;        ///< vector of diffusion coeff (e.g. mu)
        std::vector<double>                 lambda;       ///< thermal conductivity
        std::vector<double>                 phase;        ///< vector of phase identity

        int                                 nprops;       ///< number of props (vars) in props 
        std::vector<std::vector<double>* >  props;        ///< pointer to a vector of pointers (species should come last)
        std::vector<std::vector<double> >   bcprops;      ///< boundary conditions of each prop in props
        std::vector<double>                *uptr;         ///< a pointer to streamwise velocity for spatial formulation
        std::vector<std::string>            propNames;    ///< props should contain quantities per unit mass like Yi
                                                          ///< since rho*props is used as the conserved quantity
                                                          ///< so, things like u,v,w for momentum, mixf, Yi
        int                                 nspc;         ///< number of chemical species
        bool                                LpropsHasYi;  ///< Flag (see splitCell routine)
        int                                 iPropsYi;     ///< index location of Yi in props

        odtParam                           *odtP;         ///< pointer to odtParam object

        double                              upBrho;       ///< upper and lower bounds used in cell splitting
        double                              loBrho;       ///< e.g., mixture fraction should stay in 0:1
        double                              upBmolec;     ///< like upBrho
        double                              loBmolec;     ///< like loBrho
        double                              upBlambda;    ///< like upBrho
        double                              loBlambda;    ///< like loBrho
        double                              upBuptr;      ///< like upBrho
        double                              loBuptr;      ///< like loBrho
        std::vector<double>                 upBprops;     ///< like upBrho
        std::vector<double>                 loBprops;     ///< like loBrho

        std::vector<double>                 adptVar;    ///< variable to adapt on

        ////////////////////// MEMBER FUNCTIONS  /////////////////////

        int  linePositionToIndex(double position, bool LowSide=0);
        void insertEddy(anyline &eddyLine, int &iStart, int &iEnd,
                        double &leftEdge, double &rightEdge, bool Lwrap=false);
        void eraseCell(int ieras);
        void merge2cells(int imrg);
        void mergeSeveralCells(int ifirst, int ncells);
        void splitCell(int isplt, int nsplt, std::vector<double> &interCellPosns, 
                       bool Linterp=false); 
        void split2mergeMiddle(int i, double f1, double f2);

        virtual void outputProperties(std::string fname);
        virtual void outputProperties2(std::string fname);
        virtual void readProperties(std::string fname);

        void copyProps(int from, int to);

	void addL(double inL);
	
    private:

        void splitCellTemporal(int isplt, int nsplt, std::vector<double> &interCellPosns);
        void splitCellSpatial( int isplt, int nsplt, std::vector<double> &interCellPosns);

        int getABC(double &a, double &b, double &c, int isplt, 
                   double y1, double y2, double y3, bool Llinear=false);
        void interpTheVar(std::vector<double> &var, double &a, double &b, double &c, 
                          double y1, double y2, double y3, 
                          int isplt, int iCheck, std::vector<double> interCellPosns);


        ////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////

    public:

        anyline(int npts=5, double Ld=1.0);               // constructor
        virtual ~anyline() { }  // destructor


};

#endif
