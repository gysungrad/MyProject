/**
 * @file table.cc
 * Source file for class table
 */

#include "table.h"
#include "TableReader.h"
#include "processor.h"
#include "odtline.h"
#include "ETA.h"
#include <math.h>

using namespace std;

extern processor proc;

///////////////////////////////////////////////////////////////////////////////

/** Constructor to use.
 *  @param _odtl the odtline object to use.
 */

table::table(odtline* _odtl) {

	odtl 		   = _odtl;                                  ///< A pointer to the odtline

    TableReader tr;
    if(odtl->odtP->ItableLookup == odtl->odtP->TT_EQUILIBRIUM)
        tr = TableReader("../input/lookupTable/eqProf.conf", odtl);
    else if(odtl->odtP->ItableLookup == odtl->odtP->TT_FLAMELET){
        tr = TableReader("../input/flamelet/data/flmlttable.conf", odtl);
    }
    
    mixf           = tr.mixf;
	enth           = tr.enth;
    hsensVals      = tr.hsens;
    eqTables       = tr.Table;
    hlosses        = tr.hlosses;
    keys           = tr.keys;
    numColumns     = tr.numColumns;
    numRows        = tr.numRows;
    index          = tr.index;
    numSpecies     = tr.numspecies;

	uniform        = checkUniformity();                      

    //Make sure that the numSpecies is equal to the num Species in the odtl;
    //TODO: Make Vector of indexes so that the vector updates correctly
    if(numSpecies != (int)odtl->yspc.size()) {
        cout << __FILE__ << " (Line " <<__LINE__<< ") - WARNING! DIFFERENT NUMBER OF SPECIES\n(" 
                << numSpecies << " != " << odtl->yspc.size() << ")" << endl;
        
        for(int i = 0; i < numSpecies; i++) {
            if(index.species_mass_index > 0)
                cout << " " << keys[i+index.species_mass_index];
            else
                cout << " " << keys[i+index.species_mole_index];
        }
        cout << endl;
        for(int i = 0; i < (int)odtl->yspc.size(); i++)
            cout << " " << odtl->spNames[i];
        cout << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////

/** Given a mixture fraction and enthalpy value, this function interpolates
 *  the 3D table to return a vector of an estimated values.
 *
 * @param mf    \input The desired mixture fraction
 * @param h     \input The desired enthalpy
 * @return a vector containing the interpolated values
 */

vector<double> table::getVals(const double mf, const double h) {
    indices mfi;                //Indices of the mixture fractions bounding mf
    indices hli;                //Indices of the mixture fraction bounding the hloss
    double hloss;               //Heat loss used in table

    mfi = getIndices(mf, mixf);
    hloss = calcHLossFromIndices(mf, h, mfi);
    hli = getIndices(hloss, hlosses);

    //OK Now we have both the heat loss and mixf,
    //All that's left to do is interpolate all the
    //values between the 4 hloss,mixf values
    vector<double> yVals(numColumns);
    for (int i = 0; i < numColumns; i++)
            yVals[i] = bilinearInterpolate(mfi.d, mfi.u, hli.d, hli.u, mf, hloss, i);
    
    return yVals;
}

///////////////////////////////////////////////////////////////////////////////


/** Estimates the mole fractions (Y) of the species based on a given
 * mixture fraction and enthalpy
 *
 * @param mf    \input Mixture fraction
 * @param h     \input Enthalpy
 * @return a vector containing the interpolated Y values of the species
 */

vector<double> table::getYVals(const double mf, const double h){

    indices mfi;                //Indices of the mixture fractions bounding mf
    indices hli;                //Indices of the mixture fraction bounding the hloss
    double hloss;               //Heat loss used in table

    mfi = getIndices(mf, mixf);
    hloss = calcHLossFromIndices(mf, h, mfi);
    hli = getIndices(hloss, hlosses);

    //OK Now we have both the heat loss and mixf,
    //All that's left to do is interpolate all the
    //values between the 4 hloss,mixf values
    vector<double> yVals(numSpecies);
    for (int i = 0; i < numSpecies; i++)
        yVals[i] = bilinearInterpolate(mfi.d, mfi.u, hli.d, hli.u, mf, hloss, i + index.species_mole_index);

    return yVals;
}

///////////////////////////////////////////////////////////////////////////////
/** Returns the mole species at all the points in the grid.
 *  @return a vector<double> containing the mole species at all points in the grid.
 */

vector<vector<double> > table::getXMoleSp(){
    vector<vector<double> > ret;
    ret.resize(odtl->ngrd);
    for(int i = 0; i < odtl->ngrd; i++) {
        ret[i] = getSpeciesByMoleAtGridPoint(i);
    }

    return ret;
}

///////////////////////////////////////////////////////////////////////////////

/** Returns a vector containing species by mole fraction at a given position
 *      along the grid
 *
 *  @param iPos \input the position along the grid to be considered
 *  @return a vector<double> containing the values of all the species at the given
 *              position
 */

vector<double> table::getSpeciesByMoleAtGridPoint(int iPos) {
    double  mf     =  odtl->eta[0][iPos];     //The mixture fraction at the current position
    double  h      =  odtl->eta[1][iPos];     //The enthalpy at the current position
    indices mfi;                            //Indices of the mixture fractions bounding mf
    indices hli;                            //Indices of the mixture fraction bounding the hloss
    double  hloss;                          //Heat loss used in table

    mfi = getIndices(mf, mixf);
    hloss = calcHLossFromIndices(mf, h, mfi);
    hli = getIndices(hloss, hlosses);

    //OK Now we have both the heat loss and mixf,
    //All that's left to do is interpolate all the
    //values between the 4 hloss,mixf values
    vector<double> yspc(numSpecies);
    for(int i = 0; i < numSpecies; i++) {
        yspc[i] = bilinearInterpolate(mfi.d, mfi.u, hli.d, hli.u, mf, hloss, index.species_mole_index + i);
    }
    return yspc;
}

///////////////////////////////////////////////////////////////////////////////

/** Reads the absorption Coefficients from the table. Made to be compatible with
 *  the radiation class and use a table lookup for calculated kabs values instead.
 *
 *  @return a vector<double> containing the computed absorption coefficients.
 */

//todo throw an error if there is no index.kabs_index

vector<double> table::getAbsorptionCoeffs(){
    vector<double> kabs(odtl->ngrd);
    for(int i = 0; i < odtl->ngrd; i++) {
        double  mf     =  odtl->eta[0][i];     //The mixture fraction at the current position
        double  h      =  odtl->eta[1][i];     //The enthalpy at the current position
        indices mfi;                            //Indices of the mixture fractions bounding mf
        indices hli;                            //Indices of the mixture fraction bounding the hloss
        double  hloss;                          //Heat loss used in table

        mfi = getIndices(mf, mixf);
        hloss = calcHLossFromIndices(mf, h, mfi);
        hli = getIndices(hloss, hlosses);

        //OK Now we have both the heat loss and mixf,
        //All that's left to do is interpolate all the
        //values between the 4 hloss,mixf values
        kabs[i] = bilinearInterpolate(mfi.d, mfi.u, hli.d, hli.u, mf, hloss, index.kabs_index);
    }
    return kabs;

}

///////////////////////////////////////////////////////////////////////////////

/** Estimates the value of a specified property at a given grid point.
 *
 * @param iPos      \input the grid point to be considered
 * @param prop_indx \input the index of the property/species to be estimated.
 * @return the estimated value of the property prop_indx at the given grid point
 */


double  table::getValAtGridPoint(int iPos, int prop_indx){

    double  mf     =  odtl->eta[0][iPos];     //The mixture fraction at the current position
    double  h      =  odtl->eta[1][iPos];     //The enthalpy at the current position
    indices mfi;                            //Indices of the mixture fractions bounding mf
    indices hli;                            //Indices of the mixture fraction bounding the hloss
    double  hloss;                          //Heat loss used in table

    mfi = getIndices(mf, mixf);
    hloss = calcHLossFromIndices(mf, h, mfi);
    hli = getIndices(hloss, hlosses);

    //OK Now we have both the heat loss and mixf,
    //All that's left to do is interpolate all the
    //values between the 4 hloss,mixf values
    return bilinearInterpolate(mfi.d, mfi.u, hli.d, hli.u, mf, hloss, prop_indx);

}

///////////////////////////////////////////////////////////////////////////////

/** Estimates the value of specified property at a given grid point.
 *
 * @param iPos \input the grid point to be considered
 * @param prop \input a string containing the property/species to look up.
 * @return the calculated value at the given grid point and index
 */


double  table::getValAtGridPoint(int iPos, string prop){
    return getValAtGridPoint(iPos, getKeyIndex(prop));

}

///////////////////////////////////////////////////////////////////////////////


/** Returns an estimated value of the column using the given mixture
 *  fraction and enthalpy
 *
 * @param mf    \input The desired mixture fraction
 * @param h     \input The desired enthalpy
 * @param key   \input String containing the key of the desired column (ex: "H2")
 * @return The interpolated value from the internal table
 */

double table::getValOf(const double mf, const double h, const string key){
    return getValOf(mf,h,getKeyIndex(key));
}

///////////////////////////////////////////////////////////////////////////////


/** Returns an estimated value of the property specified by prop_indx using the given mixture
 *  fraction and enthalpy
 *
 * @param mf         \input Mixture fraction
 * @param h          \input Enthalpy
 * @param prop_indx  \input The index of the desired key
 * @return The estimate value for the property specified by prop_index
 *         based on the table.
 */

double table::getValOf(const double mf, const double h, int prop_indx){
    indices mfi;                //Indices of the mixture fractions bounding mf
    indices hli;                //Indices of the mixture fraction bounding the hloss
    double hloss;               //Heat loss used in table

    if (prop_indx < 0)
        perror("INVALID INDEX!!!");

    mfi = getIndices(mf, mixf);
    hloss = calcHLossFromIndices(mf, h, mfi);
    hli = getIndices(hloss, hlosses);
    return bilinearInterpolate(mfi.d, mfi.u, hli.d, hli.u, mf, hloss, prop_indx);
}

///////////////////////////////////////////////////////////////////////////////


/** Returns the property index of the specified key based on the table's
 *  internal key vector
 *
 *  @param key  \input String containing the key of the desired column (ex: "H2")
 *  @return the key's property index
 */

int table::getKeyIndex(const string key){
    int lastIndex = keys.size() - 1;             //The last key's index

    //Iterate through keys
    for (int i = 0; i <= lastIndex; i++){

        //Check each key to see if it is equal
        if(keys[i].compare(key) == 0){

            //If it is return it
            return i;
        }

    }

    //Key not found

    return -1;
}

///////////////////////////////////////////////////////////////////////////////


/** Given a mixture fraction and heat loss value, this function interpolates
 * the 3D table to return a vector of an estimated values.
 *
 * @param mf    \input mixture fraction
 * @param hloss \input heat loss
 * @return a vector containing the estimated values for the species and properties
 *   at the given mixture fraction and heat loss
 */

vector<double> table::getValsGivenHeatLoss(const double mf, const double hloss) {
    int mfu;                        // The upper row found from mixture fraction for interpolation ('MF Up')
    int mfd;                        // The lower row from the mixture fraction for interpolation ('MF Down')
                                    // Note: mfu < mfd because the upper row has a smaller index
    int hlu;                        // The upper row index found from heat loss for interpolation ('Hloss Up')
    int hld;                        // The lower row index found from heat loss for interpolation ('HLoss Down')
                                    // Note: hlu < hld because the upper row has a smaller index


    //find mixture fraction indices

    if (mf <= 0.001) {                 //First case, 0 or very close

        // Set both to first element
        mfd = 0;
        mfu = mfd;

    } else if (mf >= 1) {              //Second Case, 1 or greater (computational error)

        // Set both to last element
        mfd = mixf.size() - 1;
        mfu = mfd;

    } else {                           //Last case, between 0 and 1

        //find the indices of the two mixture fractions
        //such that the given mf is between mfu and mfd

        for (int i = 0; i < (int)mixf.size(); i++) {

            if (mixf[i] >= mf) {
                mfd = i;
                mfu = i - 1;
                break;

            }
        }
    }

    //find the heat loss indices

    if (hloss <= hlosses[0]){                        //First case, Below the lowest given Heat Loss

        hlu = hld = 0;

    } else if (hloss > hlosses[hlosses.size() - 1]){ //Second Case, Above the Highest Given HeatLoss

        hlu = hld = hlosses.size() - 1;

    } else {                                         //Last Case, Between given values

        //find the indices of the two hlosses
        //such that the given hloss is between hlu and hld

        for (int i = 0; i < (int)hlosses.size(); i++) {

            if (hlosses[i] >= hloss) {

                hld = i;
                hlu = i - 1;
                break;

            }
        }
    }

    //OK Now we have both the heat loss and mixf,
    //All that's left to do is interpolate all the
    //values between the 4 hloss,mixf values
    vector<double> yVals(numColumns - 3);
    for (int i = 0; i < numColumns - 3; i++) {
        yVals[i] = bilinearInterpolate(mfd, mfu, hld, hlu, mf, hloss, i);
    }
    return yVals;

}

///////////////////////////////////////////////////////////////////////////////

/**
 * Does a binary search on a vector to find the indices where the given value is bounded
 *
 *  @param vec   \input A vector containing a sorted list to be searched
 *  @param value \input The value to be bounded.
 *  @return an table::indices object that contains the indices to the values in vec that bound the value
 */

table::indices table::binaryIndicesSearch(vector<double> vec, double value) {
    indices in;                         //The indices to the values in vec that bound the value
    int     min     = 0;                //The lower bound
    int     max     = vec.size() - 1;   //The upper bound
    int     mid;                        //The midpoint of the upper and lower bounds
    double  minV    = vec[min];         //The value of the lower bound in vec
    double  maxV    = vec[max];         //The value of the upper bound in vec
    double  midV;                       //The value of the midpoint in vec

    //The goal is to get min and max 1 apart, these two values are bounding
    while(max - min > 1)
    {
        //Calculate midpoint and look up value
        mid = (max + min) / 2;
        midV = vec[mid];


        if (value > midV) {     //If the midpoint is bigger move the min up

            min = mid;
            minV = midV;

        } else {                //If the midpoint is smaller move the mid down

            max = mid;
            maxV = midV;

        }
    }

    in.d = max;
    in.u = min;

    return in;
}

///////////////////////////////////////////////////////////////////////////////


/** Calculates the hLoss using the given mixture fraction indices
 *
 *  @param mf   \input  mixture fraction for the hloss to be calculated for
 *  @param h    \input  enthalpy for the hloss to be calculated for
 *  @param mfi  \input  An indices struct containing the indices of the mixfs that bound the mf
 *  @return The calculated hloss
 */

double table::calcHLossFromIndices(double mf, double h, indices mfi){

    double hsens;                   //Sensible heat for given mf
    double had;                     //Adiabatic enthalpy for given mf

    //Check if interpolation for hsens is needed

    if (mfi.u == mfi.d)
        hsens = hsensVals[mfi.d];
    else
        hsens = linearInterpolate(mixf[mfi.u], hsensVals[mfi.u], mixf[mfi.d],
                hsensVals[mfi.d], mf);


    //Adiabatic heat
    had = odtl->strm->h0 * (1 - mf) + odtl->strm->h1 * mf;

    //Heat Loss
    return (had - h) / hsens;

}

///////////////////////////////////////////////////////////////////////////////


/** Checks if mixf and hlosses values are uniform
 *
 *  @return <DFN>True</DFN>  iff mixf and hlosses are uniform
 */

bool table::checkUniformity(){

    //Check mixf
    if(!uniformVector(mixf))
        return false;

    //Check hlosses
    if(!uniformVector(hlosses))
        return false;

    return true;

}

///////////////////////////////////////////////////////////////////////////////


/** Searches a vector to find the indices of the values that the given
 *      value is between
 *
 *  @param value    \input The value to be searched for
 *  @param vec      \input The vector which the value correlates with
 *  @return an indices object containing the upper and lower bound indices
 */

table::indices table::getIndices(double value, vector<double> vec){
    indices in;                 //The indices for the values in vec that bound the value
    in.d = 0;
    in.u = 0;

   //find the heat loss indices

    if (value <= vec[0] + lbError){                        //First case, Below the lowest given Heat Loss

        in.u = in.d = 0;

    } else if (value >= vec[vec.size() - 1]){ //Second Case, Above the Highest Given HeatLoss

        in.u = in.d = vec.size() - 1;

    } else {                                        //Last Case, Between given values

        if(uniform){

               double dVal=0;            //Difference between value and the first element in the vector
               double di=0;              //Change between value in indices
               double i=0;               //How many steps to go from first value to the desired value

               //Calculate values
               dVal = value - vec[0];
               di   = vec[1] - vec[0];
               i = dVal / di;
               in.u = (int)floor(i);
               in.d = in.u + 1;
               if(vec[in.u] > value){
                   --in.u;
                   --in.d;
               }
               if(vec[in.d] < value){
                   ++in.u;
                   ++in.d;
               }

       } else {
            in = binaryIndicesSearch(vec, value);

        }
    }

    return in;

}

///////////////////////////////////////////////////////////////////////////////

/** This function linearly interpolates data between two 2D points.
 *
 *  @param x1     \input The first position along the x
 *  @param x1_val \input The value at x1
 *  @param x2     \input The second position along the x axis
 *  @param x2_val \input The value at x2
 *  @param x      \input The x position for which the value is to be estimated
 *  @return if x1 and x2 are different, it returns the interpolated value, otherwise
 *      it returns the average of {x1_val,x2_val}
 */

double table::linearInterpolate(double x1, double x1_val, double x2, double x2_val,
		double x) {

    //Case where x values are equal
	if (x2 == x1) {
		if (x2_val != x1_val)
			return (x2_val + x1_val) / 2.;
		else
			return x2_val;
	}

	//Case where x values are different
	//y = y1 + (x-x1)*dy/dx

	return x1_val + (x - x1) * (x2_val - x1_val) / (x2 - x1);
}

///////////////////////////////////////////////////////////////////////////////


/** Bilinear interpolation using previously
 *      found hloss and mf indices
 *
 *<tt><pre>
 *   hlu O * * * Y O
 *    |  * * * * | *
 * hloss --------X *
 *    |  * * * * | *
 *    |  * * * * | *
 *   hld O * * *-Z-O
 *      mfu - - -|mfd
 *           mf
 *</pre></tt>
 *     Each \c 'O' has a given value, the \c 'X' is the desired value so this performs
 *     an interpolation to find out an estimated value of \c X. In this implementation,
 *     This is done first by interpolating horizontally and vertically to find the
 *     values of \c 'Y' and \c 'Z' and then interpolating
 *     between them to find <code>X</code>
 *
 *  @param mfd          \input The upper row index that the desired mf falls between
 *  @param mfu          \input The lower row index that the desired mf falls between
 *  @param hld          \input The upper row index that the desired hl falls between
 *  @param hlu          \input The lower row index that the desired hl falls between
 *  @param mf_actual    \input The actual mf
 *  @param hloss_actual \input The actual hloss
 *  @param index        \input The index of the species/property to be calculated
 *  @return The interpolated value based on the mixture fractions and heat loss values
 *      of the indexed species/property
 */

double table::bilinearInterpolate(int mfd, int mfu, int hld, int hlu,
		double mf_actual, double hloss_actual, int index) {
    
    /****************--For extrapolation--*******************/
    if(hld==hlu){ 
        if(hld==(int)hlosses.size()-1)
            hlu=hlosses.size()-2;
    }
    if(mfd==mfu){
         if(mfd==(int)mixf.size()-1)
            mfu=mixf.size()-1;
    }
    /***********************************************************/
	//First get the four values
	double dd = eqTables[hld][mfd][index];             //Lower Left Species Value
	double du = eqTables[hld][mfu][index];             //Lower Right Species Value
	double ud = eqTables[hlu][mfd][index];             //Upper Left Species Value
	double uu = eqTables[hlu][mfu][index];             //Upper Right Species Value

	//Now we get the know heat loss and mixture fraction values
	double mfu_val = mixf[mfu];                        //Upper mixf value
	double mfd_val = mixf[mfd];                        //Lower mixf value
	double hlu_val = hlosses[hlu];                     //Upper hloss value
	double hld_val = hlosses[hld];                     //Lower hloss value

	//Interpolate between all the points
	//First interpolate between mfs
	double hlumf = linearInterpolate(mfu_val, uu, mfd_val, ud, mf_actual); //upper interpolated mixf value (Y on pic above)
	double hldmf = linearInterpolate(mfu_val, du, mfd_val, dd, mf_actual); //lower interpolated mixf value (Z on pic above)

	// Last Interpolate between these heatLosses
	return linearInterpolate(hlu_val, hlumf, hld_val, hldmf, hloss_actual);

}

///////////////////////////////////////////////////////////////////////////////


/** Checks a vector to see if its values are uniform, that is, that
 *      the differences don't change between index pairs.
 *
 *  @param vec  \input A vector to be checked for uniformity
 *  @return <DFN>True</DFN> if the vector is uniform
 */

bool table::uniformVector(vector<double> vec){
    int    count        = vec.size();           //Number of elements in the vector
    double diffDiffSum  = 0;                    //The running sum of the difference between the differences
    double lDiff        = 0;                    //The last difference calculated
    double diff         = 0;                    //The current difference
    double diffDiff     = 0;                    //The difference between two differences
    double ddAvg;                               //Average Difference between differences

    //First calculate the first difference to prepare for looping
    lDiff = vec[1] - vec[0];

    //Check the rest
    for(int i = 2; i < count; i++){
        diff = vec[i] - vec[i-1];

        //Absolute value of error so it only grows
        diffDiff = fabs(diff-lDiff);
        //cout << "dd: " << diffDiff << endl;

        //Quit if a difference quickly becomes too large
        if(diffDiff > maxDeviation)
            return false;
        diffDiffSum += diffDiff;
    }

    //Get the average (larger data sets won't be punished)
    ddAvg = diffDiffSum / (double)count;
    //cout << "ddAvg: " << ddAvg << endl;

    //If too large, not uniform, otherwise Uniform
    if (ddAvg > maxDeviation)
        return false;

    return true;

}
