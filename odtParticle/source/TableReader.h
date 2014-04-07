/**
 * @file TableReader.h
 * Header file for class TableReader
 */

#ifndef TABLEREADER_H
#define TABLEREADER_H


////////////////////// INCLUDE ////////////////////////////
#include <string>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "odtline.h"
#include "table.h"

using namespace std;
class ETA_Table;

///////////////////////////////////////////////////////////////////////////////

/** A class that reads user-generated tables to use with ODT.
 *  The user first must describe the table in an exterior configuration
 *  file as follows:                                                                   <br>
 *  <DFN>                                                                              <br>
 *  \#COMMENTS                                                                         <br>
 *                                                                                     <br>
 *  FILES                                                                              <br>
 *  (file_1) (property)                                                                <br>
 *  (file_2) (property)                                                                <br>
 *  ...                                                                                <br>
 *  (file_n) (property)                                                                <br>
 *  END                                                                                <br>
 *                                                                                     <br>
 *  <VAR>Note: The FILES tag is required</VAR>                                         <br>
 *                                                                                     <br>
 *  HSENS_FILE                                                                         <br>
 *  (sensible enthalpy data file)                                                      <br>
 *                                                                                     <br>
 *  KEYS                                                                               <br>
 *  (key_1) (key_2) ... (key_n)                                                        <br>
 *                                                                                     <br>
 *  NUMPROPS                                                                           <br>
 *  (exact # of properties in the table)                                               <br>
 *                                                                                     <br>
 *  NUMSPECIES                                                                         <br>
 *  (exact # of species in the table)                                                  <br>
 *                                                                                     <br>
 *  MIXF_INDEX                                                                         <br>
 *  (mixture fraction column index in data files)                                      <br>
 *                                                                                     <br>
 *  TEMP_INDEX                                                                         <br>
 *  (temperature column index in data files)                                           <br>
 *                                                                                     <br>
 *  RHO_INDEX                                                                          <br>
 *  (rho column index in data files)                                                   <br>
 *                                                                                     <br>
 *  ENTH_INDEX                                                                         <br>
 *  (enthalpy column index in data files)                                              <br>
 *                                                                                     <br>
 *  LAMBDA_INDEX                                                                       <br>
 *  (lambda column index in data files)                                                <br>
 *                                                                                     <br>
 *  VISC_INDEX                                                                         <br>
 *  (viscosity column index in data files)                                             <br>
 *                                                                                     <br>
 *  HSENS_INDEX                                                                        <br>
 *  (sensible enthalpy column index in data files)                                     <br>
 *                                                                                     <br>
 *  KABS_INDEX                                                                         <br>
 *  (k_abs column index in data files)                                                 <br>
 *                                                                                     <br>
 *  DMIXF_INDEX                                                                        <br>
 *  (mixture fraction diffusivity column index in data files)                          <br>
 *                                                                                     <br>
 *  CP_INDEX                                                                           <br>
 *  (cp column index in data files)                                                    <br>
 *                                                                                     <br>
 *  TABLETYPE                                                                          <br>
 *  (&nbsp;For an Equilibrium table: omit tag key or put
 *   '0' without quotations or 'eta' lowercase and without quotation marks             <br>
 *   &nbsp;&nbsp;For a&nbsp;&nbsp;Flamelet table:&nbsp;&nbsp;&nbsp;
 *   '1' without quotations or 'flamelet' lowercase and without quotation marks&nbsp;) <br>
 *                                                                                     <br>
 *  SPECIES_MASS_INDEX                                                                 <br>
 *  (column index of the first species mass fraction)                                  <br>
 *                                                                                     <br>
 *  <VAR>Note: species mass fractions MUST be sequential</VAR>                         <br>
 *                                                                                     <br>
 *  SPECIES_MOLE_INDEX                                                                 <br>
 *  (column index of the first species mole fraction)                                  <br>
 *                                                                                     <br>
 *  <VAR>Note: species mole fractions MUST be sequential</VAR>                         <br>
 *                                                                                     <br>
 *  </DFN>                                                                             <br>
 *  The '#' symbol can be used within the configuration file
 *  for comments.                                                                      <br>
 *                                                                                     <br>
 *  The USER calls the constructor TableReader::TableReader(string, odtline*)
 *  and then, to retrieve the table, must
 *  call TableReader::getETA() which returns a pointer to an ETA_Table object.                  <br>
 *                                                                                     <br>
 *  <STRONG>Limitations:</STRONG>                                                      <br>
 *  <UL>
 *  <LI>The configuration file MUST be '/odt/input/eqProf.conf'</LI>
 *  <LI>Multiple files must be divided by the property of heat loss.</LI>                   <br>
 *  </UL>
 *  @author Ryan Hintze
 */

class TableReader {

    ////////////////////// CONSTANTS/ENUMS /////////////////////
public:
    //! The type of table being used. Defined in odtParam
    enum tableType {
                     /// Equilibrium Table
                     TT_EQUILIBRIUM,
                     /// Flamelet Table
                     TT_FLAMELET     };

	////////////////////// DATA MEMBERS /////////////////////

public:

    table::prop_indices index;

    vector<vector<vector<double> > >   Table;                  ///< All info stored here. Indexed by [heatloss][mixf][Value]
    vector<double>                     hlosses;                ///< A vector containing the heat loss values
	vector<double>                     hsens;                  ///< A vector containing the hsens corresponding to the mixf values
    vector<double>                     mixf;                   ///< A vector containing the values of the mixf for beginning
	vector<double>					   enth;                   ///< A vector containing the values of enthalpy based on mixf vector    
	vector<string>                     keys;                   ///< A vector containing key values for each column such as 'Temp' or 'H2'

    int        numColumns;             ///< A counter for the number of Columns in Data Set
    int        numRows;                ///< A counter for the number of Rows in Data set
	int        fileCount;              ///< The number of files read into the table
	int        numprops;               ///< The number of columns that are properties
    int        numspecies;             ///< The number of columns that are species
	int        mixf_index;             ///< The index to the column that contains mixture fraction
	int        temp_index;             ///< The index to the column that contains temperature
	int        rho_index;              ///< The index to the column that contains rho
	int        enth_index;             ///< The index to the column that contains enthalpy
	int        lambda_index;           ///< The index to the column that contains lambda
	int        visc_index;             ///< The index to the column that contains viscosity
	int        hsens_index;            ///< The index to the column that contains sensible enthalpy
	int        kabs_index;             ///< The index to the column that contains K_abs
	int        Dmixf_index;            ///< The index to the column that contains Dmixf
	int        cp_index;               ///< The index to the column that contains heat capacity
	int        species_mass_index;     ///< The index to the column where data containing species mass starts
    int        species_mole_index;     ///< The index to the column where data containing species mole starts
    bool       keysSet;                ///< True if the keys have been set, False otherwise
	bool       hsensFile;              ///< Specifies a file where hsens information is
	int        keyCount;               ///< The number of keys contained in the file. Should be equal to numColumns
    tableType  type;                   ///< The Type of Table that is being read.
	odtline*   odtl;                   ///< A pointer to the odtline object


private:
    //! A simple structure that relates a heat loss to a file. The
    //! '<' operator is overloaded for sorting.
    struct indexedFile {
        string fileName;                                ///< File name
        double hloss;                                   ///< Heat loss
        bool operator < (const indexedFile& j) const {return (hloss < j.hloss);}   ///< Overloaded for sorting by heat loss value
    };

    ////////////////////// MEMBER FUNCTIONS  /////////////////////
	vector<indexedFile>                files;                  /// A vector containing the parsed file Names

	void init();
    string         getNextLine(ifstream&);
	vector<string> getHeaderInfo(ifstream&);
	int countColumns(string in, string delimiter);

    double getMixf(string in);
    double getEnth(string in);
    double getHSens(string s);
    double getTemp(string s);
    double getColumn(string s, int col);
    vector<double> getRowData(string s);
    void readConfFile(string confFile);
    void readFiles();
    void addColumns();
    string getNextTag(ifstream&);
    void handleTag(string, ifstream&);
    /////// HandleTag Helper Functions ///////
    void handleTag_FILES(ifstream&);
    void handleTag_KEYS(ifstream&);
    void handleTag_HSENS_FILE(ifstream&);
    indexedFile    getIndexedFile(string);
    //////////////////////////////////////////

    void           tokenizer(string, vector<string>&);

public:

	////////////////////// CONSTRUCTOR FUNCTIONS  /////////////////////


	TableReader(string dir, odtline * odt);
        //default constructor
        TableReader();

};

#endif

