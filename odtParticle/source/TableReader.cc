/**
 * @file TableReader.cc
 * Implementation of TableReader.h
 */

///////////////////////////////////////////////////////////////////////////////

/** Class for parsing output of equilTable data files.
 *  
 *  @author David O. Lignell
 */

////////////////////// INCLUDE ////////////////////////////

#include "TableReader.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstdlib>
using namespace std;

////////////////////////////////////////////////
/** Default Constructor. Does not initialize anything.
 */
TableReader::TableReader(){
}

/** Constructor function to use. Automatically loads and sets up the table based
 * on the configuration file.
 *
 * @param confFile \input path to the configuration file that helps set up the table
 * @param _odtl    \input odtline object, set pointer with
 */
TableReader::TableReader(string confFile, odtline * _odtl) {
    odtl = _odtl;
    init();
    readConfFile(confFile);
    readFiles();
    if(odtl->odtP->ItableLookup == odtl->odtP->TT_FLAMELET)
        addColumns();
    /*cout << files[0].fileName << files[1].fileName << files[2].fileName << files[3].fileName <<
            files[4].fileName << files[5].fileName << files[6].fileName << files[7].fileName <<
            files[8].fileName << files[9].fileName << files[10].fileName << files[11].fileName <<
            files[12].fileName << files[13].fileName << files[14].fileName << endl;*/

    index.Dmixf_index   = Dmixf_index;
    index.enth_index    = enth_index;
    index.hsens_index   = hsens_index;
    index.kabs_index    = kabs_index;
    index.lambda_index  = lambda_index;
    index.mixf_index    = mixf_index;
    index.numprops      = numprops;
    index.rho_index     = rho_index;
    index.temp_index    = temp_index;
    index.visc_index    = visc_index;
    index.species_mass_index = species_mass_index;
    index.species_mole_index = species_mole_index;

//    if(odtl->odtP->Lprxn)
        index.cp_index      = cp_index;
}

///////////////////////////////////////////////////////////////////////////////
/** Initializes the TableReader to default values
 *
 */

void TableReader::init(){

    numRows      =  0;
    numColumns   =  0;
    fileCount    =  0;
    numprops     =  0;
    numspecies   =  0;
    mixf_index   = -1;
    temp_index   = -1;
    rho_index    = -1;
    enth_index   = -1;
    lambda_index = -1;
    visc_index   = -1;
    hsens_index  = -1;
    kabs_index   = -1;
    Dmixf_index  = -1;
    if(odtl->odtP->Lprxn)
        cp_index       = -1;
    species_mass_index = -1;
    species_mole_index = -1;
    keysSet      = false;
    hsensFile    = false;
    if(odtl->odtP->ItableLookup == odtl->odtP->TT_EQUILIBRIUM)
        type         = TT_EQUILIBRIUM;
    else if(odtl->odtP->ItableLookup == odtl->odtP->TT_FLAMELET)
        type         = TT_FLAMELET;
}


///////////////////////////////////////////////////////////////////////////////

/** Reads a given configuration file to read in the table properly.
 *  This calls handleTag(string, ifstream&) and may change any
 *      of the variables listed under its description.
 * @param confFileName \input Path and name of the configuration file to be
 *      used.
 * @returns void
 */
void TableReader::readConfFile(string confFileName) {

    // Open File
    ifstream confFile;                              //File object for configuration file
    confFile.open(confFileName.c_str(), ios::in);

    //Make sure file is opened
    if (!confFile.is_open())
        perror("Error opening File");

    //Read in each tag and parse it
    while(!confFile.eof())
    {
        string nextTag = getNextTag(confFile);      // String containing the next tag Name
        if(nextTag.length() > 0)
        {
            //cout << "Returned \'" << nextTag << "\'" << endl;
            handleTag(nextTag, confFile);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////

/** Parses each tag appropriately.
 * May change the following TableReader variables: <br>
 * <ul>
 * <li>numprops</li>             <li>numspecies</li>          <li>index_mixf</li>
 * <li>index_temp</li>           <li>index_rho</li>           <li>index_enth</li>
 * <li>index_lambda</li>         <li>index_visc</li>          <li>index_kabs</li>
 * <li>index_Dmixf</li>          <li>index_cp</li>            <li>index_species_mass</li>
 * <li>index_species_mole</li>   <li>type</li>
 * <li> Any of the variables listed under <ul>
 *     <li>handleTag_FILES</li><li>handleTag_KEYS</li>
 *     <li>handleTag_HSENS_FILE</li>
 * </li>
 * </ul>
 * </ul>
 *
 *  @param tag        \input string - the tag to be parsed
 *  @param confFile   \input ifstream - the file being parsed
 */
void TableReader::handleTag(string tag, ifstream& confFile){
    if (tag.compare("FILES") == 0) {
        handleTag_FILES(confFile);
    } else if(tag.compare("KEYS")          == 0){
        handleTag_KEYS(confFile);
    } else if(tag.compare("HSENS_FILE")    == 0){
        handleTag_HSENS_FILE(confFile);
    } else if(tag.compare("NUMPROPS")      == 0){
        numprops      = atof(getNextTag(confFile).c_str());
    } else if(tag.compare("NUMSPECIES")    == 0){
        numspecies    = atof(getNextTag(confFile).c_str());
    } else if(tag.compare("MIXF_INDEX")    == 0){
        mixf_index    = atof(getNextTag(confFile).c_str());
    } else if(tag.compare("TEMP_INDEX")    == 0){
        temp_index    = atof(getNextTag(confFile).c_str());
    } else if(tag.compare("RHO_INDEX")     == 0){
        rho_index     = atof(getNextTag(confFile).c_str());
    } else if(tag.compare("ENTH_INDEX")    == 0){
        enth_index    = atof(getNextTag(confFile).c_str());
    } else if(tag.compare("LAMBDA_INDEX")  == 0){
        lambda_index  = atof(getNextTag(confFile).c_str());
    } else if(tag.compare("VISC_INDEX")    == 0){
        visc_index    = atof(getNextTag(confFile).c_str());
    } else if(tag.compare("HSENS_INDEX")   == 0){
        hsens_index   = atof(getNextTag(confFile).c_str());
    } else if(tag.compare("KABS_INDEX")    == 0){
        kabs_index    = atof(getNextTag(confFile).c_str());
    } else if(tag.compare("DMIXF_INDEX")   == 0){
        Dmixf_index   = atof(getNextTag(confFile).c_str());
    } else if(tag.compare("CP_INDEX")      == 0){
        cp_index      = atof(getNextTag(confFile).c_str());
    } else if(tag.compare("TABLETYPE")     == 0){
        string type_str = getNextTag(confFile);
        if(type_str.compare("0") || type_str.compare("eta")      == 0)
            type = TT_EQUILIBRIUM;
        if(type_str.compare("1") || type_str.compare("flamelet") == 0)
            type = TT_FLAMELET;
    } else if(tag.compare("SPECIES_MASS_INDEX") == 0){
        species_mass_index = atof(getNextTag(confFile).c_str());
    } else if(tag.compare("SPECIES_MOLE_INDEX") == 0){
        species_mole_index = atof(getNextTag(confFile).c_str());
    } else {
        cout << "Invalid Tag - \'" << tag << "\'" << endl;
        perror("Check configuration file");
    }

}

///////////////////////////////////////////////////////////////////////////////

/** Parses a HSENS_FILE tag and its contents.
 *
 *  The hsens file is used to include a column that contains the corresponding hsens
 *  data into the table.
 *  @param confFile \input the ifstream from which the HSENS_FILE tag was found
 *          (the configuration file's ifstream).
 */

void TableReader::handleTag_HSENS_FILE(ifstream& confFile){
    //Get the file name and info
    string hsi = getNextTag(confFile);
    vector<string> fi;
    tokenizer(hsi, fi);

    // Open File
    ifstream hsensStream;                              //File object for hsens file
    hsensStream.open(fi[0].c_str(), ios::in);
    int col = atoi(fi[1].c_str());                   //Column that contains hsens info
    //Make sure file is opened
    if (!hsensStream.is_open())
       perror("Error opening File");

    //Read in each tag and parse it
    numColumns = col+1;
    while(!hsensStream.eof())
    {
       string nextLine = getNextTag(hsensStream);      // String containing the next row
       if(nextLine.length() > 0)
       {
           //cout << nextLine << endl;
           hsens.push_back(getColumn(nextLine,col));
           //cout << hsens[hsens.size()-1] << endl;

       }
    }
    numColumns = 0;
    hsensStream.close();


    hsensFile = true;

}

///////////////////////////////////////////////////////////////////////////////

/** Parses a KEYS tag and its contents. This only parses one line after the
 * KEYS tag and all keys must be on that line.<br>
 * <br>
 * ex:<br>
 * <br><DFN>
 * KEYS<br>
 * mixf temp ...</DFN>
 *
 *  @param confFile  \input the ifstream from which the KEYS tag was found
 *          (the configuration file's ifstream).
 */
void TableReader::handleTag_KEYS(ifstream& confFile){

    //Get the keys
    string ks = getNextTag(confFile);
    tokenizer(ks, keys);
    keysSet = true;
    keyCount = keys.size();

}

///////////////////////////////////////////////////////////////////////////////
/** Tokenizes the string and adds them to the given vector.
*       Tokens are added in order and they MUST be delimited by
 *      the space character.
 *
 *  @param str \input string to be tokenized
 *  @param vect \output vector to which tokens are to be added
 */

void TableReader::tokenizer(string str, vector<string>& vect) {

    string delimiters = " ";
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        vect.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

///////////////////////////////////////////////////////////////////////////////
/** Parses a FILES tag and its contents from a configuration file
 *  stream. Increments fileCount and adds the stream to the files
 *  vector.
 *
 *  @param confFile \input the ifstream from which the FILES tag was found
 *          (the configuration file's ifstream).
 */
void TableReader::handleTag_FILES(ifstream& confFile){
    string s = getNextTag(confFile);
    fileCount = 0;
    while(s.compare("END") != 0){
        fileCount++;
        files.push_back(getIndexedFile(s));
        s = getNextTag(confFile);
    }
    sort(files.begin(),files.end());
}

///////////////////////////////////////////////////////////////////////////////
/** Creates an indexed file from a line read in by TableReader::handleTag_FILES
 *
 *  This method collects two values from a single line: filename and heat loss.
 *
 *  The line must start with the filename (no preceding whitespace)
 *  and the heat loss must appear on the same line, separated by whitespace
 *  that does not end the line.<br>
 *  <br>
 *  ex: <br>
 *  <br><DFN>
 *  myfile.dat 0.05
 *  </DFN>
 *
 *  @param fileName \input string containing filename and heat loss value.
 *  @return indexedFile from parsing
 */

TableReader::indexedFile TableReader::getIndexedFile(string fileName){
    indexedFile ret;
    //Find File Name
    int sI = 0;                         //Start Index
    int eI = fileName.find_first_of(" \t"); //End Index
    ret.fileName = fileName.substr(sI, eI-sI);
    //cout << "fileName - " << ret.fileName << endl;

    //Hloss value
    sI = fileName.find_first_not_of(" \t",eI);
    eI = fileName.length();
    ret.hloss = atof(fileName.substr(sI, eI-sI).c_str());
    //cout << "hloss - " << ret.hloss << endl;
    return ret;
}

///////////////////////////////////////////////////////////////////////////////
/** Searches the given stream for the next non-comment tag, extra whitespace removed.
 *
 * @param file \input ifstream from which to find the next non-comment line. This should
 *       be the stream of the configuration file.
 */

string TableReader::getNextTag(ifstream& file)
{
    string s;
    int nonWSIndex = -1;                    // Non-WhiteSpace Index
    while(!file.eof())
    {
        s = getNextLine(file);
        //cout << "getNextTag - \'" << s << "\'" << endl;
        //Make sure that the line is not blank
        nonWSIndex = s.find_first_not_of(" \t\r\n");
        if(nonWSIndex >= 0)
            return s.substr(nonWSIndex, (s.find_last_not_of(" \t\r\n") - nonWSIndex) + 1);
    }
    return "";
}


///////////////////////////////////////////////////////////////////////////////

/*
 * @param dir input, string, Tells which directory to search for eqProfile
 * @param odtln input, odtline object, set pointer with
 */
/*TableReader::TableReader(std::string dir, odtline *odtln) {

    cout << "\n# ************* Parsing Directory \'" << dir << "\'\n" << endl;

    odtl           =       odtln;				// Pointer to odtline
    numColumns     =       0;					// Column Counter
    numRows        =       0;					// Row Counter

	////////////////////////////// Find Files in Directory ///////////////////////////////////////////
	ifstream       files[50];					        // Array of file to be searched
	string         fNames[50];					        // Files' names, set after files are opened
	int            numFiles    =  0;	                // Counter for number of files

	// Create directory entry
	struct dirent  *entry;								// A single entry inside a directory
	DIR*           dp          =  opendir(dir.c_str());	// Directory reader

	// Make sure directory is found
	if (dp == NULL)
		perror("Opening Directory");

	// Found! search for eqProfile files by searching all entries
	while ((entry = readdir(dp))) {

		string fName = entry->d_name;					        //A single entry's file name

		if (fName.find("eqProfile_") == 0) {

			//eqProfile file found

			//Get the eqProfile hloss index from the filename
			int fI = getFileIndex(fName);				//Index of eqProfile given by the fileName

			//Open into file array
			cout << "# opening " << dir  + fName << endl;
			files[fI].open((dir + "/" + fName).c_str(), ios::in);

			//Make sure file is opened
			if (!files[fI].is_open())
				perror("Error opening File");

			//Success
			fNames[fI] = fName;
			numFiles++;
		}

	}
	closedir(dp);
 */
 
/** Reads in the data from all the files found in the configuration files.
 * This checks that all files have the same number of columns and rows
 * for consistency.
 *
 */
void TableReader::readFiles()
{
	////////////////////////////// READ FILES INTO THE TABLE ///////////////////////////////////////////
    // Open Files
	ifstream*    inFiles = new ifstream[fileCount];
    for(int i =0; i< fileCount; i++) {
        inFiles[i].open(files[i].fileName.c_str(), ios::in);

        //Make sure file is opened
        if(!inFiles[i].is_open()) {
            cout << "Could not find \'" << files[i].fileName << "\'" << endl;
            perror("Error Opening File");
        }
    }

	// Files are Open, Read them into a table
	bool columnsNumSet = false;                         //For Verifying column consistency
	bool rowsNumSet    = false;                         //For Verifying row consistency

	hlosses.resize(fileCount);
	if(!keysSet)
	    keys = getHeaderInfo(inFiles[0]);

	//Do each file individually
	for (int i = 0; i < fileCount; i++) {

		//First Get the Heat Loss Value
		hlosses[i] = files[i].hloss;				//File's heat loss value

		//Next Find how many columns there are in the file
		string s;
		s = getNextLine(inFiles[i]);
		int colCount = countColumns(s, " ");

		//Testing Files for consistency
		if (!columnsNumSet) {

			columnsNumSet = true;
			numColumns = colCount;

		} else {

			//Test if they are equal
			if (numColumns != colCount) {

				//If not throw an error.
				perror("Differing Column Lengths in Files");

			}
		}

		//Start reading in data (starting with s)
		vector<vector<double> > fData;                    //Data from the file in 2D array

        //We need to include the already read row
		if (!rowsNumSet) {
			numRows++;

			//Grab mixf and corresponding hsens
			mixf.push_back(getMixf(s));
			if(hsens_index >=0)
			    hsens.push_back(getHSens(s));
			if(enth_index >=0)
				enth.push_back(getEnth(s));
		}

		//Add row to file data
		fData.push_back(getRowData(s));

		//Read in the rest of the file
		while (!inFiles[i].eof()) {
			s = getNextLine(inFiles[i]);
			if(s.length() > 0) {

                //Count number of Rows if not set
                if (!rowsNumSet) {
                    numRows++;
                    //Grab mixf and corresponding hsens (Only put into the Vector once)
                    mixf.push_back(getMixf(s));
                    if(hsens_index >=0)
                        hsens.push_back(getHSens(s));
					if(enth_index >=0)
						enth.push_back(getEnth(s));
                }

                //Finally add the data from the row
                fData.push_back(getRowData(s));
			}
		}
		//number of rows set after first
		rowsNumSet = true;
		//Put data into the table
		Table.push_back(fData);

	}
	//cout << "# Done" <<endl;
	delete[] inFiles;
}

///////////////////////////////////////////////////////////////////////////////
/** Add data not included in Table
 *
 *  @return void
 */
void TableReader::addColumns(){

    int numToAdd = 5 + numspecies;          //number of columns to all
    vector<double>          molysp(numspecies, 0.0);
    vector<double>          rho(numRows,0.0);
    vector<double>          visc(numRows,0.0);
    vector<double>          lambda(numRows,0.0);
    vector<double>          kabs(numRows,0.0);
    vector<double>          Dmixf(numRows,0.0);
    vector<vector<double> > radCoefs;  ///< [spc][coef]
    vector<int>    iRadIndx;  ///< radiation species indicies: ch4 co2 h2o co: negative if not present

    radCoefs.resize(4,vector<double>(6,0)); // kp (=) 1/atm*m
    radCoefs[0][0] =  1.017015E+1;          // ch4; kp=2.798 at 1150K
    radCoefs[0][1] = -7.947312E-03;
    radCoefs[0][2] =  4.342446E-7;
    radCoefs[0][3] =  1.048611E-9;
    radCoefs[0][4] = -2.287861E-13;
    radCoefs[0][5] =  0.000000E+0;
    radCoefs[1][0] =  3.24442E+1;           // co2; kp=29.197 at 925K
    radCoefs[1][1] =  7.537513E-02;
    radCoefs[1][2] = -1.535140E-04;
    radCoefs[1][3] =  9.48794E-8;
    radCoefs[1][4] = -2.509259E-11;
    radCoefs[1][5] =  2.447995E-15;
    radCoefs[2][0] =  6.86948E+1;           // h2o; kp=4.474  at 1119K
    radCoefs[2][1] = -1.52349E-01;
    radCoefs[2][2] =  1.417848E-04;
    radCoefs[2][3] = -6.620996E-8;
    radCoefs[2][4] =  1.52415E-11;
    radCoefs[2][5] = -1.373456E-15;
    radCoefs[3][0] =  1.56536E+0;           // co; kp=2.501 at 1007 K
    radCoefs[3][1] =  1.483914E-02;
    radCoefs[3][2] = -2.656035E-05;
    radCoefs[3][3] =  1.68798E-8;
    radCoefs[3][4] = -4.674473E-12;

   iRadIndx.resize(4);
   bool fmissing = false;
   int isp;
   int lastIndex = keys.size() - 1;         //The last key's index

   isp = -1;
   for (int i = 0; i <= lastIndex; i++)
       if(keys[i].compare("CH4") == 0)
           isp = i;
   iRadIndx[0] = isp;
   if(isp < 0) fmissing = true;

   isp = -1;
   for (int i = 0; i <= lastIndex; i++)
       if(keys[i].compare("CO2") == 0)
           isp = i;
   iRadIndx[1] = isp;
   if(isp < 0) fmissing = true;

   isp = -1;
   for (int i = 0; i <= lastIndex; i++)
       if(keys[i].compare("H2O") == 0)
           isp = i;
   iRadIndx[2] = isp;
   if(isp < 0) fmissing = true;

   isp = -1;
   for (int i = 0; i <= lastIndex; i++)
       if(keys[i].compare("CO") == 0)
           isp = i;
   iRadIndx[3] = isp;
   if(isp < 0) fmissing = true;
   if(fmissing) 
       cout << endl << "Warning one or more radiating species missing from mechanism" << endl;

    for(int k = 0; k < fileCount; k++){
        for(int i = 0; i < numRows; i++){
            vector<double> yspc(numspecies,0.0);
            for(int j = 0; j < numspecies; j++){
                yspc[j] = Table[k][i][species_mass_index+j];
            }
            odtl->gas->setState_TPY(Table[k][i][temp_index], odtl->pres, &yspc[0]);
            lambda[i] = odtl->tran->thermalConductivity();
            visc[i] = odtl->tran->viscosity();
            rho[i] = odtl->gas->density();
            odtl->gas->getMoleFractions( &molysp[0]);
            Dmixf[i] = lambda[i]/rho[i]/odtl->gas->cp_mass();

            double Kabs;
            double Ktot = 0.0;

            for(int n=0; n<4; n++) {
                if(iRadIndx[n] < 0)         // this radiation species k is not in the mechanism
                    continue; 
                Kabs = radCoefs[n][5];
            for(int j=4; j>=0; j--)
                Kabs = Kabs * Table[k][i][temp_index] + radCoefs[n][j];

            Ktot += molysp[iRadIndx[n]]*odtl->pres/101325.*Kabs;
            }
            kabs[i] = Ktot;

            Table[k][i].resize(numColumns+numToAdd);
            Table[k][i][numColumns] = rho[i];
            rho_index = numColumns;
            Table[k][i][numColumns + 1] = lambda[i];
            lambda_index = numColumns+1;
            Table[k][i][numColumns + 2] = visc[i];
            visc_index = numColumns + 2;
            Table[k][i][numColumns + 3] = Dmixf[i];
            Dmixf_index = numColumns + 3;
            Table[k][i][numColumns + 4] = kabs[i];
            kabs_index = numColumns + 4;
            for(int j = 0; j < numspecies; j++)
                Table[k][i][numColumns + 5 + j];
            species_mole_index = numColumns + 5;

        }
    }
    numColumns += numToAdd;
}

///////////////////////////////////////////////////////////////////////////////

/** Parse out the mixture fraction from the given row
 *  The variable index_mixf must be set (most likely by readConfFile) or this will
 *  throw a fatal error.
 *
 *  @param  in \input string containing an unparsed line from a data file.
 *  @return the mixture fraction for the given row
 */

double TableReader::getMixf(string in) {
	if(mixf_index <0) {
	    perror("MIXF_INDEX IS NOT DEFINED, PLEASE DEFINE IN CONFIGURATION FILE");
	    exit(1);
	}
	return getColumn(in, mixf_index);
}
////////////////////////////////////////////////////////////////////////////
/** Parse out the enthalpy from the given row
 *
 *  @param in input: string containing an unparsed line from a eqProfile file.
 *  @return the Enthalpy for this row
 */

double TableReader::getEnth(string in) {
	if(enth_index <0) {
	    perror("ENTH_INDEX IS NOT DEFINED, PLEASE DEFINE IN CONFIGURATION FILE");
	    exit(1);
	}
	return getColumn(in, enth_index);
}

///////////////////////////////////////////////////////////////////////////////

/** Parse out sensible enthalpy (hsens) from the given row.
 *  The variable index_hsens must be set (most likely by readConfFile) or this will
 *  throw a fatal error.
 *
 *  @param  s \input string containing an unparsed line from a data file.
 *  @return the HSens value for the given row
 */

double TableReader::getHSens(string s) {
    if(hsens_index <0) {
        perror("HSENS_INDEX IS NOT DEFINED, PLEASE DEFINE IN CONFIGURATION FILE");
        exit(1);
    }
	return getColumn(s, hsens_index);
}

///////////////////////////////////////////////////////////////////////////////

/** Parse out the Temperature from the given row.
 **  The variable index_temp must be set (most likely by readConfFile) or this will
 *  throw a fatal error.
 *  @param s \input string containing an unparsed line from a data file.
 *  @return the temperature value for the given row
 */

double TableReader::getTemp(string s) {
    if(temp_index <0) {
        perror("TEMP_INDEX IS NOT DEFINED, PLEASE DEFINE IN CONFIGURATION FILE");
        exit(1);
    }
	return getColumn(s, temp_index);
}

///////////////////////////////////////////////////////////////////////////////

/* Get's the HLVal from a standard eqProfile name
 *
 *  The file's name should be formatted like 'eqProfile_1_-0.15.dat'. Here the 1
 *  corresponds to the file index and the -0.15 is the heat loss value for the data
 *  within the file.
 *
 *  @param fName input: a string containing a name of a eqProfile file.
 *  @return Retrieved heat loss value.
 */

/*double TableReader::getHLVal(string fName) {

    int sIndex;         //Indexes to the start of the heatloss value
    string num;         //contains the hloss after parsing
    sIndex = fName.find('_', 10) + 1;
    num = fName.substr(sIndex, fName.find(".dat", 10) - sIndex);
    return atof(num.c_str());
}*/

///////////////////////////////////////////////////////////////////////////////

/* Get's the file's index from a standard eqProfile name
 *
 *  The file's name should be formatted like 'eqProfile_1_-0.15.dat'. Here the 1
 *  corresponds to the file index and the -0.15 is the heat loss value for the data
 *  within the file.
 *
 *  @param fName input: a string containing a name of a eqProfile file.
 *  @return Retrieved file index.
 *

int TableReader::getFileIndex(string fName) {
    string num = fName.substr(10, fName.find('_', 10) - 10);
    return atoi(num.c_str());
}*/

///////////////////////////////////////////////////////////////////////////////

/** Parses a string to find the value of a certain column
 *
 *  @param s   \input string containing an unparsed line from a eqProfile file.
 *  @param col \input column to look up
 *  @return the value of the column col in string s as a double
 */

double TableReader::getColumn(string s, int col) {

    //Make sure the column is not out of bounds
	if (col >= numColumns)
		perror("Invalid Format, check columns");

	int start = 0;
	int next  = 0;

	//Finding Column but checking for non-whitespace characters
	for (int i = 0; i < col + 1; i++) {
		start = s.find_first_not_of(' ', next);
		next = s.find(" ", start);
	}

	if (start < 0)
		perror("Problem Parsing File");

	string sub = s.substr(start, next - start);

	return atof(sub.c_str());

}

///////////////////////////////////////////////////////////////////////////////

/** Parses a row of data to return a vector of doubles that contains
 *  the data contained in the line. Data is assumed to be separated
 *  by spaces.
 *
 *  @param s \input A string containg a row data separated by spaces
 *  @return a vector of doubles containing the parsed row data
 */

vector<double> TableReader::getRowData(string s) {
	vector<double> rowData(numColumns);    // A vector for storing the row data
	int            start = 0;              // Indexes to the start of the next data point in the string
	int            next  = 0;              // Indexes to the end of the data point after start
	//Add all values of y
	for (int i = 0; i < numColumns; i++) {
		start = s.find_first_not_of(' ', next);
		next = s.find(" ", start);

		if (start < 0)
			perror("Problem Parsing File");

		string sub = s.substr(start, next - start);
		//cout << sub << endl;
		rowData[i] = atof(sub.c_str());

	}
	return rowData;

}

///////////////////////////////////////////////////////////////////////////////

/** A function that counts the number of columns separated by a delimiter
 *
 *  @param line  \input string to be counted
 *  @param delimiter \input string containing the delimiter that separates the columns
 *  @return The number of columns separated by the delimiter
 */

int TableReader::countColumns(string line, string delimiter) {
	int count = 0;                                 // current column Count
	int start = line.find_first_not_of(delimiter, 0);      // indexes to first value not of the delimiter
	int next  = 0;                                 // indexes to the next occurrence of the delimiter
	if (start < 0)
		perror("Invalid File Format");
	while (start >= 0) {
		next = line.find(delimiter, start);
		start = line.find_first_not_of(delimiter, next);
		count++;
	}
	//cout << "Col. Count - " << count << endl;
	return count;

}

///////////////////////////////////////////////////////////////////////////////

/** Gets the next line that is not a comment (doesn't start with #)
 *  and is not blank (not just spaces)
 *
 *  @param file \input an input file stream from which to retrieve the next line
 *  @return the next non-comment line from the ifstream
 */

string TableReader::getNextLine(ifstream& file) {
	string s = "#";                    // Stores the next line

	//While not the end of the file and not a comment or blank line
	while ((s.length() == 0 || s.at(0) == '#') && s.find_first_not_of(" \n\r\t") >= 0
	        && !file.eof()) {

		string nLine;                  // Empty address for reading in the next line
		getline(file, nLine);
		s = nLine;
	}
	if (file.eof() && s.compare("#") == 0)
		perror("End of File Reached Early");
	//cout << "Returning " << s << endl;
	return s;
}

///////////////////////////////////////////////////////////////////////////////

/** Parses an data files header to give key values to each column.
 *  The header should be on the first line and the header must be formatted as:
 *  <br> <DFN>'#      1_mf      2_hsens     ...     n_key'</DFN> <br>
 *  <br>
 *  This function would parse the preceding as:
 *  <br> <DFN>\<'mf','hsens', ... , 'key'\></DFN>
 *
 *  @param file \input a ifstream for a file that contains a relevant header.
 *  This will only read one line of the file.
 *  @return  A vector of the parsed key names.
 */

vector<string> TableReader::getHeaderInfo(ifstream& file) {
    string         s;                       // A String for dumping the row data
    vector<string> info;                    // A vector to store the key values from Header
    int            sIndex;                  // Indexes to the start of a key value in s
    int            eIndex;                  // indexes to the end of a key value in s
    getline(file, s);
    if(s.at(0) != '#'){
        //No Header, Index by number
        for(int i = 0; i < numColumns; i++)
        {
            stringstream ss;
            ss << i;
            info.push_back(ss.str());
        }
    } else {
        //Header, try to index

        //Find first _ which precedes key value
        sIndex = s.find('_', 1)+1;

        //Find space after sIndex which follows the key value
        eIndex = s.find(' ', sIndex);

        while(sIndex > 0)
        {
            //Parse Key and add to vector
            string key = s.substr(sIndex,eIndex-sIndex);
            info.push_back(key);

            //Find next key values
            sIndex = s.find('_', sIndex)+1;
            eIndex = s.find(' ', sIndex);
        }
    }
    return info;
}

