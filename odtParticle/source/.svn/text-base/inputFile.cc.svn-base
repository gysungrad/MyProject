/**
  * @file inputFile.cc
  * Source file for inputFile class
  */

#include "inputFile.h"
#include <vector>
#include <map>

using namespace std;


/////////////////////////////////////////////////////
/** Constructor function */

inputFile::inputFile ()
{
    fileSet_ = false;
    fileParsed_ = false;
}

/////////////////////////////////////////////////////

/** Constructor function 
 *  @param fileName \input input file to open
 */
inputFile::inputFile (string fileName) : filename_(fileName)
{
    fileSet_ = true;
    parseInput();
}



///////////////////////////////////////////////////////

/** @details    This will re-set the input file associated with this inputFile object
  *             and will parse the newly-designated input file.
  * @param fileName \input the new file name to set with.
  */
void inputFile::setFile (string fileName) {
    filename_ = fileName;
    fileSet_ = true;
    fileParsed_ = false;

    parseInput();
}

//////////////////////////////////////////////////////

/** @details  This function can be used to parse an arbitrary number of input files.
  *           The new variables from each input file are simply added to the input_map_
  *           map.  This facilitates the ability to either split variables into groups
  *           to organize input files, or combine all variables into a single input file.
  *           
  *           Vectors are parsed as follows:                                            \n
  *           \vc{
  *           --------------------------------------------------------------------------
  *           1         (someVariableName)  This is the name associated with the vector
  *           v
  *           80        (element1)
  *           82        (element2)
  *           87        (element3)
  *           v
  *           }
  *
  *           The vector is put into input_map_ with the variable name \c someVariableName_vector#
  *           (where # is the element number of that vector); so \c 80 has the variable name
  *           \c someVariableName_vector1, \c 82 has the variable name \c someVariableName_vector2,
  *           and so on.  Vectors may be any arbitrary size.
  *
  *           The names of each vector element are stored in name_map_, with the same key:
  *           \c someVariableName_vector#. This way the name of each element and the value of each
  *           element can be obtained easily (and separately, in case element names are unimportant).
  *
  *           Also note that the line starting with a "v" does not have to contain ONLY "v"; it
  *           is treated as a comment line and everything after the "v" character is ignored.
  *
  */
void inputFile::parseInput () {
  // open an ifstream
  ifstream ifile(filename_.c_str());
  if(!ifile) {
    ostringstream errmsg;
    errmsg << "ERROR: inputFile::parseInput(): could not find the specified file: " << filename_ << endl;
    throw runtime_error( errmsg.str() );
  }

  string varname;
  string vectorname;
  string s1;
  
  int lineNumber = 0;
  int z = 0;

  bool vectorMode   = false;  ///< Flag: does this line belong to a vector
  
  while(!ifile.eof() ) {
    
    getline(ifile, s1);
    ++lineNumber;
 
    // -----------------------------

    // check the line for special prefixes:
    if( s1.empty() || s1.substr(0,1)=="#" || s1.substr(0,2)=="--" ) {
      // ignore blank lines and lines starting with "#" or "-"
      continue;

    } else if (s1.substr(0,1)=="v" ) {
      
      // vector tag (starts w/ "v" )
      if (vectorMode == false) {
        // (re)initialize counter to measure vector size
        z = 0;
      }

      // turn on vector mode flag
      vectorMode    = (!vectorMode);
      
      // save previous varname as the vector name:
      vectorname = varname;

    } else {
      // deal with other/non prefixes
    }

    // ---------------------------------

    // Grab keyword from input file, put it in varname
    const unsigned int pos1 = s1.find("(",0);
    const unsigned int pos2 = s1.find(")",0);

    // check to make sure there are a set of ()'s first! (& in left-to-right order)
    if ( pos1 == s1.size() || pos2 < pos1 ) {
      ostringstream errmsg;
      errmsg << "ERROR: inputFile::parseInput(): could not find matching delimiters '(' and ')' in line "
             << lineNumber <<" of input file.  Indicate comment lines using '-' or '#'." << endl;
      throw runtime_error( errmsg.str() );
    }
    varname = s1.substr(pos1+1,pos2-pos1-1);

    if (vectorMode) {
      if (z > 0) {
        // otherwise you are still on the line with the opening vector tag!
        
        // make the right vector name: varname_vector[1,2,3,...]
        stringstream out;
        out << vectorname << "_vector" << z;
        string varname_temp;
        varname_temp = out.str();

        // push the line number into the map
        input_map_[varname_temp] = lineNumber;

        // if element has a name, push element name into the map
        if ( !varname.empty() ) {
          name_map_[varname_temp] = varname;
        }
      }
      ++z;
    } else {
    // -------------------------
    // now deal with everything else

      //save this information into a map
      input_map_[varname] = lineNumber;
    }

  } //endwhile

  fileParsed_ = true;

  // want something from the map? try:
  // input_map_.find('search_key')->second;

  ifile.close();

}



