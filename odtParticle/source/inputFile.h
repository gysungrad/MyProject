/**
  * @file inputFile.h
  * Header file for inputFile class
  */

#ifndef INPUTFILE_H
#define INPUTFILE_H

#include <string>
#include <vector>
#include <map>
#include <typeinfo>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>

using namespace std;

///////////////////////////////////////////////////////

/** 
  * @class      inputFile
  * 
  * @brief      This class wraps an input file and provides access functions for variables in the file.
  * 
  * @details    Each instance of the inputFile class is associated with a particular input file.
  *             This should be managed using a factory, so that multiple inputFile objects are
  *             not created for a single input file.
  *
  * @author     Charles Reid
  *
  */

///////////////////////////////////////////////////////

class inputFile {

    public:

        /** @brief    Default constructor: instantiates the class without setting the
          *           input file it is assocated with (that's done by setFile).
          */
        inputFile ();

        /** @brief    Grabs the variable name and corresponding line number from the input file
          *           and puts it into a map
          * 
          * @param fileName   The name of the file to parse 
          */
        inputFile (string fileName);



    private:

        /** @brief
          */
        void parseInput();


        string                         filename_;    ///< Name of input file being parsed & stored
        map<string, int>          input_map_;   ///< Map containing variable names and their line numbers
        map<string, string>  name_map_;    ///< Map containing names of individual vector elements
        bool                                fileSet_;     ///< Boolean: is there an input file associated with this inputFile object?
        bool                                fileParsed_;  ///< Boolean: is the input file parsed?



    public:

        /** @brief    This function re-sets the input file this inputFile object is associated with
          *
          * @param fileName   The filename of the new input file
          */
        void setFile (string fileName);

        /** @brief      Requests the variable value at the given line number and puts it into double* variable,
          *             using a default value if the variable is unavailable. 
          *
          * @param variableName     The name of the variable given in the input file (between parentheses) 
          * @param variable         The variable in which to put the value given in the input file
          * @param defaultValue     If the variable is not found in the input file, this is the value assigned it
          */ 
        template<class T>
        void getParameter (string variableName, T* variable, T defaultValue) {
          if ( !fileSet_ || !fileParsed_ ) {
            ostringstream errmsg;
            errmsg  << "WARNING: inputFile::getParameter(): the input file has not been set " 
                    << "or has not been parsed. " << endl;
            throw runtime_error( errmsg.str() );
          }
 
          // if a variable cannot be found in the map, use default
          map<string, int>::iterator iLine = input_map_.find(variableName);
          if (iLine == input_map_.end() ) {
            (*variable) = defaultValue;
            return;
          }
          int lineNumber = iLine->second;

          // skip forward to line of interest
          ifstream ifile( filename_.c_str() );
          string s1;
          for (int i = 1; i < lineNumber; ++i) {
            getline(ifile, s1);
          }

          ifile >> (*variable);   getline(ifile, s1);

          ifile.close();
        }



        /** @brief      Requests the variable value at the given line number and puts it into T* variable,
          *             printing a warning if the variable is unavailable. 
          *
          * @param variableName     The name of the variable given in the input file (between parentheses) 
          * @param variable         The variable in which to put the value given in the input file
          */ 
        template<class T>
        void getParameter (string variableName, T* variable) {
          if ( !fileSet_ || !fileParsed_ ) {
            ostringstream errmsg;
            errmsg  << "WARNING: inputFile::getParameter(): the input file has not been set " 
                    << "or has not been parsed. " << endl;
            throw runtime_error( errmsg.str() );
          }
 
          // if a variable cannot be found in the map, there IS no default... so print a warning
          map<string, int>::iterator iLine = input_map_.find(variableName);
          if (iLine == input_map_.end() ) {
            cout  << "WARNING: inputFile::getParameter(): could not find requested variable " << variableName
                       << " in the list of values from the input file." << endl;
            return;
          }
          int lineNumber = iLine->second;

          // skip forward to line of interest
          ifstream ifile( filename_.c_str() );
          string s1;
          for (int i = 1; i < lineNumber; ++i) {
            getline(ifile, s1);
          }

          ifile >> (*variable);   getline(ifile, s1);

          ifile.close();
        }



        /** @brief      Requests the variable value at the given line number and puts it into T* variable,
          *             throwing an error if the variable is unavailable. 
          *
          * @param variableName     The name of the variable given in the input file (between parentheses) 
          * @param variable         The variable in which to put the value given in the input file
          */ 
        template<class T>
        void requireParameter (string variableName, T* variable) {
          if ( !fileSet_ || !fileParsed_ ) {
            ostringstream errmsg;
            errmsg  << "WARNING: inputFile::requireParameter(): the input file has not been set " 
                    << "or has not been parsed. " << endl;
            throw runtime_error( errmsg.str() );
          }
 
          // if a variable cannot be found in the map, there IS no default... so print a warning
          map<string, int>::iterator iLine = input_map_.find(variableName);
          if (iLine == input_map_.end() ) {
            ostringstream errmsg;
            errmsg  << "ERROR: inputFile::requireParameter(): could not find the requested variable " << variableName
                    << " in the list of values from the input file." << endl;
            throw runtime_error( errmsg.str() );
          }
          int lineNumber = iLine->second;

          // skip forward to line of interest
          ifstream ifile( filename_.c_str() );
          string s1;
          for (int i = 1; i < lineNumber; ++i) {
            getline(ifile, s1);
          }

          ifile >> (*variable);   getline(ifile, s1);

          ifile.close();
        }



        /** @brief      This grabs a vector from the input file (more precisely, from the list/map of input file variables).
          *
          * @param vectorName     The base name of the vector (variable name given on the input file line before
          *                       the vector starts).
          * @param variable       Pointer to an array of doubles in which to push values obtained from the input file.
          */
        template<class T>
        void getVector (string vectorName, vector<T>* variable) {
          //assert(fileSet_);
          //assert(fileParsed_);
          if ( !fileSet_ || !fileParsed_ ) {
            ostringstream errmsg;
            errmsg  << "WARNING: inputFile::getVector(): the input file has not been set " 
                    << "or has not been parsed. " << endl;
            throw runtime_error( errmsg.str() );
          }
 
          bool atEnd = false;
          int z = 0;
          int lineNumber;
          vector<int> lines;

          while (!atEnd) {
            // reconstruct vector name
            stringstream out;
            out << vectorName << "_vector" << (z+1);
            string vecname;
            vecname = out.str();
        
            map<string, int>::iterator iLine = input_map_.find(vecname);
        
            // if not even the first element can be found, print warning
            if ( z == 0 && iLine == input_map_.end() ) {
              ostringstream errmsg;
              errmsg  << "WARNING: inputFile::getVector(): the requested vector " << vecname
                      << " does not exist in the input file." << endl;
              throw runtime_error( errmsg.str() );
            } else if ( iLine == input_map_.end() ) {
              atEnd = true;
              break;
            }
        
            // now grab the line number
            lineNumber = iLine->second;
            lines.push_back(iLine->second);
            ++z;
          } //end while
        
          (*variable).clear();
          for ( vector<int>::iterator ilines = lines.begin();
                ilines != lines.end(); ++ilines) {
            
            // skip forward to line of interest
            lineNumber = *ilines;
            ifstream ifile( filename_.c_str() );
            string s1;
            for (int i = 1; i < lineNumber; ++i) {
              getline(ifile, s1);
            }
            
            // extract and assign the variable value
            double tempvariable;
            ifile >> tempvariable;    // the vector element's value
                                      // (we don't care about the label)
            (*variable).push_back(tempvariable);

            ifile.close();
      
          }

        }

        /** @brief      This grabs a vector from the input file (more precisely, from the list/map of input file variables).
          *
          * @param vectorName     The base name of the vector (variable name given on the input file line before
          *                       the vector starts).
          * @param variable       Pointer to an array of doubles in which to push values obtained from the input file.
          * @param varNames       Pointer to an array of strings in which to push names of vector elements
          */
        template<class T>
        void getVector (string vectorName,
                        vector<T>* variable,
                        vector<string>* varNames) {
          //assert(fileSet_);
          //assert(fileParsed_);
          if ( !fileSet_ || !fileParsed_ ) {
            ostringstream errmsg;
            errmsg  << "WARNING: inputFile::getVector(): the input file has not been set " 
                    << "or has not been parsed. " << endl;
            throw runtime_error( errmsg.str() );
          }
 
          bool atEnd = false;
          int z = 0;
          int lineNumber;
          vector<int> lines;


          // grab all the line numbers that are needed
          while (!atEnd) {
            // reconstruct vector name
            stringstream out;
            out << vectorName << "_vector" << (z+1);
            string vecname;
            vecname = out.str();
        
            map<string, int>::iterator iLine = input_map_.find(vecname);
        
            // if not even the first element can be found, print warning
            if ( z == 0 && iLine == input_map_.end() ) {
              ostringstream errmsg;
              errmsg  << "WARNING: inputFile::getVector(): the requested vector " << vecname
                      << " does not exist in the input file." << endl;
              throw runtime_error( errmsg.str() );
            } else if ( iLine == input_map_.end() ) {
              atEnd = true;
              break;
            }
        
            // now grab the line number
            lineNumber = iLine->second;
            lines.push_back(iLine->second);
            ++z;
          } //end while


          // actually seek out the line numbers to get vector elements' 
          // values and labels/names
          (*variable).clear();
          (*varNames).clear();
          for ( vector<int>::iterator ilines = lines.begin();
                ilines != lines.end(); ++ilines) {
            
            // skip forward to line of interest
            lineNumber = *ilines;
            ifstream ifile( filename_.c_str() );
            string s1;
            for (int i = 1; i < lineNumber; ++i) {
              getline(ifile, s1);
            }
            
            // extract and assign the variable value
            double tempvariable;
            ifile >> tempvariable;      // the vector element's value
            ifile >> s1;                // the vector element's label
            
            // get keyword from between ()'s
            const unsigned int pos1 = s1.find("(",0);
            const unsigned int pos2 = s1.find(")",0);
            if (pos1 == s1.size() || pos2 < pos1 ) {
              ostringstream errmsg;
              errmsg << "ERROR: inputFile::readStreams(): Could not find matching delimiters '(' and ')' in line "
                     << lineNumber << "of input file. Indicate comment lines using '-' or '#'." << endl;
              throw runtime_error( errmsg.str() );
            }
            string elementName = s1.substr(pos1+1, pos2-pos1-1);

            // save it all
            (*variable).push_back(tempvariable);
            (*varNames).push_back(elementName);

            ifile.close();

          } //end for
      
        } //end method


 
};

#endif

