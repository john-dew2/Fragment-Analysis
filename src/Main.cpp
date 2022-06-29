/*
 *  This file is part of esynth.
 *
 *  esynth is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  esynth is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with esynth.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdio>
#include <sstream>
#include <cstdlib>
#include <cctype>
#include <mcheck.h>

//
// Open Babel
//
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/generic.h>
#include <openbabel/atom.h>
#include <openbabel/bond.h>
#include <openbabel/groupcontrib.h>


//
// This project molecular representation
//
#include "Atom.h"
#include "Bond.h"
#include "Molecule.h"
#include "Brick.h"
#include "Linker.h"

//
// File processing in / out.
//
#include "OBWriter.h"
#include "Options.h"

//
// Synthesis-Based Functionality
//
#include "EdgeAnnotation.h"

#include "Utilities.h"
#include "IdFactory.h"
#include "Constants.h"
#include "FragmentAnalysis.h"


//
// Global set of linkers and bricks read from the input files, and another global set to store the molecules
// Global string and molecule to represent the file name and actual molecular object that will be the subject of distribution analysis
//
std::vector<Linker*> linkers;
std::vector<Brick*> bricks;
std::vector<OpenBabel::OBMol*> brickmol;
std::vector<OpenBabel::OBMol*> linkmol;
std::string subject_file;
OpenBabel::OBMol* subject = new OpenBabel::OBMol();

void Cleanup(std::vector<Linker*>& linkers, std::vector<Brick*>& bricks);

//
// Split the molecule into atoms/bonds as well as connectivity information
//
bool splitMolecule(std::ifstream& infile, std::string& name,
                   std::string& prefix, std::string& suffix)
{
    prefix = "";
    suffix = "";

    std::string line = "";

    // Eat #### in large files (if it exists)
    eatWhiteLines(infile); 
    if (infile.peek() == '#')
    {
        getline(infile, line);
        name = line;
        eatWhiteLines(infile);
    }

    getline(infile, line);
    prefix += line + '\n';

    // Nothing left to read...
    if (infile.eof() || infile.fail()) return false;

    // Read the prefix (end indicated by END)
    while(line.find("END") == std::string::npos)
    {
        getline(infile, line);
        prefix += line + '\n';
    }

    // Add '$$$$' to the prefix.
    // prefix += "\n$$$$";

    // Set suffix equal to remainder of the file
    while (line.find("$$$$") == std::string::npos)
    {
        if (!getline(infile, line)) break;
        suffix += line + '\n';
    }

    return true;
}

//
// Given a molecule (mol )and its extraneous information (mtype, name, suffix), convert it to a
// local linker or brick, and feed the information as data for the actual molecular object (mol)
//
/*
Molecule* createLocalMolecule(OpenBabel::OBMol* mol, MoleculeT mType,
                              const std::string& name, std::string& suffix)
{
    //
    // Add the suffix as comment data to the actual OBMol object.  
    //
    OpenBabel::OBCommentData* cData = new OpenBabel::OBCommentData();
    cData->SetAttribute("Comment");
    cData->SetData(suffix);
    mol->SetData(cData);

    //
    // Create this particular molecule type based on the name of the file.
    //
    if (mType == LINKER)
    {
        return new Linker(mol, name);
    }
    else if (mType == BRICK)
    {
        return new Brick(mol, name);
    }
    
    return 0;
}
*/

//
// Takes a molecule and its type and appends it to either the linker, bricks, or the catch all category
//
/*
void addLocalMolecule(char type, Molecule* molecule)
{
    // linker
    if (type == 'l')
    {
        linkers.push_back(static_cast<Linker*>(molecule));        
    }
    // brick or brick
    else if (type == 'r' || type == 'b')
    {
        bricks.push_back(static_cast<Brick*>(molecule));        
    }
    else
    {
        std::cerr << "Unrecognized file prefix: "
                  << type << "Expected: l->linker, b->brick" << std::endl;
    }
}
*/

//
// Takes a molecule and its type and appends it to either the linker, bricks, or the catch all category
//
void addOBMolecule(char type, OpenBabel::OBMol* molecule)
{
	//if we are intrested in distribution analysis, then the first file will be our subject
	if (Options::IS_SUBJECT_FILE)
	{
		Options::IS_SUBJECT_FILE = false;
		subject = molecule;			
	}
	if (type == 'l')
	{
		linkmol.push_back(molecule);
	}
	// brick or brick
	else if (type == 'r' || type == 'b')
	{
		brickmol.push_back(molecule);        
	}
}

void readMoleculeFile(const char* fileName)
{
	//
    // Input parser conversion functionality for Open babel
    //
    OpenBabel::OBConversion obConversion;
    obConversion.SetInFormat("SDF");

    //
    // Open the file, split the current molecule into Molecule Data (prefix)
    // and Our Data (Suffix)
    //
    std::ifstream infile;
    infile.open(fileName);

    std::string name = "UNKNOWN";
    std::string prefix = "";
    std::string suffix = "";
    
    while(splitMolecule(infile, name, prefix, suffix))
    {
        //
        // If the name of molecule is not given, overwrite it with
        // the name of the file.
        //
        if (name == "UNKNOWN")
        {
           name = "####   ";
           name += fileName;
           name += "    ####";
        }

        if (g_debug_output) std::cerr << "Name: " << std::endl << name << std::endl;
        if (g_debug_output) std::cerr << "Prefix: " << std::endl << prefix << std::endl;
        if (g_debug_output) std::cerr << "Suffix: " << std::endl << suffix << std::endl;

        // Create and parse using Open Babel
        OpenBabel::OBMol* mol = new OpenBabel::OBMol();
        //bool notAtEnd = 
		obConversion.ReadString(mol, prefix);

        // Assign all needed data to the molecule (comment data)
        //Molecule* local = createLocalMolecule(mol,
         //                                     tolower(fileName[0]) == 'l' ? LINKER : BRICK,
         //                                     name, suffix);


        // calculate the molecular weight, H donors and acceptors and the plogp
        //local->openBabelPredictLipinski();

        // add to logfile
        if (Molecule::isOpenBabelLipinskiCompliant(*mol))
        {
            std::ofstream logfile("synth_log_initial_fragments_logfile.txt",
                                  std::ofstream::out | std::ofstream::app); // append
            logfile << fileName << "\nMolWt = " << local->getMolWt() << "\n";
            logfile << "HBD = " << local->getHBD() << "\n";
            logfile << "HBA1 = " << local->getHBA1() << "\n";
            logfile << "logP = " << local->getlogP() << "\n";
            logfile << std::endl;
            logfile.close();
        }
        else std::cerr << "Main: predictLipinski failed somehow!" << endl;

        if (g_debug_output) std::cout << "Local: " << *local << "|" << std::endl;
    
        // Add to the linker or brick list as needed.
        //addLocalMolecule(tolower(fileName[0]), local); 
		
		addOBMolecule(tolower(fileName[0]), mol); 

			
		//COME BACK AND FIX
		//delete mol;
	}
}


//
// Parse each input data files
//
bool readInputFiles(const Options& options)
{
	if (Options::DISTRIBUTION_ANALYSIS){
		subject_file = options.inFiles.front();
		options.inFiles.erase(options.inFiles.begin());
		
	}
	
    for (std::vector<std::string>::const_iterator it = options.inFiles.begin();
         it != options.inFiles.end(); it++)
    {
        char charPrefix = tolower((*it)[0]); 
        if (charPrefix != 'l' && charPrefix != 'r' && charPrefix != 'b')
        {
            cerr << "Unexpected file prefix: \'" << (*it)[0]
                 << "\' with file " << *it << endl;
            return false;
        }

        readMoleculeFile((*it).c_str());
    }

    return true;
}


int main(int argc, char** argv)
{	
	if (argc < 2)
	{
		std::cerr << "Usage: <program> [SDF-file-list] -o <output-file> -v <validation-file>"
				  << " -pool <#obgen-threads>" << std::endl;
		return 1;
	}

    //
    // Global options object.
    //
	Options options(argc, argv);
	if (!options.parseCommandLine())
	{
		std::cerr << "Command-line parsing failed; exiting." << std::endl;
		return 1;
	}

    // 
    // Output command-line option information
    //
	if (Options::THREADED) std::cout << "Threaded execution :P." << std::endl;
	else if (Options::SERIAL) std::cout << "Serial execution." << std::endl;
	else
	{
		std::cerr << "Neither serial nor threaded specified; exiting." << std::endl;
		return 1;
	}

    // std::cout << "SMI Comparison Level: " << Options::SMI_LEVEL_BOUND << std::endl;
	std::cout << "Probability Filtration Level: "
			  << Options::PROBABILITY_PRUNE_LEVEL_START << std::endl;

    // Printing the specified Tanimoto value to the user.
    // std::cerr << "Tanimoto Coefficient Threshold Specified: "
    //           << Options::TANIMOTO << std::endl;
	if (!Options::SMI_ONLY)
	{
		std::cerr << "OBGEN output thread pool size: "
				  << Options::OBGEN_THREAD_POOL_SIZE << std::endl;
	}
	if (!readInputFiles(options)) 
	{
		std::cout << "Did not read input files" <<std:: endl;
		return 1;
	}    
	
	//Added 6/23 - Dewey
	FragmentAnalysis analyzer(linkmol, brickmol, subject);
	analyzer.doAnalysis();
      
	Cleanup(linkers, bricks);
    
	return 0;
}

void Cleanup(std::vector<Linker*>& linkers, std::vector<Brick*>& bricks)
{
    for (unsigned ell = 0; ell < linkers.size(); ell++)
    {
        delete linkers[ell];
    }

    for (unsigned r = 0; r < bricks.size(); r++)
    {
        delete bricks[r];
    }
}
