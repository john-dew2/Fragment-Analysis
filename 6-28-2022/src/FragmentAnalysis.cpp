

//Frequency Analysis

#include <vector>
#include<iterator> // for iterators
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdio>
#include <sstream>
#include <cstdlib>
#include <cctype>
//#include <mcheck.h>


//
// Open Babel
//
//#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
//#include <openbabel/generic.h>
//#include <openbabel/atom.h>
//#include <openbabel/bond.h>
//#include <openbabel/groupcontrib.h>
#include <openbabel/fingerprint.h>


//
// This project molecular representation
//
//#include "Atom.h"
//#include "Bond.h"
//#include "Molecule.h"
//#include "Brick.h"


//
// File processing in / out.
//
#include "OBWriter.h"
#include "Options.h"
//#include "Validator.h"



//
// Synthesis-Based Functionality
//

//#include "EdgeAnnotation.h"
//#include "Instantiator.h"

//#include "Utilities.h"
//#include "IdFactory.h"
#include "Constants.h"



#include "FragmentAnalysis.h"
	FragmentAnalysis::FragmentAnalysis(std::vector<OpenBabel::OBMol*>& linkers, std::vector<OpenBabel::OBMol*>& bricks) : _linkers(linkers), _bricks(bricks){
		//FragmentAnalysis::_linkers = linkers;
		//FragmentAnalysis::_bricks = bricks;
	}
	
	double FragmentAnalysis::tanimoto_calc(OpenBabel::OBMol* mol1, OpenBabel::OBMol* mol2)
	{
		//create two vectors
		std::vector<unsigned int> vector1;
		std::vector<unsigned int> vector2;
		
		//create an object to convert our fragments
		//OBConversion conv;
		//conv.SetInFormat("smi");
		
		//convert our fragments to a molecular object
		//OpenBabel::OBMol* frag1 = new OpenBabel::OBMol();
		//OpenBabel::OBMol* frag2 = new OpenBabel::OBMol();
		
		//conv.ReadString(frag1, brick1);
		//conv.ReadString(frag2, brick2);

		//create a fingerprint
		OpenBabel::OBFingerprint* fpType1 = OpenBabel::OBFingerprint::FindFingerprint("");
		
		fpType1->GetFingerprint(mol1, vector1);
		fpType1->GetFingerprint(mol2, vector2);

		// calculate tc
		double tanimoto = OpenBabel::OBFingerprint::Tanimoto(vector1, vector2);
		
		return tanimoto;
		
	}

	void FragmentAnalysis::write_out()
	{
		
		string printStr;

		// Create and open a text file
		ofstream brickFile("BrickAnalysis.txt");

		// Write to the file
		//for (pair<char*, char**> entry : brickMap ) {
			//printStr += entry.first + entry.second;
		//}
		brickFile << printStr;
		printStr = "";

		// Close the file
		brickFile.close();
		
		// Create and open a text file
		ofstream linkerFile("LinkerAnalysis.txt");

		// Write to the file
			//for (pair<char*, char**> entry : linkerMap ) {
			//printStr += entry.first + entry.second;
		//}
		linkerFile << printStr;

		// Close the file
		linkerFile.close();

	}

	map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> FragmentAnalysis::freqAnalysis(std::vector<OpenBabel::OBMol*>& fragments)
	{
		
		std::vector<OpenBabel::OBMol*> similarities;
		double tc;
		double TC_THRESHOLD = 0.8;
		map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> frequencyMap;
		
		//for each brick, compare it to every brick
		for (unsigned i = 0; i < fragments.size(); i++)
		{
			//clear the similarity list for each iteration
			similarities.clear();
			
			//grab one brick
			OpenBabel::OBMol* frag1 = fragments[i];
			for (unsigned j = 0; j < fragments.size(); j++)
		{
				//grab a second brick to compare
				OpenBabel::OBMol* frag2 = fragments[j];
				
				//find the tc and add the brick to the list if the tc value is > the threshold
				tc = tanimoto_calc(frag1, frag2);
				if (tc > TC_THRESHOLD){
					similarities.push_back(frag2);
				}
		}
			//add a new map entry before the next iteration
			frequencyMap.insert(pair<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>>(frag1, similarities));
		}
		
		return frequencyMap;
		
	}
	
	void FragmentAnalysis::printMap(map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> map){
		/* for (const auto& [key, value] : map) {
			std::cout << '[' << key << "] = " << value << "; ";
		} */
		
		for (const auto& n : map) {
			const char* name = n.first->GetTitle();
			int length = n.second.size();
			std::cout << name << " = " << length << "; ";
		}
    }


	
	void FragmentAnalysis::doFragmentAnalysis(){
		
		map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> brickMap;
		map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>>linkerMap;
		
		brickMap  = freqAnalysis(FragmentAnalysis::_linkers);
		linkerMap = freqAnalysis(FragmentAnalysis::_bricks);
		
		//std::cout << "Here is the brickMap: " << std::endl << brickMap
		
		std::cout<<"Brick Map contents: " << std::endl;
		toString(brickMap);
		
		std::cout<<"Linker Map contents: " << std::endl;
		toString(linkerMap);
		//std::cout << "Here is the linkerMap: " << std::endl << linkerMap
		
	}