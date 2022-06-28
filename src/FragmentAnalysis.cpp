

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
#include <cmath>
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
	
	double FragmentAnalysis::tanimotoCalc(OpenBabel::OBMol* mol1, OpenBabel::OBMol* mol2)
	{
		//create two vectors
		std::vector<unsigned int> vector1;
		std::vector<unsigned int> vector2;

		//create a fingerprint
		OpenBabel::OBFingerprint* fpType1 = OpenBabel::OBFingerprint::FindFingerprint("");
		fpType1->GetFingerprint(mol1, vector1);
		fpType1->GetFingerprint(mol2, vector2);

		//calculate the tc of the two vectors
		double tanimoto = OpenBabel::OBFingerprint::Tanimoto(vector1, vector2);
		
		std::cout << tanimoto << std::endl;
		
		return tanimoto;
		
	}

	void FragmentAnalysis::writeReport(std::string fileName, map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> map)
	{
		const char* name;
		int numSimilar;
		std::string similarNames;
		std::vector<OpenBabel::OBMol*> molVector;
		
		
		// Create and open a text file
		ofstream reportFile(fileName);
		
		for (const auto& pair : map) 
		{
			similarNames = "";
			name = pair.first->GetTitle();
			molVector = pair.second;
			numSimilar = molVector.size();
			
			for (unsigned i = 0; i < numSimilar; i++)
			{
				similarNames += molVector[i]->GetTitle();
				if (!(i+1 >= numSimilar)) similarNames += ", ";
			}
			reportFile << name << " has " << numSimilar << " similar elements. Their names are: " << similarNames << std::endl;
		}
		// Close the file
		reportFile.close();
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
				tc = tanimotoCalc(frag1, frag2);
				if (tc > TC_THRESHOLD){
					similarities.push_back(frag2);
				}
			}
			//add a new map entry before the next iteration
			frequencyMap.insert(pair<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>>(frag1, similarities));
		}
		
		return frequencyMap;
		
	}
	
	map<double, int> FragmentAnalysis::distributionAnalysis(OpenBabel::OBMol* subject, std::vector<OpenBabel::OBMol*>& fragments){
		std::vector<double> coefficients{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
		map<double, int> distribution;
		double tc;
		OpenBabel::OBMol* fragment;
		
		//populate the map
		for (unsigned x = 0; x < coefficents.size(); x++)
		{
			distribution.insert(pair<double, int>(coefficents[x], 0));
		}
		
		//compare each brick to the subject and note the tc value between the two molecules
		for (unsigned y = 0; y < fragments.size(); y++)
		{
			fragment = fragments[y];
			tc = tanimotoCalc(subject, fragment);
			
			tc *= 10;
			int temp = floor(tc);
			tc = temp * 0.10;
			
			temp = (distribution.find(tc)) + 1
			distribution.insert(pair<double, int>(tc, temp));
		}
		
		return distribution;
		
	}
	
	void FragmentAnalysis::printMap(map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> map)
	{
		//for every pair in the map, print its name and the nummber of elements similar
		for (const auto& pair : map) 
		{
			const char* name = pair.first->GetTitle();
			int length = pair.second.size();
			std::cout << name << ": " << length << " fragments similar, ";
		}
    }


	
	void FragmentAnalysis::doFragmentAnalysis(){
		
		//create two maps
		//map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> brickMap;
		//map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>>linkerMap;
		
		//run a frequency analysis on the two sets
		//brickMap  = freqAnalysis(FragmentAnalysis::_bricks);
		//linkerMap = freqAnalysis(FragmentAnalysis::_linkers);
		
		//print the contents
		//std::cout<<"Brick Map contents: ";
		//printMap(brickMap);
		//std::cout << std::endl;
		
		//std::cout<<"Linker Map contents: ";
		//printMap(linkerMap);
		//std::cout << std::endl;
		
		distributionBrick  = distributionAnalysis(_bricks[0], _bricks);
		distributionLinker = distributionAnalysis(_linkers[0], _linkers);
		
		printMap(distributionBrick);
		printMap(distributionLinker);
		
		//writeReport("brick-report.txt", brickMap);
		//writeReport("linker-report.txt", linkerMap);
		
	}