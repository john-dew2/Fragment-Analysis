

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


	FragmentAnalysis::FragmentAnalysis(std::vector<OpenBabel::OBMol*>& linkers, std::vector<OpenBabel::OBMol*>& bricks) : _linkers(linkers), _bricks(bricks)
	{
		//FragmentAnalysis::_linkers = linkers;
		//FragmentAnalysis::_bricks = bricks;
	}
	
	
	//
	// calculate the tanimoto coefficient of two moelcules
	//
	double FragmentAnalysis::tanimotoCalc(OpenBabel::OBMol* mol1, OpenBabel::OBMol* mol2)
	{
		//create two vectors to assign a fingerprint to
		std::vector<unsigned int> vector1;
		std::vector<unsigned int> vector2;

		//create a fingerprint object and calculate the fingerprint for the two molecules
		OpenBabel::OBFingerprint* fpType1 = OpenBabel::OBFingerprint::FindFingerprint("");
		fpType1->GetFingerprint(mol1, vector1);
		fpType1->GetFingerprint(mol2, vector2);

		//calculate the tc of the two vectors
		double tanimoto = OpenBabel::OBFingerprint::Tanimoto(vector1, vector2);
		
		return tanimoto;
	}
	
	
	//
	// write a report to a file with the name of the molecule and what other molecules are similar to that one
	//
	void FragmentAnalysis::writeReport(std::string fileName, map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> map)
	{
		const char* name;
		int numSimilar;
		std::string similarNames;
		std::vector<OpenBabel::OBMol*> molVector;
		
		// Create and open a text file
		ofstream reportFile(fileName);
		
		//for each entry in the map
		for (const auto& pair : map) 
		{
			//grab the name and the size of the vector
			name = pair.first->GetTitle();
			molVector = pair.second;
			numSimilar = molVector.size();
			
			//for each element in the list
			for (unsigned i = 0; i < numSimilar; i++)
			{
				//grab the names of the list
				similarNames += molVector[i]->GetTitle();
				if (!(i+1 >= numSimilar)) similarNames += ", ";
			}
			//concatenate the report and clear the similarNames variable
			reportFile << name << " has " << numSimilar << " similar elements. Their names are: " << similarNames << std::endl;
			similarNames = "";
		}
		// Close the file
		reportFile.close();
	}
	
	
	//
	// pritns the ouput of a map in the format of pair<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>>
	//
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
	
	void FragmentAnalysis::printDistribution(map<double, int> map)
	{
		//for every pair in the map, print its name and the nummber of elements similar
		for (const auto& pair : map) 
		{
			double category = pair.first;
			int amount = pair.second;
			std::cout << "Tanimoto Coefficient: " << category << ", contains " << amount << " fragments." << std::endl;
		}
	}


	//
	// Runs an analysis on all fragments and identifies which fragments are similar of the same
	//
	map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> FragmentAnalysis::fragmentAnalysis(std::vector<OpenBabel::OBMol*>& fragments)
	{
		std::vector<OpenBabel::OBMol*> similarities;
		double tanimoto;
		double TC_THRESHOLD = 0.8;
		map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> frequencyMap;
		
		//for each brick, compare it to every brick
		for (unsigned i = 0; i < fragments.size(); i++)
		{
			//grab one brick
			OpenBabel::OBMol* frag1 = fragments[i];
			for (unsigned j = 0; j < fragments.size(); j++)
			{
				//grab a second brick to compare
				OpenBabel::OBMol* frag2 = fragments[j];
				
				//find the tc and add the brick to the list if the tc value is > the threshold
				tanimoto = tanimotoCalc(frag1, frag2);
				if (tanimoto > TC_THRESHOLD) similarities.push_back(frag2);
			}
			//add a new map entry before the next iteration, and clear the similarity list after each iteration
			frequencyMap.insert(pair<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>>(frag1, similarities));
			similarities.clear();
		}
		
		return frequencyMap;
	}
	
	
	//
	// Has 11 categories representing the tenths place between 0 and 1 inclusive (0.0, 0.1, ..., 1.0) and collects data on how frequently 
	// a fragment will have a certain tanimoto coefficent when being comapred to a single molecule
	//	
	map<double, int> FragmentAnalysis::distributionAnalysis(OpenBabel::OBMol* subject, std::vector<OpenBabel::OBMol*>& fragments)
	{
		std::vector<double> coefficients{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
		map<double, int> distribution;
		double tanimoto;
		OpenBabel::OBMol* fragment;
		int temp;
		
		//populate the map with default values so we can grab the entrys later
		for (unsigned x = 0; x < coefficients.size(); x++)
		{
			distribution.insert(pair<double, int>(coefficients[x], 0));
		}
		
		//for all the fragments
		for (unsigned y = 0; y < fragments.size(); y++)
		{
			//grab a fragment and compare them to the subject
			fragment = fragments[y];
			tanimoto = tanimotoCalc(subject, fragment);
			
			//get rid of the hundreths place and beyond because we are only intrested in the tenths place
			tanimoto = floor(tanimoto * 10) * 0.10;
			
			//add one to the map entry and reinsert it
			temp = 0;
			temp = distribution[tanimoto] + 1;
			distribution[tanimoto] = temp;
		}
		return distribution;
	}


	//
	// Main
	//
	void FragmentAnalysis::doAnalysis(){
		
		//create two maps
		map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> brickMap;
		map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>>linkerMap;
		
		//run a frequency analysis on the two sets
		brickMap  = fragmentAnalysis(FragmentAnalysis::_bricks);
		linkerMap = fragmentAnalysis(FragmentAnalysis::_linkers);
		
		//print the contents
		std::cout<<"Brick Map contents: ";
		printMap(brickMap);
		std::cout << std::endl;
		
		std::cout<<"Linker Map contents: ";
		printMap(linkerMap);
		std::cout << std::endl;
		
		map<double, int> distributionBrick  = distributionAnalysis(_bricks[0], _bricks);
		map<double, int> distributionLinker = distributionAnalysis(_linkers[0], _linkers);
		
		printDistribution(distributionBrick);
		printDistribution(distributionLinker);
		
		//writeReport("brick-report.txt", brickMap);
		//writeReport("linker-report.txt", linkerMap);
		
	}