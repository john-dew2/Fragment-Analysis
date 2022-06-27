#ifndef _ANALYSIS_GUARD
#define _ANALYSIS_GUARD 1

#include <map>
#include <vector>

#include <openbabel/mol.h>


class Analysis{
	public:
		static double tanimoto_calc(OpenBabel::OBMol* mol1, OpenBabel::OBMol* mol2);
		static void write_out();
		static map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> freqAnalysis(std::vector<OpenBabel::OBMol*>& fragments);
		static void doFrequencyAnalysis(std::vector<OpenBabel::OBMol*>& linkerList, std::vector<OpenBabel::OBMol*>& brickList);
};

// void Analysis::doFrequencyAnalysis(std::vector<OpenBabel::OBMol*>& linkerList, std::vector<OpenBabel::OBMol*>& brickList){
		
		// map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> brickMap;
		// map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>>linkerMap;
		
		// brickMap  = freqAnalysis(brickList);
		// linkerMap = freqAnalysis(linkerList);
		
	// }

#endif