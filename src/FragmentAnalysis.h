#ifndef _ANALYSIS_GUARD
#define _ANALYSIS_GUARD 1

#include <map>
#include <vector>

#include <openbabel/mol.h>


class FragmentAnalysis{
	public:
		FragmentAnalysis(std::vector<OpenBabel::OBMol*>& linkers, std::vector<OpenBabel::OBMol*>& bricks);
		double tanimoto_calc(OpenBabel::OBMol* mol1, OpenBabel::OBMol* mol2);
		void write_out();
		map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> freqAnalysis(std::vector<OpenBabel::OBMol*>& fragments);
		void doFragmentAnalysis();
		
	private:
		std::vector<OpenBabel::OBMol*>& linkers;
		std::vector<OpenBabel::OBMol*>& bricks;
};

// void Analysis::doFrequencyAnalysis(std::vector<OpenBabel::OBMol*>& linkerList, std::vector<OpenBabel::OBMol*>& brickList){
		
		// map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> brickMap;
		// map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>>linkerMap;
		
		// brickMap  = freqAnalysis(brickList);
		// linkerMap = freqAnalysis(linkerList);
		
	// }

#endif