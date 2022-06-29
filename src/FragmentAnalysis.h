#ifndef _ANALYSIS_GUARD
#define _ANALYSIS_GUARD 1

#include <map>
#include <vector>

#include <openbabel/mol.h>


class FragmentAnalysis{
	public:
		FragmentAnalysis(std::vector<OpenBabel::OBMol*>& linkers, std::vector<OpenBabel::OBMol*>& bricks, OpenBabel::OBMol* subject);
		double tanimotoCalc(OpenBabel::OBMol* mol1, OpenBabel::OBMol* mol2);
		void writeReport(std::string fileName, map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>>);
		void printMap(map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> map);
		void printDistribution(map<double, int> map);
		map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> fragmentAnalysis(std::vector<OpenBabel::OBMol*>& fragments);
		void doAnalysis();
		map<double, int> distributionAnalysis(OpenBabel::OBMol* subject, std::vector<OpenBabel::OBMol*>& fragments);
		
	private:
		std::vector<OpenBabel::OBMol*>& _linkers;
		std::vector<OpenBabel::OBMol*>& _bricks;
		OpenBabel::OBMol* _subject;
};

// void Analysis::doFrequencyAnalysis(std::vector<OpenBabel::OBMol*>& linkerList, std::vector<OpenBabel::OBMol*>& brickList){
		
		// map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> brickMap;
		// map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>>linkerMap;
		
		// brickMap  = freqAnalysis(brickList);
		// linkerMap = freqAnalysis(linkerList);
		
	// }

#endif