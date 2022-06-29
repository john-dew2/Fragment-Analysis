#ifndef _ANALYSIS_GUARD
#define _ANALYSIS_GUARD 1

#include <map>
#include <vector>

#include <openbabel/mol.h>


class FragmentAnalysis{
	
	public:
		//general functions
		FragmentAnalysis(std::vector<OpenBabel::OBMol*>& linkers, std::vector<OpenBabel::OBMol*>& bricks, OpenBabel::OBMol* subject);
		void doAnalysis();
		
		//utility functions
		double tanimotoCalc(OpenBabel::OBMol* mol1, OpenBabel::OBMol* mol2);
		
		//output functions
		void writeFreqReport(std::string fileName, map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>>);
		void writeDistReport(std::string fileName, map<double, int> map)
		void printMap(map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> map);
		void printDistribution(map<double, int> map);
		
		//analysis functions
		map<OpenBabel::OBMol*, std::vector<OpenBabel::OBMol*>> fragmentAnalysis(std::vector<OpenBabel::OBMol*>& fragments);
		map<double, int> distributionAnalysis(OpenBabel::OBMol* subject, std::vector<OpenBabel::OBMol*>& fragments);
		
	private:
		std::vector<OpenBabel::OBMol*>& _linkers;
		std::vector<OpenBabel::OBMol*>& _bricks;
		OpenBabel::OBMol* _subject;
};


#endif