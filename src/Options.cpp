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

#include <string>
#include <vector>
#include <iostream>
#include <cstring>
#include <cstdlib>

#include <sys/stat.h>

#include "Options.h"
#include "Constants.h"


double Options::TANIMOTO = 0.95;
std::string Options::shmPath = "/run/shm";
std::string Options::writerPath = "./";
bool Options::THREADED = false;
bool Options::SERIAL = true;
bool Options::OPENBABEL = true;
bool Options::SMI_ONLY = false;
bool Options::USE_LIPINSKI = false;
bool FREQUENCY_ANALYSIS = false;
bool DISTRIBUTION_ANALYSIS = false;
bool IS_SUBJECT_FILE = false;
bool VALIDATE = false;
unsigned Options::OBGEN_THREAD_POOL_SIZE = 15;
//unsigned Options::SMI_LEVEL_BOUND = 3;
unsigned Options::PROBABILITY_PRUNE_LEVEL_START = 5;
std::string Options::OUTPUT_DIR_SUFFIX = "";

double Options::SA_THRESHOLD = 5.0;
std::string Options::PYTHON_MODULE_NAME;
std::string Options::PYTHON_FUNCTION_NAME;

Options::Options(int argCount, char** vals) : argc(argCount), argv(vals)
{
    // Default values
    outFile = "molecules.sdf";
    outFileSMI = "molecules.smi";
    validationFile = "";

    Options::TANIMOTO = 0.95;
    Options::THREADED = false;
}


bool Options::acquireEnvironmentVariable(const std::string& variable,
                                         const std::string& suffix, std::string& value)
{
    //
    // Acquire the environment variablei value
    //
    char* path = getenv(variable.c_str());
    if (!path)
    {
        std::cerr << variable << " environment variable not specified."
                  << std::endl; 
        value = "";
        return false;
    }
    else
    {
        std::cerr << variable << " environmental path: " << path << std::endl;
    }

    // Create the path of the actual file.
    std::string filePath = path;

    if (filePath[filePath.size() - 1] != '/') filePath += '/'; 

    // Assign to the return variable.
    value = filePath;

    filePath += suffix;

    // Check that such a file exists.
    struct stat buffer;   
    if (stat (filePath.c_str(), &buffer) != 0)
    {
        std::cerr << "Specified path: " << filePath << " does not exist." << std::endl;
        value = "";
        return false;
    }

    return true;
}

bool Options::AnalyzeEnvironment()
{
    bool returnVal = acquireEnvironmentVariable("COMPLIANT_WRITER", COMPLIANT_EXE, Options::writerPath);
    returnVal = acquireEnvironmentVariable("SHM_PATH", "", Options::shmPath) && returnVal;

    return returnVal;
}

bool Options::parseCommandLine()
{
    for (int i = 1; i < argc; i++)
    {
        if (argv[i][0] == '-')
        {
            if (!handleOption(i)) return false;
        }
        else
        {
            inFiles.push_back(argv[i]);
        }
    }

    return true;
}

//
// Deal with the actual options specified on the command-line.
//
bool Options::handleOption(int& index)
{
    if (strcmp(argv[index], "-o") == 0)
    {
        outFile = argv[++index];
        return true;
    }
    if (strcmp(argv[index], "-v") == 0)
    {
        VALIDATE = true;
        validationFile = argv[++index];
        return true;
    }
    if (strncmp(argv[index], "-tc", 3) == 0)
    {
        // not directly following; e.g. -tc 0.95
        if (strcmp(argv[index], "-tc") == 0)
        {
            Options::TANIMOTO = atof(argv[++index]);
        }
        // -tc followed by a double value; e.g. -tc0.95
        else
        {
            Options::TANIMOTO = atof(&argv[index][3]);
        }
        return true;
    }

    if (strncmp(argv[index], "-mw", 3) == 0)
    {
        if (strcmp(argv[index], "-mw") == 0)
            MOLWT_UPPERBOUND = atof(argv[++index]);
        else
            MOLWT_UPPERBOUND = atof(&argv[index][3]);
        return true;
    }
    if (strncmp(argv[index], "-sa", 3) == 0)
    {
        if (strcmp(argv[index], "-sa") == 0)
            SA_THRESHOLD = atof(argv[++index]);
        else
            SA_THRESHOLD = atof(&argv[index][3]);
        return true;
    }
    if (strncmp(argv[index], "-hd", 3) == 0)
    {
        if (strcmp(argv[index], "-hd") == 0)
            HBD_UPPERBOUND = atof(argv[++index]);
        else
            HBD_UPPERBOUND = atof(&argv[index][3]);
        return true;
    }
    if (strncmp(argv[index], "-ha", 3) == 0)
    {
        if (strcmp(argv[index], "-ha") == 0)
            HBA1_UPPERBOUND = atof(argv[++index]);
        else
            HBA1_UPPERBOUND = atof(&argv[index][3]);
        return true;
    }
    if (strncmp(argv[index], "-lp", 3) == 0)
    {
        if (strcmp(argv[index], "-lp") == 0)
            LOGP_UPPERBOUND = atof(argv[++index]);
        else
            LOGP_UPPERBOUND = atof(&argv[index][3]);
        return true;
    }
    if (strncmp(argv[index], "-hl", 3) == 0)
    {
        if (strcmp(argv[index], "-hl") == 0)
            HIERARCHICAL_LEVEL_BOUND = atoi(argv[++index]);
        else
            HIERARCHICAL_LEVEL_BOUND = atoi(&argv[index][3]);
        return true;
    }
/*
    if (strncmp(argv[index], "-smi-level", 10) == 0)
    {
        if (strcmp(argv[index], "-smi-level") == 0)
            SMI_LEVEL_BOUND = atoi(argv[++index]);
        else
            SMI_LEVEL_BOUND = atoi(&argv[index][10]);
        return true;
    }
*/
    if (strncmp(argv[index], "-prob-level", 11) == 0)
    {
        if (strcmp(argv[index], "-prob-level") == 0)
            PROBABILITY_PRUNE_LEVEL_START = atoi(argv[++index]);
        else
            PROBABILITY_PRUNE_LEVEL_START = atoi(&argv[index][11]);
        return true;
    }
    if (strncmp(argv[index], "-smi-only", 9) == 0)
    {
        Options::SMI_ONLY = true;
        return true;
    }
    if (strncmp(argv[index], "-serial", 7) == 0)
    {
        Options::THREADED = false;
        Options::SERIAL = true;
        return true;
    }
    if (strncmp(argv[index], "-threaded", 9) == 0)
    {
        Options::THREADED = true;
        Options::SERIAL = false;
        return true;
    }
    if (strncmp(argv[index], "-nopen", 6) == 0)
    {
        Options::OPENBABEL = false;
        return true;
    }
    if (strncmp(argv[index], "-pool", 5) == 0)
    {
        if (strcmp(argv[index], "-mw") == 0)
            OBGEN_THREAD_POOL_SIZE = atoi(argv[++index]);
        else
            OBGEN_THREAD_POOL_SIZE = atoi(&argv[index][5]);
        return true;
    }
    if (strncmp(argv[index], "-odir", 5) == 0)
    {
        if (strcmp(argv[index], "-odir") == 0)
            OUTPUT_DIR_SUFFIX = argv[++index];
        return true;
    }

    if (strncmp(argv[index], "-lip", 4) == 0)
    {
        if (strcmp(argv[index], "-lip") == 0)
            USE_LIPINSKI = true;
        return true;
    }
	
	if (strncmp(argv[index], "-freq", 4) == 0)
	{
		if (strcmp(argv[index], "-freq") == 0)
            FREQUENCY_ANALYSIS = true;
			IS_SUBJECT_FILE = false;
        return true;
	}
		
	if (strncmp(argv[index], "-dist", 4) == 0)
	{
		if (strcmp(argv[index], "-dist") == 0)
            DISTRIBUTION_ANALYSIS = true;
			IS_SUBJECT_FILE = true;
        return true;
	}

/*

    //
    // Check for an error: last command-line argument is an option (with no subsequent file)
    //
    if (index + 1 >= argc)
    {
        std::cerr << "Specified option " << argv[index]
                  << " not followed by the option value." << std::endl;
        return false;
    }
*/

    return false;
}

