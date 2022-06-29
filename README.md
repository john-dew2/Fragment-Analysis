# Fragment Analysis

Fragment Analysis is a analytic tool that conducts TC comparisons on molecular fragments.
This program branches heavily off of the work of eSynth, and uses its fragment processing functions to generate the sets used in analyzation. If you find this tool useful be sure to cite Dr. Chris Alvin and Computational Systems Biology Group:

Naderi, Misagh, Chris Alvin, Yun Ding, Supratik Mukhopadhyay, and Michal Brylinski. "A graph-based approach to construct target-focused libraries for virtual screening." Journal of Cheminformatics 8, no. 1 (2016): 14.

This README file is written by Johnathan Dewey. Last Update: 06/29/2022

# Requirements
  - OpenBabel 2.3.1 or newer (OpenBabel is no longer supported and wont recieve updates)
  - C++
  - Makefile
  
# Usage
  1. First, OpenBabel must be downloaded, this take approximately 11 minutes, so don't get worried if it takes longer
  2. Then clone this github and run these commands:
      - !mkdir Fragment-Analysis/src/obj
      - !make -C /content/Fragment-Analysis/src/
      - !rm -r /content/Fragment-Analysis/src/synth_output_dir
  3. Enter input
      - %cd /content/Fragment-Analysis/src/ (cd the source file so we can access any input files)
      - ! ./efreq -nopen -serial -smi-only [types of analysis] -prob-level 2 [input file names]
            - use "-freq" to use frequency analysis
            - use "-dist" to use distribution analysis
# Output
  All output is done to the directory /Fragment-Analysis/src/
