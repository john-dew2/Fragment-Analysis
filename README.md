# Fragment Analysis

Fragment Analysis is a analytic tool that conducts TC comparisons on molecular fragments.
This program branches heavily off of the work of eSynth, and uses its fragment processing functions to generate the sets used in analyzation. If you find this tool to be useful, be sure to cite the creators of eSynth:

Naderi, Misagh, Chris Alvin, Yun Ding, Supratik Mukhopadhyay, and Michal Brylinski. "A graph-based approach to construct target-focused libraries for virtual screening." Journal of Cheminformatics 8, no. 1 (2016): 14.

This README file is written by Johnathan Dewey. Last Update: 06/29/2022

# Requirements
  - GSL and ZLIB
  - OpenBabel 2.4.1 or newer (OpenBabel is no longer supported and wont recieve updates)
  - C++
  - Makefile
  
# Usage
  1. First, OpenBabel must be downloaded, this take approximately 11 minutes, so don't get worried if it takes longer. The Openbabel Download script is located in /Fragment-Analysis/src
  2. Then clone this github and run these commands:
      - !mkdir Fragment-Analysis/src/obj
      - !make -C /content/Fragment-Analysis/src/
      - !rm -r /content/Fragment-Analysis/src/synth_output_dir
  3. Enter input
      - %cd /content/Fragment-Analysis/src/ (cd the source file so we can access any input files)
      - ! ./efreq -nopen -serial -smi-only [types of analysis] -prob-level 2 [input file names]
      -   use tag "-freq" to use frequency analysis (frequency analysis takes sets of fragments)
      -   use tag "-dist" to use distribution analysis (note that distribution analysis takes an input file of one molecule first, and sets of fragments)

# Analysis
  1. Frequency Analysis: Given brick and linker fragments, this analysis takes every fragment and runs a Tanimoto Coefficient (TC) comparison on every fragment in the input. If the TC value reveals a similarity between the two fragments, then the relationship is recorded
      - Format:  [File Name].sdf: [Number of Fragments Similar],...
      - Example: b-CHEMBL175476.mol2-001.sdf: 1 fragments similar,...
  2. Distribution Analysis: Given a fragment, named the subject, and a set of brick or linker fragments, this analysis will run TC comparisons on the given set by comparing each fragment to the subject. This analysis will tell you how similar a specific fragment is to a set of fragments.
      - Format: Tanimoto Coefficient: [TC Group Value (0.0-1.0)], contains [Number of Fragments with TC Group Value]...
      - Example: Tanimoto Coefficient: 0.1, contains 4 fragments...
# Output
  All output is done to the directory /Fragment-Analysis/src/
  The types of files:
  - brick-frequency-report.txt     | Reports the outcome of a frequency analysis on a set of bricks
  - linker-frequency-report.txt    | Reports the outcome of a frequency analysis on a set of linkers
  - brick-distribution-report.txt  | Reports the outcome of a distribution analysis on a set of bricks
  - linker-distribution-report.txt | Reports the outcome of a distribution analysis on a set of linkers
