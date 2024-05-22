# Murine_vs_Human_HEV_FDP

This Repository comprehends the different files that mount the whole pipeline. Inside this it is located 2 different folders,
containing the appropiate **scripts**.

## Folders:

### Old_Version:
Contains inside the previous version of the files, used for other project that was performed. This were used in order to generate
counts of _Pig and human_ data. It has less versions for the different functions.

Essentially the different documents inside here does the same as in the newer version but changing some core functionality.

### Murine_FDP:
Contains the actual files and how does it works in the sake of analysis of Murine and Human samples.

The files are the following ones:

- Fastqdumprun_3: Source script containing the majority of functions related to the analysis, being python script plus calls on bash.
- FastTrimm: This is a script for doing a quick cut of the sequences and generate a quality control for troubleshooting.
- Counting_seq: Some troubleshooting, doing instead of the whole pipeline, only one step at a time. Then following on it has the
  rest of the test, trying to apply the pipeline for doing countng check, arriving to the aounts of the genes.
- CovergeCheck_v2.py: Script that performs the coverage of the different samples in respect to the different reference genomes, this
  contains several examples to do it in different ways.


In the following steps I will explain what some of the key functions from **Fastqdumprun_3.py** are important
or develop the function in order to know  how it is used:

