# Murine vs Human HEV FDP

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

-----
### Main Functions

In the following steps I will explain what some of the key functions from **Fastqdumprun_3.py** are important
or develop the function in order to know  how it is used:

- Fastqdump: extract the information from the GEO database for obtaining the different samples.
- run_Fastqc: Performs the basic analyisis for quality control for the whole folder.
- run_Multiqc: Perform analysis on multi, all together joint inside a single file.
- run_Trimmomatic: Does the _Trimmommatic_ execution, inside it has base functions, but apart has a flag (*SE_K*)
  that is used for change and do it for single ended files (if seted up on True).
- run_Hisat_extract: Extract the relevant exons or splice sites depending if the flag exon it's setted to
  True or False.
- run_Hisat_build: Function that generates the reference file as it is required for doing the acual alignement.
  It has a flag, one for doing the basic execution (simple) which will make the execution much shorter but still
  good, then if it's not selected as simple the exons and splice sites can be setted down which will help to define a more
  accurate reference.
- run_Hisat_v2: Runs the hisat alignement, with the reference previously build, having a SE flag for single end files or not. Using
  the new files for the murine samples
- run_Hisat: Works as the previous one but it has another source, but it is configured for files with different names, for example previous
  project.
- run_Hisat_old: Same as the previous before but the SE flag does not exist.
