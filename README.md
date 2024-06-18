# Murine vs Human HEV FDP

This repository contains the files that make up the entire pipeline for analyzing murine and human HEV FDP samples. It is organized into two main folders, each containing the relevant scripts.

## Folders

### Old_Version
This folder contains the previous version of the files used in another project focused on generating counts for pig and human data. The functions in these files are similar to those in the newer version but with some core differences.

### Murine_FDP
This folder contains the current version of the files used for analyzing murine and human samples. The key scripts in this folder are:

- **Fastqdumprun_3.py**: This Python script, which includes several Bash calls, contains the majority of functions related to the analysis.
- **FastTrimm.py**: A script for trimming sequences and performing quality control.
- **Counting_seq.py**: A troubleshooting script that allows for step-by-step execution of the pipeline, focusing on generating gene counts.
- **CoverageCheck_v2.py**: A script that checks the coverage of samples against different reference genomes, with several examples provided.
- **FC.sh**: A Bash script for running the FeatureCounts tool on an entire folder.
- **retriver_adaptors.sh**: A Bash script that iterates over quality control files to extract overrepresented sequences and store them in a new file.
- **Analysis_Murine_Human_Filtered_v3.Rmd**: Rmd script containing every visualization, DGEA, and GO analysis for human and murine samples.
- **Coverage_Samples.Rmd**: Rmd script containing visualization of coverage plot coming from the mapping against p6 genome. Shows the coverage plots with an additional scheme of p6 genome down for recognizing where it's the coverage present.
- **Counting_seq_Feature.py**: A troubleshooting script that allows for step-by-step execution of the pipeline, focusing on generating gene counts, through the usage of FeatureCounts.

## Main Functions

Here are some key functions from **Fastqdumprun_3.py** and their descriptions:

- **Fastqdump**: Extracts information from the GEO database to obtain samples.
- **run_Fastqc**: Performs basic quality control analysis for the entire folder.
- **run_Multiqc**: Compiles multi-sample QC reports into a single file.
- **run_Trimmomatic**: Executes Trimmomatic for trimming sequences. It includes a flag (*SE_K*) to process single-end files if set to True.
- **run_Hisat_extract**: Extracts relevant exons or splice sites based on the *exon* flag.
- **run_Hisat_build**: Generates the reference file required for alignment. It includes a flag for simple execution and options for including exons and splice sites for a more accurate reference.
- **run_Hisat_v2**: Runs the HISAT alignment for murine samples using the generated reference, with an option for single-end files.
- **run_Hisat**: Similar to *run_Hisat_v2* but configured for files from previous projects.
- **run_Hisat_old**: An older version of the HISAT alignment function without the single-end file option.
- **Samtools Functions**: Various functions using Samtools for:
  - Sorting
  - Removing and retrieving unmapped reads
  - Converting between SAM and BAM, and BAM and FASTQ formats
  - Indexing aligned files to generate BAI files
  - Checking coverage and generating coverage plots
- **run_Picard_V2**: Runs Picard tools for analysis, with flags to remove duplicates and sequencing duplicates.
- **run_Picard**: Similar to _run_Picard_V2_ but with different sample naming conventions.
- **run_Picard_fast_new**: Performs the same analysis as the previous two functions but requires pre-sorted and filtered data.
- **run_HTseq**: Counts sequences using a _GTF/GFF_ file for reference.

In addition to HTseq and FeatureCounts, the pipeline supports multiple mappers, including BWA, HISAT, and Tanoti.
