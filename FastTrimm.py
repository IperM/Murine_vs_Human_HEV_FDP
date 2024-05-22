from Fastqdumprun_3 import *

TRIM = run_Trimmomatic(input_Path= "/mnt/Viro_NGS_Data/raw_files/HEV_RNAseqNicola/2024-03_Viro_Nicola_Richard_01/2024-03_Viro_Nicola_Richard_01", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2", adapter_file = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Assistant_Functions/Adaptors_sumed/adapter_found.fasta" )
FAS2 =run_Fastqc(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Fastqc3Res")
