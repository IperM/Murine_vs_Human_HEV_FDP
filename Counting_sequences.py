from Fastqdumprun import *

'''
This it's a short script that it's in order to check how well does perform the different
counting system generated, in order to make it working on the best block, for several reasons
Should be performed step by step but we can make it directly the whole pipeline. In order to 
make more understanding when executing on the terminal we must put inside a file the output so 
we end up getting a file with the time stamps in seconds per run(in this case will only measure
the ammount of time spended on the build and the mapping)
'''

#============================================================================================================================#
#                                                   PREPARING THE DATA                                                       #
#============================================================================================================================#

'''
Example of use, in order to use this you can just import the functions inside the other scripts or directly
apply the usage . Another options could be ust changing the paths of the bottom functions, this will create the files and 
store them in the different folders, notice that some of them will be requiring to add a Workfile generation of the propper files.

Download the full set of data on the 6 different states for the control and the afected samples 
See its initial behaviour on the overall sense. 
'''

#DMP = Fastqdump(input_path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/FasqDumpRes/Workfiles", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/FasqDumpRes", old=False, paired=True)
#FAS = run_Fastqc(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/FasqDumpRes", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/Fastqc1Res")
#MUL = run_Multiqc(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/Fastqc1Res", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/MultiqcRes")


'''Check the trimmed versions, separing the ones that are Paired, Unpaired, etc, then we will sort and work with them'''

#TRIM = run_Trimmomatic(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/FasqDumpRes", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes" )
#FAS2 =run_Fastqc(input_Path= TRIM, savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/Fastqc2Res")
#MUL2 = run_Multiqc(input_Path= FAS2, savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/MultiqcRes")


#============================================================================================================================#
#                                              TESTING MAPPING TOOLS                                                         #
#============================================================================================================================#

BUILD = run_HiSat2_build(input_Path="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles", reference= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/GCF_000001405.40_GRCh38.p14_genomic.fna", simple= True )
HS2 = run_HiSat2(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res",
21                  reference= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/GCF_000001405.40_GRCh38.p14_genomic" )
PIC = run_Picard(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res" , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", rm= True, rm_seq_dup= False)
HTS =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/genomic.gtf")
move_files("/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes","/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes/HiSat2", "counts")

'''
We will start with running over STAR, so we can do the most direct comparation (consider that we're indexing )
Just before inside the run_STAR function.
'''

#STAR = run_STAR(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/STARRes",
#                  reference= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/GCF_000001405.40_GRCh38.p14_genomic.fna",  gtf_file= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/genomic.gtf")
#PIC = run_Picard(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/STARRes" , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", rm= True, rm_seq_dup= False)
#HTS =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/genomic.gtf")
#move_files("/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes","/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes/STAR", "_counts")


''' Over this set of different files we won't need at first the propper annotation file of the genomic sequence, but further on it will be used.'''

#Bow = run_Bowtie(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/BowtieRes",
#                  reference= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/GCF_000001405.40_GRCh38.p14_genomic.fna")
#PIC = run_Picard(input_Path= Bow , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", rm= True, rm_seq_dup= False)
#HTS =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/genomic.gtf")
#move_files("/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes","/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes/Bowtie2", "counts")


#BWA = run_BWA(input_Path="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/BowtieRes",
#              reference="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/GCF_000001405.40_GRCh38.p14_genomic.fna")
#PIC = run_Picard(input_Path= BWA , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", rm= True, rm_seq_dup= False)
#HTS =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/genomic.gtf")
#move_files("/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes","/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes/BWA", "counts")


#Kal = run_Kallisto(input_Path="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/BowtieRes",
#              reference="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/GCF_000001405.40_GRCh38.p14_genomic.fna", gtf_file = "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/genomic.gtf")
#PIC = run_Picard(input_Path= BWA , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", rm= True, rm_seq_dup= False)
#HTS =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/genomic.gtf")
#move_files("/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes","/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes/Kallisto", "counts")


#Nov = run_NovoAlign(input_Path="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/NovoRes",
#              reference="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/GCF_000001405.40_GRCh38.p14_genomic.fna")
#PIC = run_Picard(input_Path= Nov , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", rm= True, rm_seq_dup= False)
#HTS =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/genomic.gtf")
#move_files("/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes","/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes/NovoAllign", "counts")


#============================================================================================================================#
#                                       TESTING SINGLE ENVIROMENT MAPPING TOOLS                                              #
#============================================================================================================================#

'''
This set of tools have different errors inside compatibilities so we will avoid this problems
 by generating a new enviroment with just the necessay information
'''


#SOAP2 = run_SOAP2(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/SOAP2Res",
#                  reference= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/GCF_000001405.40_GRCh38.p14_genomic.fna")              
#PIC = run_Picard(input_Path= SOAP2 , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", rm= True, rm_seq_dup= False)
#HTS =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/genomic.gtf")
#move_files("/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes","/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes/SOAP2", "counts")


#MPS = run_MapSplicer(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/SOAP2Res",
#                  reference= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/GCF_000001405.40_GRCh38.p14_genomic.fna")
#PIC = run_Picard(input_Path= MPS , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", rm= True, rm_seq_dup= False)
#HTS =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/genomic.gtf")
#move_files("/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes","/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes/MapSplicer", "counts")
