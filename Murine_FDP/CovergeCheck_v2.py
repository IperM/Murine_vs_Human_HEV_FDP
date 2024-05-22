from Fastqdumprun_3 import *


'''
In this case we will use HiSat and BWA in order to find the mapping of the different 
files ending up with the coverage of each sample with the p6 HEV genome. So we will 
know how much infected genome have inside. 

This will be done through different loops, one in PPH data and another in PHH data.

'''

reference = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE_P6/p6_genome.fasta"
MURINE_FILES = file_list("/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Murine_samples")
HUMAN_FILES = file_list("/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Human_samples")
Mapping_Save_folder = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Mapping_p6_genome"

print("MURINE Files:"+"\n")
print(MURINE_FILES)
print("===================================")
print("Human Files:"+"\n")
print(HUMAN_FILES)

### Execution of the building the new HiSat reference system for both samples, human and pig.
### Does not to be re-generated cause it has already done on previous steps.
#run_HiSat2_build(input_Path= "/mnt/Viro_Data/Mitarbeiter/Leyla/HEV_PPH/Dataset/p6HEV_reference_genome", reference= reference, simple=True)


print("=====================================================")

###Here starts the preparation for looking after the coverage on the tools, by trimming and mapping letting the data ready for further executions.

#Murine = run_HiSat2_v2(input_Path = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Murine_samples", savepath = Mapping_Save_folder+"/Murine_p6", reference = reference)
#print("Execution finished - check the path")

#Human = run_HiSat2_v2(input_Path="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Human_samples", savepath= Mapping_Save_folder+"/Human_p6", reference= reference)
#print("Execution finished - check the path")

#PIC_Murine = run_Picard(input_Path= Mapping_Save_folder+"/Murine_p6" , savepath= Mapping_Save_folder+"/Murine_Picard_p6", rm= True, rm_seq_dup= False)
#print("Picard Murine - To check")

#PIC_Human = run_Picard(input_Path= Mapping_Save_folder+"/Human_p6" , savepath= Mapping_Save_folder+"/Human_Picard_p6", rm= True, rm_seq_dup= False)
#print("Picard Human - To check")

#run_rm(PIC_Murine, prefix="marked_dup_metrics_")
#run_rm(PIC_Human, prefix="marked_dup_metrics_")


print("=====================================================")

os.chdir(Mapping_Save_folder+"/Murine_Picard_p6")

for file in file_list(Mapping_Save_folder+"/Murine_Picard_p6"):
    to_move = Sam_Coverage(file)
    Sam_Summary(file)

os.chdir(Mapping_Save_folder+"/Human_Picard_p6")

for file in file_list(Mapping_Save_folder+"/Human_Picard_p6"):
    to_move = Sam_Coverage(file)
    Sam_Summary(file)

print("txt files generated, check to see further information")

print("=====================================================")

os.chdir("/mnt/Viro_Data/Mitarbeiter/Leyla/HEV_PPH")

###Afterwards we will have to run picard and the gene counting over
###the different samples to  then make some plots on further steps if required ( normally it wont be
###  required for basic analysis)

#HTS_PPH =  run_Htseq(input_Path= PIC_PPH, savepath= Mapping_Save_folder+"/Murine_HTseq_p6", annotated="/mnt/Viro_Data/Mitarbeiter/Leyla/HEV_PPH/Dataset/p6HEV_reference_genome/annotation_p6_genome.gff3")
#print("Counts PPH - To check")

#HTS_PHH =  run_Htseq(input_Path= PIC_PHH, savepath= Mapping_Save_folder+"/Human_HTseq_p6", annotated="/mnt/Viro_Data/Mitarbeiter/Leyla/HEV_PPH/Dataset/p6HEV_reference_genome/annotation_p6_genome.gff3")
#print("Counts PHH - To check")


#================================================================================================#
### Now it'S the turn for BWA which it'S one of the most widely used mappers, but also has high detection of every count.

BWA_Murine = run_BWA_v2(input_Path  = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Murine_samples", savepath = Mapping_Save_folder+"/Murine_BWA_p6", reference = reference)
print("Execution finished - check the path")

BWA_Human = run_BWA_v2(input_Path="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Human_samples", savepath= Mapping_Save_folder+"/Human_BWA_p6", reference= reference)
print("Execution finished - check the path")

PIC_BWA_Murine = run_Picard(input_Path= Mapping_Save_folder+"/Murine_BWA", savepath= Mapping_Save_folder+"/Murine_Picard_BWA_p6", rm= True, rm_seq_dup= False)
print("Picard PPH - To check")

PIC_BWA_Human = run_Picard(input_Path= Mapping_Save_folder+"/Human_BWA" , savepath= Mapping_Save_folder+"/Human_Picard_BWA_p6", rm= True, rm_seq_dup= False)
print("Picard PHH - To check")

run_rm(PIC_BWA_Murine, prefix="marked_dup_metrics_")
run_rm(PIC_BWA_Human, prefix="marked_dup_metrics_")

print("=====================================================")

os.chdir(Mapping_Save_folder+"/PPH_Picard_BWA_p6")

for file in file_list(Mapping_Save_folder+"/PPH_Picard_BWA_p6"):
    to_move = Sam_Coverage(file)
    Sam_Summary(file)

os.chdir(Mapping_Save_folder+"/PHH_Picard_BWA_p6")

for file in file_list(Mapping_Save_folder+"/PHH_Picard_BWA_p6"):
    to_move = Sam_Coverage(file)
    Sam_Summary(file)


#================================================================================================#
### Few steps into it, we can start with the data of Tanoti, which it'S specially designed  to treat with viral data

print("=====================================================")

###Here starts the preparation for looking after the coverage on the tools, by trimming and mapping letting the data ready for further executions.

Murine_tan = run_Tanoti_v2(input_Path = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Murine_samples", savepath = Mapping_Save_folder+"/Murine_tanoti_p6", reference = reference)
print("Execution finished - check the path")

Human_tan = run_Tanoti_v2(input_Path="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Human_samples", savepath= Mapping_Save_folder+"/Human_tanoti_p6", reference= reference)
print("Execution finished - check the path")

PIC_tan_Murine = run_Picard(input_Path= Mapping_Save_folder+"/Murine_tanoti" , savepath= Mapping_Save_folder+"/Murine_Picard_tanoti_p6", rm= True, rm_seq_dup= False)
print("Picard Murine - To check")

PIC_tan_Human = run_Picard(input_Path= Mapping_Save_folder+"/Human_tanoti" , savepath= Mapping_Save_folder+"/Human_Picard_tanoti_p6", rm= True, rm_seq_dup= False)
print("Picard Human - To check")

run_rm(PIC_tan_Murine, prefix="marked_dup_metrics_")
run_rm(PIC_tan_Human, prefix="marked_dup_metrics_")


print("=====================================================")

os.chdir(Mapping_Save_folder+"/Murine_tanoti_Picard_p6")

for file in file_list(Mapping_Save_folder+"/Murine_tanoti_Picard_p6"):
    to_move = Sam_Coverage(file)
    Sam_Summary(file)

os.chdir(Mapping_Save_folder+"/Human_tanoti_Picard_p6")

for file in file_list(Mapping_Save_folder+"/Human_tanoti_Picard_p6"):
    to_move = Sam_Coverage(file)
    Sam_Summary(file)

    
print("=====================================================")
