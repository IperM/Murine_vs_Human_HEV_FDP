from Fastqdumprun import *


'''
In this case we will use HiSat and BWA in order to find the mapping of the different 
files ending up with the coverage of each sample with the p6 HEV genome. So we will 
know how much infected genome have inside. 

This will be done through different loops, one in PPH data and another in PHH data.

'''
reference = "/mnt/Viro_Data/Mitarbeiter/Leyla/HEV_PPH/Dataset/p6HEV_reference_genome/p6HEV_genome.fasta"
PPH_FILES = file_list("/mnt/Viro_Data/Mitarbeiter/Leyla/HEV_PPH/Results/Trimmed_Fasta")
PHH_FILES = file_list("/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes")
Mapping_Save_folder = "/mnt/Viro_Data/Mitarbeiter/Leyla/HEV_PPH/Results/Mapping_p6_genome"

print("PPH Files:"+"\n"+PPH_FILES)
print("===================================")
print("PHH Files:"+"\n"+PHH_FILES)


run_HiSat2_build(input_Path= "/mnt/Viro_Data/Mitarbeiter/Leyla/HEV_PPH/Dataset/p6HEV_reference_genome", reference= reference, simple=True)

run_HiSat2(input_Path = "/mnt/Viro_Data/Mitarbeiter/Leyla/HEV_PPH/Results/Trimmed_Fasta", savepath = Mapping_Save_folder+"/PPH", reference = reference)
print("Execution finished - check the path")

run_HiSat2(input_Path="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes", savepath= Mapping_Save_folder+"/PHH", reference= reference)
print("Execution finished - check the path")
