from Fastqdumprun_3 import *

#BWA_Murine = run_BWA_v2(input_Path  = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Murine_samples", savepath ="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/BWA/Murine",
#    reference = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa")

#BWA_Human = run_BWA_v2(input_Path="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Human_samples", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/BWA/Human", 
#    reference= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HUMAN_REFERENCE/Homo_sapiens.GRCh38.dna.primary_assembly.fa")

work_folder = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/BWA"
save_folder = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Unmaped_Reads"
reference = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE_P6/p6_genome.fasta"
os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/")

#for file in file_list(work_folder+"/Human"):
#    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/BWA/Human")
#    Sam_unmap(file, BAM = False, Keep_Unmaped = True)
#    move_files(work_folder+"/Human", save_folder+"/HUMAN", '.bam')

#for file in file_list(save_folder+"/HUMAN"):
#    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Unmaped_Reads/HUMAN")
#    Sam_To_Fastq(file)
#    move_files(save_folder+"/HUMAN",save_folder+"/HUMAN/FASTQ" , '.fastq')



#for file in file_list(work_folder+"/Murine"):
#    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/BWA/Murine")
#    Sam_unmap(file, BAM = False, Keep_Unmaped = True)
#    move_files(work_folder+"/Murine", save_folder+"/MURINE", '.bam')

#for file in file_list(save_folder+"/MURINE"):
#    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Unmaped_Reads/MURINE")
#    Sam_To_Fastq(file)
#    move_files(save_folder+"/MURINE",save_folder+"/MURINE/FASTQ" , '.fastq')

#os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/")

#run_Tanoti_v2(input_Path = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Unmaped_Reads/HUMAN/FASTQ", savepath = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Unmaped_Reads/Tanoti/Human", reference = reference)
#run_Tanoti_v2(input_Path = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Unmaped_Reads/MURINE/FASTQ", savepath = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Unmaped_Reads/Tanoti/Murine", reference = reference)

#PIC = run_Picard_v2(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Unmaped_Reads/Tanoti/Human", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Unmaped_Reads/Picard/Human", rm= True, rm_seq_dup= False)
#PIC_2 = run_Picard_v2(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Unmaped_Reads/Tanoti/Murine" , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Unmaped_Reads/Picard/Murine", rm= True, rm_seq_dup= False)

os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Unmaped_Reads/Picard/Murine")

for file in file_list("/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Unmaped_Reads/Picard/Murine"):
    to_move = Sam_Coverage(file)
    Sam_Summary(file)

os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Unmaped_Reads/Picard/Human")

for file in file_list("/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Unmaped_Reads/Picard/Human"):
    to_move = Sam_Coverage(file)
    Sam_Summary(file)