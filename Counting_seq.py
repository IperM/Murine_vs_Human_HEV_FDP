from Fastqdumprun_3 import *

'''
This it's a short script that it's in order to check how well does perform the different
counting system generated, in order to make it working on the best block, for several reasons
Should be performed step by step but we can make it directly the whole pipeline. In order to 
make more understanding when executing on the terminal we must put inside a file the output so 
we end up getting a file with the time stamps in seconds per run(in this case will only measure
the ammount of time spended on the build and the mapping)
'''


def run_Picard_Singles(input_Path, savepath, rm = True, rm_seq_dup = False):

    os.chdir(input_Path)

    #In this case Picard requires BAM files that also needs to be sorted in order to work as fast as
    #possible. (can also recieve SAM theoretically: https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)

#    files = file_list(input_Path)
#    print(files)

#    for file in files:
#        if file == "Workfiles":
#            continue
        #Remove the whole set of unmapped reads, then sort the file and remove the 
#        new_file = Sam_unmap(file)
#        Sam_sorting(new_file, BAM= True)
#        os.remove(input_Path+"/"+new_file)
        

#    move_files(input_Path, savepath+'/Workfiles/', '.bam')
    
    files = file_list(savepath+'/Workfiles')
    print("SAM FILES READY:")
    print(files)

    Previous_Files = file_list(savepath)
    print(Previous_Files)

    os.chdir(savepath+'/Workfiles')

    if rm == True and rm_seq_dup == False:
        #In this case we will only move one folder before not as previously, as we want to keep only the sorted and maped
        #reads so we can just avoid the whole set of repeated files and the overfilling of memory.

        for file in files:
            file_check = "../marked_dup_"+ file[:-14] +".bam"
            if file_check in Previous_Files:
                print("Not missing anything :)")
                print(file)
                continue
            
            else:
                command = "picard MarkDuplicates REMOVE_DUPLICATES=true I={} O={} ".format(file, "../marked_dup_"+ file[:-14] +".bam")
                print(command)
                print(file)

                runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    elif rm_seq_dup == True and rm == False:
        for file in files:
            file_check = "../marked_dup_"+ file[:-14] +".bam"
            if file_check in Previous_Files:
                print("Not missing anything :)")
                print(file)
                continue
                
            else:
                command = "picard MarkDuplicates REMOVE_SEQUENCING_DUPLICATES=true I={} O={}".format(file, "../marked_dup_"+ file[:-14] +".bam")
                print(command)
                print(file)

                runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    elif rm_seq_dup == True and rm == True:   
        for file in files:
            file_check = "../marked_dup_"+ file[:-14] +".bam"
            if file_check in Previous_Files:
                print("Not missing anything :)")
                print(file)
                continue
            else:
                    
                command = "picard MarkDuplicates REMOVE_SEQUENCING_DUPLICATES=true REMOVE_DUPLICATES=true I={} O={} ".format(file, "../marked_dup_"+ file[:-14] +".bam")
                print(command)
                print(file)

                runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    else:
        for file in files:
            file_check = "../marked_dup_"+ file[:-14] +".bam"
            if file_check in Previous_Files:
                print("Not missing anything :)")
                print(file)
                continue

            else:

                command = "picard MarkDuplicates I={} O={} ".format(file, "../marked_dup_"+ file[:-14] +".bam")
                print(command)
                print(file)

                runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    return savepath


def run_Htseq_Singles(input_Path, savepath, annotated):

    files = file_list(input_Path)
    Previous_Files = file_list(savepath)
    print(Previous_Files)
    os.chdir(input_Path)

    for file in files:

        if file == "Workfiles":
            continue

        elif file in savepath:
            print("Not missing anything :)")
            print(file)
            continue

        #Dont accept CRAM files,
        #Generates the complete set of files in order to indexing
        #Then to work the overalll  we will be focusing on the generation of the indexed woy.
        #Sam_index(file)

        #You must be carefull when using the tool, you will have to set it up the most accurate to the real tools
        command = "htseq-count -f {} -r pos -t exon {} {} > {}".format( file[-3:],  file, annotated , savepath + '/' + file[:-4] + '_counts.txt')
        print(command)
        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)


    #Inside of this will become moving the files to propper places in order to not miss the whole errors
    #So we will be able inside the knowledge and start working with the different files.            
    move_files(input_Path, savepath, ".txt")
    move_files(input_Path, savepath + "/Workfiles", ".bai")
    move_files(input_Path, savepath + "/Workfiles", ".bam")


    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    return savepath

#Just a selection of the patterns that wants to become separated so it will move the files to the propper folders, desired
def file_detector(path, Savepath_1, Savepath_2, Savepath_3 = "", elefant = "HepG", mice_1 = "MLT", mice_2 = "56D"):

    files = file_list(path)
    unfound_files = []
    for file in files:

        index_1 = file.find(elefant)
        index_2 = file.find(mice_1)
        index_3 = file.find(mice_2)

        if index_1 != -1:
            #Detected as a human file, just a joke as being an elefant, but still usable for entering 
            #just move the file and continue
            move_files(path, Savepath_1, file)

        elif index_2 != -1 or index_3 != -1:
            #Detected as a murine file, move the file into the propper folder
            if Savepath_3 != "":
                move_files(path, Savepath_3, file)
            else:
                move_files(path, Savepath_2, file)

        else:
            #Unable to locate where it belongs, feels like home!
            #just keeps the in the folder and further on printed to screen
            unfound_files.append(file)

    print("This are your lost childs behave well with them next time:")
    print(unfound_files)
    return unfound_files

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

#DMP = Fastqdump(input_path= "/mnt/Viro_NGS_Data/raw_files/HEV_RNAseqNicola/2024-03_Viro_Nicola_Richard_01/2024-03_Viro_Nicola_Richard_01", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Fastqdump", old=False, paired=True)
#FAS = run_Fastqc(input_Path= "/mnt/Viro_NGS_Data/raw_files/HEV_RNAseqNicola/2024-03_Viro_Nicola_Richard_01/2024-03_Viro_Nicola_Richard_01", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Fastqc1Res")
#MUL = run_Multiqc(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Fastqc1Res", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/MultiqcRes")


'''Check the trimmed versions, separing the ones that are Paired, Unpaired, etc, then we will sort and work with them'''

#TRIM = run_Trimmomatic(input_Path= "/mnt/Viro_NGS_Data/raw_files/HEV_RNAseqNicola/2024-03_Viro_Nicola_Richard_01/2024-03_Viro_Nicola_Richard_01", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes", adapter_file = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE/adapter_found.fasta" )
#FAS2 =run_Fastqc(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Fastqc2Res")
#MUL2 = run_Multiqc(input_Path= FAS2, savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/MultiqcRes")


#============================================================================================================================#
#                                              TESTING MAPPING TOOLS                                                         #
#============================================================================================================================#


#BE CAREFULL, THIS EXPECTS 3 CHARACTERS IN THE BEFORE THE DOT; PROVIDING POSSIBLE ERRORS
#BUILD = run_HiSat2_build(input_Path="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE", reference= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa", simple= True )
#BUILD2= run_HiSat2_build(input_Path="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HUMAN_REFERENCE", reference= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HUMAN_REFERENCE/Homo_sapiens.GRCh38.dna.primary_assembly.fa", simple= True, fna = False)

#file_detector(path = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2", Savepath_1 = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Human_samples", Savepath_2="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Murine_samples")

#HS2_1 = run_HiSat2_v2(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Human_samples", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HisatRes/Human",
#                  reference= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HUMAN_REFERENCE/Homo_sapiens.GRCh38.dna.primary_assembly")
#HS2_2 = run_HiSat2_v2(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Murine_samples", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HisatRes/Murine",
#                  reference= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE/Mus_musculus.GRCm39.dna_sm.primary_assembly")


#PIC = run_Picard_v2(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HisatRes/Human" , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Picard/Human", rm= True, rm_seq_dup= False)
#HTS =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Picard/Human", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Htseq/Human", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HUMAN_REFERENCE/Homo_sapiens.GRCh38.111.gtf")

#PIC_2 = run_Picard_v2(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HisatRes/Murine" , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Picard/Murine", rm= True, rm_seq_dup= False)
#HTS_2 =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Picard/Murine", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Htseq/Murine", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE/Mus_musculus.GRCm39.111.gtf")

#============================================================================================================================#
#                                            DOING WITHOUT ANY TRIMM                                                         #
#============================================================================================================================#

#/mnt/Viro_NGS_Data/raw_files/HEV_RNAseqNicola/2024-03_Viro_Nicola_Richard_01/2024-03_Viro_Nicola_Richard_01

#BUILD = run_HiSat2_build(input_Path="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE", reference= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE/Mus_musculus.GRCm39.dna_sm.primary_assembly.fa", simple= True )
#BUILD2= run_HiSat2_build(input_Path="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HUMAN_REFERENCE", reference= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HUMAN_REFERENCE/Homo_sapiens.GRCh38.dna.primary_assembly.fa", simple= True, fna = False)

#file_detector(path = "/mnt/Viro_NGS_Data/raw_files/HEV_RNAseqNicola/2024-03_Viro_Nicola_Richard_01/2024-03_Viro_Nicola_Richard_01", Savepath_1 = "/mnt/Viro_NGS_Data/raw_files/HEV_RNAseqNicola/2024-03_Viro_Nicola_Richard_01/2024-03_Viro_Nicola_Richard_01/Human", Savepath_2="/mnt/Viro_NGS_Data/raw_files/HEV_RNAseqNicola/2024-03_Viro_Nicola_Richard_01/2024-03_Viro_Nicola_Richard_01/Murine")

#HS2_1 = run_HiSat2_v2(input_Path= "/mnt/Viro_NGS_Data/raw_files/HEV_RNAseqNicola/2024-03_Viro_Nicola_Richard_01/2024-03_Viro_Nicola_Richard_01/Human", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Untrimmed_procedure/HisatRes/Human",
#                  reference= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HUMAN_REFERENCE/Homo_sapiens.GRCh38.dna.primary_assembly")
#HS2_2 = run_HiSat2_v2(input_Path= "/mnt/Viro_NGS_Data/raw_files/HEV_RNAseqNicola/2024-03_Viro_Nicola_Richard_01/2024-03_Viro_Nicola_Richard_01/Murine", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Untrimmed_procedure/HisatRes/Murine",
#                  reference= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE/Mus_musculus.GRCm39.dna_sm.primary_assembly")


#PIC = run_Picard_v2(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Untrimmed_procedure/HisatRes/Human" , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Untrimmed_procedure/Picard/Human", rm= True, rm_seq_dup= False)
#HTS =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Untrimmed_procedure/Picard/Human", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Untrimmed_procedure/Htseq/Human", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HUMAN_REFERENCE/Homo_sapiens.GRCh38.111.gtf")

#PIC_2 = run_Picard_v2(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Untrimmed_procedure/HisatRes/Murine" , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Untrimmed_procedure/Picard/Murine", rm= True, rm_seq_dup= False)
#HTS_2 =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Untrimmed_procedure/Picard/Murine", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Untrimmed_procedure/Htseq/Murine", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE/Mus_musculus.GRCm39.111.gtf")

#============================================================================================================================#
#                                            HUMAN AGAINST MURINE                                                            #
#============================================================================================================================#

HS2_1 = run_HiSat2_v2(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Human_samples", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HisatRes/HumanVSMurine",
                 reference= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE/Mus_musculus.GRCm39.dna_sm.primary_assembly")


PIC = run_Picard_v2(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/HisatRes/HumanVSMurine" , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Untrimmed_procedure/Picard/HumanVSMurine", rm= True, rm_seq_dup= False)
HTS =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Untrimmed_procedure/Picard/HumanVSMurine", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Untrimmed_procedure/Htseq/HumanVSMurine", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE/Mus_musculus.GRCm39.111.gtf")

#============================================================================================================================#
#                                                  AGAINST P6 GENOME                                                         #
#============================================================================================================================#

reference = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE_P6/p6_genome.fasta"
MURINE_FILES = file_list("/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Murine_samples")
HUMAN_FILES = file_list("/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Human_samples")
Mapping_Save_folder = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/Mapping_p6_genome"

#run_HiSat2_build(input_Path= "/mnt/Viro_Data/Mitarbeiter/Leyla/HEV_PPH/Dataset/p6HEV_reference_genome", reference= reference, simple=True)

Murine = run_HiSat2_v2(input_Path = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Murine_samples", savepath = Mapping_Save_folder+"/Murine_p6", reference = reference)

Human = run_HiSat2_v2(input_Path = "/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/TrimmRes2/Human_samples", savepath= Mapping_Save_folder+"/Human_p6", reference= reference)

PIC_Murine = run_Picard(input_Path= Mapping_Save_folder+"/Murine_p6" , savepath= Mapping_Save_folder+"/Murine_Picard_p6", rm= True, rm_seq_dup= False)
PIC_Human = run_Picard(input_Path= Mapping_Save_folder+"/Human_p6" , savepath= Mapping_Save_folder+"/Human_Picard_p6", rm= True, rm_seq_dup= False)

HTS_Murine =  run_Htseq(input_Path= PIC_Murine, savepath= Mapping_Save_folder+"/Murine_HTseq_p6", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE_P6/p6_genome.gff3")
HTS_Human=  run_Htseq(input_Path= PIC_Human, savepath= Mapping_Save_folder+"/Human_HTseq_p6", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/Murine_FDP/REFERENCE_P6/p6_genome.gff3")
