import sys, os
import subprocess as sp
import shutil
from pathlib import Path
import time

#============================================================================================================================#
#                                                EXECUTION   FUNCTION                                                        #
#============================================================================================================================#


def file_list(path):
    listed = os.listdir(path)
    return(listed)


def move_files(source_folder, destination_folder, file_extension):
    # Ensure that both the source and destination folders exist
    if not os.path.exists(source_folder):
        print(f"Source folder '{source_folder}' does not exist.")
        return ''

    if not os.path.exists(destination_folder):
        print(f"Destination folder '{destination_folder}' does not exist. Creating it.")
        os.makedirs(destination_folder)

    # Get a list of all files in the source folder with the specified extension
    files_to_move = [file for file in os.listdir(source_folder) if file.endswith(file_extension)]

    # Move each file to the destination folder
    for file_name in files_to_move:
        source_path = os.path.join(source_folder, file_name)
        destination_path = os.path.join(destination_folder, file_name)
        shutil.move(source_path, destination_path)
        print(f"Moved '{file_name}' from '{source_folder}' to '{destination_folder}'.")

    return ''


def Fastqdump(input_path, savepath, paired = False, old = True):
    

    os.chdir(input_path)


    for file in file_list(input_path):
        if old == True:
            command = "fastq-dump  --split-3 {} ".format(file)
            print(command)

            runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

        elif paired == True:
            items = sp.run("cat {}".format(file), shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)
            SRR = items.stdout.split()
            for SR in SRR:
                command = "fastq-dump --split-3 {} ".format(SR)
                print(command)

                runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

        else:
            items = sp.run("cat {}".format(file), shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)
            SRR = items.stdout.split()
            for SR in SRR:
                command = "fastq-dump {} ".format(file)
                print(command)

                runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    move_files(input_path, savepath, ".fastq")

    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")


    return savepath


def run_Fastqc(input_Path, savepath):
    os.chdir(input_Path)
    archives = file_list(input_Path)
    separator= ' ' 
    archive = separator.join(archives)

    command = ["fastqc --extract --delete {} -o {}".format( archive, savepath)]
    print(command)

    runs = sp.run(command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, text=True)

    out = runs.stdout.splitlines()

    # Returns the path where the actual fastqc runs are saved
    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    return savepath


def run_Multiqc(input_Path, savepath , options=''):
    os.chdir(input_Path)

    command = ["multiqc{} -o {} {} ".format(' '+ options if options else '', savepath, input_Path)]
    print(command)

    runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    # Returns the path where the actual fastqc runs are saved and changes the working directory to the source of the project
    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")  
    return savepath


def PathMatcher(path):
    #Aux Function
    # Get the list of files from the specified path
    files = file_list(path)
    grouped_files = {}

    files = files[1:]
    for file_name in files:

        identifier = file_name[:-11]  # Extract identifier

        if identifier not in grouped_files:
            grouped_files[identifier] = []
        grouped_files[identifier].append(file_name)

    # Create pairs by taking elements 2 by 2 from each group
    pairs = [grouped_files[identifier][i:i+2] for identifier in grouped_files for i in range(0, len(grouped_files[identifier]), 2)]

    return pairs


def run_Trimmomatic(input_Path , 
                    savepath,
                     adapter_file ="/mnt/Viro_Data/Mitarbeiter/Ian/a_module/adapters/NEBnextUltra.fa"):


    #Implement the ILLUMINACLIP:adapters   : must be the full length path

    #In this case as we're using Multiqc we can see that all are affected by Nextera adapters, and by looking after the 
    #supplementary information we can enclose and just delete the ones that are inside this file

    os.chdir(input_Path)
    input_Path_archives = PathMatcher(input_Path)

    
    for set_files in input_Path_archives:
        print(set_files)

        command = "trimmomatic PE -phred33 {} {} -baseout ../TrimmomaticRes/{}  ILLUMINACLIP:{}:4:30:10 ".format(set_files[1][:-7]+'1'+set_files[1][12:], set_files[0][:-7]+'2'+set_files[0][12:], set_files[0][:-8], adapter_file)

        print(command)
        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)


    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    return savepath

def Files_Mating(path):
    files = file_list(path)
    paired = []

    for file_name in files:
        if file_name.endswith('1P'):
            paired.append(file_name)

    return paired

def gunzip(input_Path, file):
    os.chdir(input_Path)

    command = "gunzip {} ".format(file)
    print(command)
    runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    return file[:-3]


def run_HiSat2_extract(input_Path, file, exon = False):
    os.chdir(input_Path)

    if exon == True:
        
        file = gunzip(input_Path, file)

        command = "hisat2_extract_exons.py {} > {}".format(file, file[:-4] + ".exon")
        print(command)
        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

        os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
        return file[:-4] + ".exon"
    
    else :

        file = gunzip(input_Path, file)

        command = "hisat2_extract_splice_sites.py {} > {}".format(file, file[:-4]+'.ss')
        print(command)
        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

        os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
        return file[:-4]+'.ss'


def run_HiSat2_build(input_Path, reference, ss="", exon="", simple = False):
    os.chdir(input_Path)
    start = time.time()
    #They need to be renamed as .fasta .fa and uncompressed in order to work

    if simple == True: 
        #Requires 6GB memmory (20 min)
        command = "hisat2-build -p 16 {} {}".format(reference, reference[:-3])
        print(command)
        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

        os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
        end = time.time()
        print("Time elapsed on the whole execution = {} sec".format(end - start))
        return input_Path + '/' + reference
    
    else:
        #Requires 160GB memory ( 1 h)
        command = "hisat2-build -p 16 --exon {} --ss {} {} {}".format(exon, ss, reference, reference[:-3])
        print(command)
        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

        os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
        end = time.time()
        print("Time elapsed on the whole execution = {} sec".format(end - start))
        return input_Path + '/' + reference
    


def run_HiSat2(input_Path, savepath, reference):
    '''
   Be carefull when dealing with this we will put only the name of the sentence than comes before the acual 
   indexing file termination. If you selet the referencing file, it will end up generating a mess .

   Also very sensitive when talking of file types.
    '''
    start = time.time()
    os.chdir(input_Path)

    input_Path_archives = Files_Mating(input_Path)
    print(input_Path_archives)
    for set_files in input_Path_archives:

        command = "hisat2 -x {} -1 {} -2 {} -S {}".format(reference, set_files, set_files[:-2]+"2P",  set_files[:-2] + 'P_HiSat.sam')
        print(command)

        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    move_files(input_Path, savepath, ".sam")

    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    end = time.time()
    print("Time elapsed on the whole execution = {} sec".format(end - start))
    return savepath


def Sam_to_Bam(file_to_convert):
    command = "samtools view -b {} -o {}.bam".format(file_to_convert, file_to_convert[:-4])
    print(command)

    runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    return file_to_convert[:-4]+'.bam'


def Bam_to_Sam(file_to_convert):
    command = "samtools view -h {} -o {}.sam".format(file_to_convert, file_to_convert[:-4])
    print(command)

    runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    return file_to_convert[:-4]+'.sam'


def Sam_unmap(file_to_convert, BAM = False):
     
    #F in order to keep the mapped
    #f in order to keep the unmaped

    if BAM == True:
        if file_to_convert[-4:] == ".bam":

            command = "samtools view -h -b -F 4 {}>{}_mapped.bam".format(file_to_convert, file_to_convert[:-4])
            print(command)

            runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    else:
        if file_to_convert[-4:] == ".sam":

            command = "samtools view -h -F 4 {}>{}_mapped.sam".format(file_to_convert, file_to_convert[:-4])
            print(command)

            runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    return file_to_convert[:-4]+"_mapped.sam"


def Sam_sorting(file_to_sort, BAM = False):
    
    #By default it is sorting by coords so we dont need to change add new parameters

    if BAM == True:
        command = "samtools sort {} -o {}_sorted.bam".format(file_to_sort, file_to_sort[:-4])
        print(command)
        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)
            
        return file_to_sort[:-4]+"_sorted.bam"

    else:
        if file_to_sort[-4:] == ".sam":
            command = "samtools sort  {} -o {}_sorted.sam".format(file_to_sort, file_to_sort[:-4])
            print(command)

            runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

            return file_to_sort[:-4]+"_sorted.sam"


def Sam_index(file_to_index):
    #It can recieve bot BAM or SAM files without any matter, it will generate the possible indexedd
    #files
    command = "samtools index {}".format(file_to_index)
    print(command)

    runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)
    
    return


def run_rm(folder_path, files_to_remove):
    for file_name in files_to_remove:
        file_path = os.path.join(folder_path, file_name)
        try:
            os.remove(file_path)
            print(f"File {file_name} removed successfully.")
        except FileNotFoundError:
            print(f"File {file_name} not found.")
        except Exception as e:
            print(f"Error removing file {file_name}: {e}")

def run_Picard(input_Path, savepath, rm = True, rm_seq_dup = False):

    os.chdir(input_Path)

    #In this case Picard requires BAM files that also needs to be sorted in order to work as fast as
    #possible. (can also recieve SAM theoretically: https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)

    files = file_list(input_Path)[1:]
    print(files)
    for file in files:
        if file == "Workfiles":
            continue
        #Remove the whole set of unmapped reads, then sort the file and remove the 
        new_file = Sam_unmap(file)
        Sam_sorting(new_file, BAM= True)
        os.remove(input_Path+"/"+new_file)
        

    move_files(input_Path, savepath+'/Workfiles/', '.bam')
    
    files = file_list(savepath+'/Workfiles')

    os.chdir(savepath+'/Workfiles')

    if rm == True and rm_seq_dup == False:
        #In this case we will only move one folder before not as previously, as we want to keep only the sorted and maped
        #reads so we can just avoid the whole set of repeated files and the overfilling of memory.

        for file in files:
            command = "picard MarkDuplicates REMOVE_DUPLICATES=true I={} O={} M={}".format(file, "../marked_dup_"+ file[:-14] +".bam",
                                                                                           "../marked_dup_metrics_"+ file[:-14] +".bam")
            print(command)

            runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    elif rm_seq_dup == True and rm == False:
        for file in files:
      
            command = "picard MarkDuplicates REMOVE_SEQUENCING_DUPLICATES=true I={} O={} M={}".format(file, "../marked_dup_"+ file[:-14] +".bam",
                                                                                                "../marked_dup_metrics_"+ file[:-14] +".bam" )
            print(command)

            runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    elif rm_seq_dup == True and rm == True:   
        for file in files:
    
            command = "picard MarkDuplicates REMOVE_SEQUENCING_DUPLICATES=true REMOVE_DUPLICATES=true I={} O={} M={}".format(file, "../marked_dup_"+ file[:-14] +".bam",
                                                                                                "../marked_dup_metrics_"+ file[:-14] +".bam" )
            print(command)

            runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    else:
        for file in files:

            command = "picard MarkDuplicates I={} O={} M={}".format(file, "../marked_dup_"+ file[:-14] +".bam",
                                                                            "../marked_dup_metrics_"+ file[:-14] +".bam" )
            print(command)

            runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    return savepath


def run_Htseq(input_Path, savepath, annotated):

    files = file_list(input_Path)[1:]
    

    os.chdir(input_Path)

    for file in files:
 
            #Dont accept CRAM files,
            #Generates the complete set of files in order to indexing
            #Then to work the overalll  we will be focusing on the generation of the indexed woy.
            Sam_index(file)

            #You must be carefull when using the tool, you will have to set it up the most accurate to the real tools
            command = "htseq-count -f {} -r pos -t exon {} {} > {}".format( file[-3:],  file, annotated , savepath + '/' + file[:-4] + '_counts.txt')
            print(command)
            runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)


    #Inside of this will become moving the files to propper places in order to not miss the whole errors
    #So we will be able inside the knowledge and start working with the different files.            
    move_files(input_Path, savepath, ".txt")
    move_files(input_Path, savepath + "/Workfiles", ".bai")

    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    return savepath


#============================================================================================================================#
#                                                   MAPPING TOOLS                                                            #
#============================================================================================================================#


def run_STAR(input_Path, savepath, reference, gtf_file):

    '''
    Must be selected a Reference genome and it's propper ammount of information so we will
    also be carefull with the ammount of data provided, the GTF must be selected in concience

    Reference must be a directory where it is located the index of the new genome,
    then the Fasta file it's the actual reference genome, then GTF file explained before and the next value
    it's left by default.
    '''

    start = time.time()
    os.chdir(input_Path)
    
    input_Path_archives = Files_Mating(input_Path)
    print(input_Path_archives)

    command = "STAR --runThreadN 16 --runMode genomeGenerate --genomeDir {} --genomeFastaFiles {} --sjdbGTFfile {} --sjdbOverhang 100".format(savepath + "/Workfiles", reference, gtf_file )
    print(command)
    runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)


    for set_files in input_Path_archives:

        command = "STAR --genomeDir {} --runThreadN 16 -- readFilesIn {} {} --outFileNamePrefix {}".format(savepath + "/Workfiles", set_files, set_files[:-2]+"2P", set_files[:-2])
        print(command)

        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    move_files(savepath + "/Workfiles", input_Path, ".sam")

    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    end = time.time()
    print("Time elapsed on the whole execution = {} sec".format(end - start))
    return savepath


def run_Bowtie(input_Path, savepath, reference):

    start = time.time()

    os.chdir(input_Path)

    input_Path_archives = Files_Mating(input_Path)
    print(input_Path_archives)

    command = "bowtie2-build {} {}".format(reference, reference[61:-4] )
    print(command)
    runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    move_files(input_Path, reference[:61], ".bt2")

    os.chdir(reference[:61])
    for set_files in input_Path_archives:

        command = "bowtie2 -x {} -q -1 {} -2 {} -S {}".format(reference[61:-4], input_Path + "/" + set_files, input_Path + "/" + set_files[:-2]+"2P",  set_files[:-2]+".sam")
        print(command)

        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    move_files(savepath + "/Workfiles", savepath, ".sam")

    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    end = time.time()
    print("Time elapsed on the whole execution = {} sec".format(end - start))
    return savepath


###MUST CHECK THE BEHAVIOUR OF TOPHAT FUNCTION IN ORDER TO GET GET BEST ACCURACY SELECTING THE OVERALL:
def run_TopHat2(input_Path, savepath, reference, gtf_file, NO_WORKING = True):
    '''
    Needs reference and annotations
    Uses the indexing of bowtie2, needed in order to achieve faster the differenct files
    We have to look after how to install in a correct way inside mamba enviroment.
    '''

    start  = time.time()
    os.chdir(input_Path)

    input_Path_archives = Files_Mating(input_Path)
    print(input_Path_archives)

    command = "bowtie2-build {} {}".format(reference, savepath + "/Workfiles/" + reference[:-4] + ".btindex" )
    print(command)
    runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)


    for set_files in input_Path_archives:

        command = "tophat -p 16 -G {} -o {} --no-novel-juncs {} {} {}".format(gtf_file, set_files[:-2]+".sam", savepath + "/Workfiles/" + reference[61:-4] + ".btindex", set_files, set_files[:-2]+"2P")
        print(command)

        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    move_files(savepath + "/Workfiles", savepath, ".sam")

    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    end = time.time()
    print("Time elapsed on the whole execution = {} sec".format(end - start))
    return savepath


def run_BWA(reference, input_Path, savepath):

    start = time.time()

    os.chdir(input_Path)

    input_Path_archives = Files_Mating(input_Path)
    print(input_Path_archives)

    command = "bwa index {}".format(reference)
    print(command)
    runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    for set_files in input_Path_archives:

        command = "bwa mem {} {} {}  > {}".format(reference,  set_files, set_files[:-2]+"2P", set_files[:-2] + ".sam")
        print(command)

        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    move_files(input_Path, savepath, ".sam")

    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    end = time.time()
    print("Time elapsed on the whole execution = {} sec".format(end - start))
    return savepath


def run_SOAP2(reference, input_Path, savepath):

    start = time.time()

    os.chdir(input_Path)

    input_Path_archives = Files_Mating(input_Path)
    print(input_Path_archives)

    command = "2bwt-builder {}".format(reference)
    print(command)
    runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    for set_files in input_Path_archives:

        command = "soap -D {} -a {} -b {} -o {}".format(reference + ".index", set_files, set_files[:-2]+"2P", set_files[:-2] + ".soap" )
        print(command)

        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

        command = "soap2sam.pl {} > {}".format(set_files[:-2] + ".soap", set_files[:-2] + ".soap.sam")
        print(command)

        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

        command = "samtools -T {} {} > {}".format(reference, set_files[:-2] + ".soap.sam", set_files[:-2] + ".sam")
        print(command)

        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    move_files(input_Path, savepath, ".sam")

    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    end = time.time()
    print("Time elapsed on the whole execution = {} sec".format(end - start))
    return savepath


def run_Kallisto(reference, input_Path, savepath, gtf_file):

    start = time.time()

    os.chdir(input_Path)

    input_Path_archives = Files_Mating(input_Path)
    print(input_Path_archives)

    command = "kallisto index -i {} {}".format(reference[61:-4] + ".idx", reference)
    print(command)
    runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    for set_files in input_Path_archives:
        
        ###REMIND THAT THE DATA MUST BE NEEDED INTO FASTQ

        command = "kallisto quant --pseudobam -i {} -o {} --gtf {} -b 100 {} {}".format(reference[-39:-4] + ".idx",set_files[:-3]  ,gtf_file ,set_files, set_files[:-2]+"2P" )
        print(command)

        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    move_files(input_Path, savepath, ".sam")

    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    end = time.time()
    print("Time elapsed on the whole execution = {} sec".format(end - start))
    return savepath


def run_NovoAlign(reference, input_Path, savepath):

    start = time.time()

    os.chdir(input_Path)

    input_Path_archives = Files_Mating(input_Path)
    print(input_Path_archives)

    command = "novoindex {} {}".format(reference[61:-4] + ".nix", reference)
    print(command)
    runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    for set_files in input_Path_archives:
        
        command = "novoalign -d {} -o {} -f {} {}".format(reference[:-4] + ".nix",set_files[:-3]  + ".sam", set_files, set_files[:-2]+"2P" )
        print(command)

        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    move_files(input_Path, savepath, ".sam")

    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    end = time.time()
    print("Time elapsed on the whole execution = {} sec".format(end - start))
    return savepath


###'Have his own enviroment due problems with different versions so we will avoid
### With a new enviroment. Problems we can see so far, as the indexing it'S done with
### Bowtie2 algorithm will probably generate similar nor equal counts further on.

def run_MapSplicer(reference, input_Path, savepath):

    start = time.time()

    os.chdir(input_Path)

    input_Path_archives = Files_Mating(input_Path)
    print(input_Path_archives)

    command = "bowtie2-build {} {}".format(reference, savepath + "/Workfiles/" + reference[:-4] + ".btindex" )
    print(command)
    runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    for set_files in input_Path_archives:
        
        command = "python mapsplice.py -p 16 -x {}  -1 {} -2 {} -o {}".format(reference, savepath + "/Workfiles/" + reference[:-4] + ".btindex", set_files, set_files[:-2]+"2P" , savepath)
        print(command)

        runs = sp.run(command, shell=True, stderr=sp.PIPE, stdout=sp.PIPE, text=True)

    move_files(input_Path, savepath, ".sam")

    os.chdir("/mnt/Viro_Data/Mitarbeiter/Ian")
    end = time.time()
    print("Time elapsed on the whole execution = {} sec".format(end - start))
    return savepath

#============================================================================================================================#
#                                                OBTAINING  RESULTS                                                          #
#============================================================================================================================#


### Example of use, in order to use this you can just import the functions inside the other scripts or directly
### apply the usage . Another options could be ust changing the paths of the bottom functions, this will create the files and 
### store them in the different folders, notice that some of them will be requiring to add a Workfile generation of the propper files.

### Download the full set of data on the 6 different states for the control and the afected samples 
### See its initial behaviour on the overall sense. 
#DMP = Fastqdump(input_path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/FasqDumpRes/Workfiles", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/FasqDumpRes", old=False, paired=True)
#FAS = run_Fastqc(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/FasqDumpRes", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/Fastqc1Res")
#MUL = run_Multiqc(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/Fastqc1Res", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/MultiqcRes")


###Check the trimmed versions, separing the ones that are Paired, Unpaired, etc, then we will sort and work with them

#TRIM = run_Trimmomatic(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/FasqDumpRes", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes" )
#FAS2 =run_Fastqc(input_Path= TRIM, savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/Fastqc2Res")
#MUL2 = run_Multiqc(input_Path= FAS2, savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/MultiqcRes")


### Do the final analysis, remove the fies , work and map the whole set of files and then cut the ones we're no interested(MarkDuplicates). Finally
### Take and count the different elements

#BUILD = run_HiSat2_build(input_Path="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles", reference= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/GCF_000001405.40_GRCh38.p14_genomic.fna", simple= True )
#HS2 = run_HiSat2(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res",
#                  reference= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/GCF_000001405.40_GRCh38.p14_genomic" )
#PIC = run_Picard(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res" , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", rm= True, rm_seq_dup= False)
#HTS =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/genomic.gtf")

#============================================================================================================================#
#                                              TESTING MAPPING TOOLS                                                         #
#============================================================================================================================#

### We will start with running over STAR, so we can do the most direct comparation (consider that we're indexing )
### Just before inside the run_STAR function.

#STAR = run_STAR(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/STARRes",
#                  reference= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/GCF_000001405.40_GRCh38.p14_genomic.fna",  gtf_file= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/genomic.gtf")
#PIC = run_Picard(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/STARRes" , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", rm= True, rm_seq_dup= False)
#HTS =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/genomic.gtf")


### Over this set of different files we won't need at first the propper annotation file of the genomic sequence, but further on it will be used.

#Bow = run_Bowtie(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/BowtieRes",
#                  reference= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/GCF_000001405.40_GRCh38.p14_genomic.fna")
#PIC = run_Picard(input_Path= Bow , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", rm= True, rm_seq_dup= False)
#HTS =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/genomic.gtf")


#SOAP2 = run_SOAP2(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/TrimmomaticRes", savepath="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/SOAP2Res",
#                  reference= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/GCF_000001405.40_GRCh38.p14_genomic.fna")
                  
#PIC = run_Picard(input_Path= SOAP2 , savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", rm= True, rm_seq_dup= False)
#HTS =  run_Htseq(input_Path= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/PicardRes", savepath= "/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HTseqRes", annotated="/mnt/Viro_Data/Mitarbeiter/Ian/TFG+Leyla/HiSat2Res/Workfiles/genomic.gtf")
