"""
@authors: Juan L. Trincado
@email: juanluis.trincado@upf.edu

format_unique_fasta_RefSeq_annotation.py: Takes a Fasta annotation and returns the same file without
 repeated ids. Salmon needs an annotation without repeated ids
 The discarded ids will be outputed in another file.
"""

import time, sys

try:

    print("Starting execution: " + time.strftime('%H:%M:%S') + "\n")

    fasta_path = sys.argv[1]
    output_path = sys.argv[2]
    output_discarded_path = sys.argv[3]

    # fasta_path = "/projects_rg/SUPPA2/general_files/refseq_hg19.fa"
    # output_path = "/projects_rg/SUPPA2/general_files/refseq_hg19_unique.fa"
    # output_discarded_path = "/projects_rg/SUPPA2/general_files/refseq_hg19_discarded.txt"


    # Open the connection to the fasta file
    fasta_file = open(fasta_path, 'r')
    # Create a new file for the output
    output_file = open(output_path, 'w')
    output_discarded_file = open(output_discarded_path, 'w')

    # Iterate over the file, saving in a dictionary each transcript id
    # If the trasncript id is not repeated, we will included in the final fasta_file
    # We will remove all the non default chromosomes
    dict = {}
    accepted_chr = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12",
                    "chr13","chr14","chr15","chr16","chr16","chr17","chr18","chr19","chr20","chr21","chr22",
                    "chrX","chrY"]
    accepted_flag = False
    for line in fasta_file:
        if(line[0]==">"):
            tokens = line.rstrip().split(" ")
            # If it doesn't exist, and the chromosome is from 1-22, X or Y, save it and change the flag
            # for saving the following sequence lines until the next transcript id
            id = tokens[0][14:]
            chr = tokens[1].split(":")[0].split("=")[1]
            if((chr in accepted_chr) and (id not in dict)):
                accepted_flag = True
                dict[id] = 0
                #Reformat the line and print it in the output file
                new_line = ">"+id+" "+" ".join(tokens[1:])
                output_file.write(new_line)
            else:
                accepted_flag = False
                #Save this repeated or discarded ids in antoher file
                output_discarded_file.write(line)
        else:
            #If accepted_flag is True, print lines in the output file
            if(accepted_flag):
                output_file.write(line)

    # Close the handlers
    output_file.close()

    print("Done. Exiting program. " + time.strftime('%H:%M:%S') + "\n")


except Exception as error:
    print('\nERROR: ' + repr(error))
    print("Aborting execution")
    sys.exit(1)
