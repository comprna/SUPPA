# -*- coding: utf-8 -*-
"""
Created on Thu May 22 11:24:33 2014

@author: Gael P Alamancos
@email: gael.perez[at]upf.edu
"""

import sys
from argparse import ArgumentParser, RawTextHelpFormatter, SUPPRESS

description = \
    "Description:\n\n" + \
    "This little script takes 1 or more fields from multiple files with\n" + \
    "a common format and at least a common field which can be used as\n" + \
    "a unique identifier.\n" + \
    "If a identifier do not appear in a file, the corresponding field will\n" + \
    "be set up to 0."
parser = ArgumentParser(description = description, formatter_class=RawTextHelpFormatter)
parser.add_argument("-i", "--input-files", nargs="+", required = True,
                    help = "spaced separated list of files to join.")
parser.add_argument("-k", "--key-field", required = True,
                    help = "common field among the input files.")
parser.add_argument("-f", "--fields", nargs="+", required = True,
                    help = "spaced separated list of fields to select.(starting in 1)")
parser.add_argument("-s", "--separator", default = "\t",
                    help = "field separator")
parser.add_argument("-o", "--output-file", required=True,
                    help = "name of the output file.")
parser.add_argument("--no-header", action="store_true",
                    help = "use it if the file has no header.")
args = parser.parse_args()
try:
    #Recovering the parsed options
    filesList = args.input_files
    keyField = int(args.key_field)
    selectedFields = []
    for i in args.fields:
        selectedFields.append(int(i))
    separator = args.separator
    outputFile = args.output_file
    if args.no_header:
        header = False
    else:
        header = True
    #To remove the key-field from the selectedFields if it is there
    for i in range(len(selectedFields)):
        if selectedFields[i] == keyField:
            del selectedFields[i]
     
    dictionary = {} #key = keyField; value=list of fields to join from different files.         
    if header:
        dictionary["header"] = {}
    #BUFFERING INPUT
    #Loop through all the files    
    for inputFile in filesList:
        seen_header = not header
        print("INFO: Reading file: %s" % inputFile)
        lineNumber = 1
        with open(inputFile, 'r') as f:
            for l in f:
                line = l.rstrip('\n').split(separator)
                if header and not seen_header and (l.startswith('# Transcript') or not(l.startswith('#'))):
                    seen_header = True
                    #Storing the header
                    #print "ASS %s" % line
                    #print selectedFields
                    dictionary["header"][inputFile] =  []
                    for field in selectedFields:
                        dictionary["header"][inputFile].append(line[(field-1)])
                    lineNumber += 1
                    continue
                if l.startswith('#'):
                    continue
                #Initialize to an empty list if the key does not exist
                dictionary.setdefault(line[(keyField-1)], {})
                #Load all the fields from the file
                dictionary[line[(keyField-1)]][inputFile] = []
                for field in selectedFields:
                    dictionary[line[(keyField-1)]][inputFile].append(
                        str(line[(field-1)]))
                lineNumber += 1
        print("INFO: File %s closed." % inputFile)
    #Generating new header
    outputHeader = []
    for inputFile in filesList:
        #aux2 = os.path.basename(inputFile)
        if header:
            aux1 = inputFile.split("/")
            aux1.pop()
            aux2 = aux1.pop()
            #Add file+fields to the header
            for x in dictionary["header"][inputFile]:
                outputHeader.append(aux2)
                #outputHeader.append(aux2 + "_" + str(x))
        else:
            #Add file + num to the header
            aux1 = inputFile.split("/")
            aux1.pop()
            aux2 = aux1.pop()
            for x in range(1, len(selectedFields) + 1):
                outputHeader.append(aux2)
                #outputHeader.append(aux2 + "_" + str(x))

    #WRITING OUTPUT
    print("INFO: Writing output to %s" % outputFile)
    f = open(outputFile, 'w')
    #Writing the header in the output file
    f.write("\t".join(outputHeader) + "\n")
    #Looping through all the unique IDs avoiding the "header" .
    for key, value in [(x,y) for x, y in dictionary.items() if x != "header"]:

        line = []
        line.append(key)    #Adding the common id
        for inputFile in filesList:
            #If the unique id exists in inputFile
            if value.get(inputFile): 
                line += value.get(inputFile)
            #If the unique id does not exist in inputFile
            else:
                #Setting up to 0 all the selectedFields of inputfile
                for x in range(len(selectedFields)):
                    line.append(0)
        #Writing the line in the output
        f.write("\t".join(line) + "\n")
    f.close()
    print("INFO: %s closed." % outputFile)
    print("INFO: Process ended succesfully.")
except:
    print("ERROR: %s" % err)
    sys.exit(1)
