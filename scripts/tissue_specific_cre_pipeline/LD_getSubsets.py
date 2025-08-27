##########################################################################################################################################
# USAGE: python LD_getSubsets.py file1 file2 file3 etc... outputDir
# Requirements: Latest version of BEDtools
# 18/11/2016
##########################################################################################################################################
import sys
import os
import string
import random
import re


def main():
    #inputDir = sys.argv[1]
    outputDir = sys.argv[len(sys.argv)-1]
    
    #--help option - HN 20210813
    if outputDir.endswith("--help") or outputDir.endswith("-h"):
        print('Usage: python3 LD_getSubsets.py [INPUT_FILES] [OPTIONS]... >> [OUTPUT_DIR}\n  Input files: BED file format, each separated by a space\n  Options:\n    --help: Show this message and exit.\n  Output directory:  A directory path where the output BED files are deposited')
        sys.exit()
        
    elif not outputDir.endswith("/"):
        sys.exit("Error: must specify output directory")
    
    files = []
    print(f"outputDir: {outputDir}")

    print(f"arguments: {sys.argv}")
    for f in sys.argv :
        if f.endswith(".bed") or f.endswith("*") :
            files.append(f)

    files = ' '.join(files)
    print(f"files are: {files}")

    #create temporary directory within output for temporary files
    tempDir = outputDir + "temp/"
    if not os.path.exists(tempDir):
        os.makedirs(tempDir)

    
    print(f"Start concatenating bed files")
    # #concatenate bed files into one single bed file
    # #os.system('cat %s | sort -n -k 2 > %s' % (inputDir + "*", tempDir + "all.txt"))
    os.system('cat %s | sort -k 1,1 -k 2,2n > %s' % (files, tempDir + "all.txt"))
        
    print(f"Reach end of file")


    # #run BEDTools' mergeBed on the concatenated file all.txt
    # #os.system('/home/tradoa/apps/BEDTools-Version-2.11.2/bin/mergeBed -nms -d 100 -i %s > %s' % (tempDir + "all.txt", tempDir + "all_merged.txt"))
    # #os.system('/Users/ramialis/Documents/MIMI_LAB/001_PROJECTS/000_PIPELINES/000_packages/bedtools2-2.19.1/bin/mergeBed -n -i %s > %s' % (tempDir + "all.txt", tempDir + "all_merged.txt")) #####nmr deprecated
    # #os.system('/Users/ramialis/Documents/MIMI_LAB/001_PROJECTS/000_PIPELINES/000_packages/bedtools2-2.19.1/bin/mergeBed -n -c4 -ocollapse -delim "_" -i %s > %s' % (tempDir + "all.txt", tempDir + "all_merged.txt"))
    os.system('mergeBed -c 4 -o collapse -delim ";" -i %s > %s' % (tempDir + "all.txt", tempDir + "all_merged.txt"))

    # #Split the merged file into subsets
    annotate_result(tempDir + "all_merged.txt", outputDir)
    

import os
import re

import re

def clean_id_string(raw):
    parts = raw.strip().split(";")
    cleaned = []
    for part in parts:
        part = re.sub(r'_\d+$', '', part)
        if part.startswith('ovary_'):
            cleaned.append('ovary')
        elif part.startswith('testes_'):
            cleaned.append('testes')
        elif part.startswith('stomach_'):
            cleaned.append('stomach')
        elif part.startswith('eye_'):
            cleaned.append('eye')
        elif part.startswith('heart_'):
            cleaned.append('heart')
        elif part.startswith('cerebrum_'):
            cleaned.append('cerebrum')
        else:
            cleaned.append(part)
    # remove duplicates, preserve order
    seen = set()
    final = []
    for x in cleaned:
        if x not in seen:
            seen.add(x)
            final.append(x)
    return ";".join(final)


def annotate_result(bedFile, outputDir):
    bedFile = open(bedFile, "r")
    print(f"bedFile: {bedFile}")
    
    for ind, line  in enumerate(bedFile):
        print(f"Iteration {ind}")
        line = line.strip()
        linel = line.split("\t")

        id = linel[3]
        print(f"id: {id}")


        # perform cleaning on id

        
        id_2 = clean_id_string(id)
            
        print(f"id2: {id_2}")
        
        # id_2 = id_2.replace('|', ';')
        # id_2 = id_2.replace('Dam', '')

        factorList = id_2.split(";")
        for i, name in enumerate(factorList):
            factorList[i] = name.strip()
        print(f"factorList: {factorList}")

        # serve the purpose of bed file name
        # Separate the 'ovary' containing elements and others
        ovary_list = sorted(set([item for item in factorList if 'ovary' in item]))
        non_ovary_list = sorted(set([item for item in factorList if 'ovary' not in item]))


        factorList = ovary_list + non_ovary_list

        print(f"factorList: {factorList}")

        id_2 = "_".join(factorList)

        newLine = "\t".join(linel)
        outputFile = outputDir + id_2 + ".bed"

        print(f"outputFile: {outputFile}")
        print(f"id_2: {id_2}")
        os.system('echo "%s" >> %s' % (newLine, outputFile))	

        # break

        # if ind > 200: 
        #     break

        
    bedFile.close()

    
main()
