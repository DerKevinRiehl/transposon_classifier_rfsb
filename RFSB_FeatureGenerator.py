############################################################################
##### Transposon Classifier RFSB - part of Transposon Ultimate #############
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

## Imports
from Bio.Blast.Applications import NcbirpstblastnCommandline
import os

## Methods
def countKmerOccurence(seqDNA, kmer):
    count = 0
    fromIndex = seqDNA.find(kmer)
    while fromIndex !=-1:
        fromIndex = seqDNA.find(kmer, fromIndex+1)
        count += 1
    return count

def createFeatures_kmer(fastaFile, kmerConfigFile, outputFile):
    # Load configuration
    f = open(kmerConfigFile, "r")
    config = f.readline().replace("\n","")
    f.close()
    head = config
    attributes = config.split(",")
    # Generate features and write output file
    f1 = open(outputFile, "w+")
    f1.write(head)
    f1.write("\n")
    f2 = open(fastaFile, "r")
    line = f2.readline()
    seqHead = ""
    seqDNA = ""
    while line!="":
        if line.startswith(">"):
            if(seqHead!=""):
                feature = ""
                for kmer in attributes:
                    if(len(kmer)==0 or len(seqDNA)==0):
                        relKmer = 0
                    else:
                        relKmer = countKmerOccurence(seqDNA, kmer)/(len(seqDNA)/len(kmer))-0.5
                    feature = feature + str( relKmer ) + ","
                f1.write(feature+"@"+seqHead.replace("@","")+"\n")
            seqHead = line.replace("\n","")
            seqDNA = ""
        elif line!="\n":
            seqDNA = seqDNA + line.replace("\n","")
        line = f2.readline()
    if(seqHead!=""):
        feature = ""
        for kmer in attributes:
            if(len(kmer)==0 or len(seqDNA)==0):
                relKmer = 0
            else:
                relKmer = countKmerOccurence(seqDNA, kmer)/(len(seqDNA)/len(kmer))-0.5
            feature = feature + str( relKmer ) + ","
        f1.write(feature+"@"+seqHead.replace("@","")+"\n")
    f1.close()
    f2.close()
        
def createFeatures_prot(fastaFile, db, outputFile1, outputFile2):
    # Run RPSTBLASTN
    #cline = NcbirpstblastnCommandline(cmd='rpstblastn', query=fastaFile, db=db[:-3], outfmt="7 sacc stitle qframe evalue bitscore qstart qend qlen sstart send", out=outputFile1)
    #cline()
    print("rpstblastn -query "+fastaFile+" -db "+db[:-3]+" -outfmt \"7 sacc stitle qframe evalue bitscore qstart qend qlen sstart send\" -evalue 0.1 -out "+outputFile1)
    os.system("rpstblastn -query "+fastaFile+" -db "+db[:-3]+" -outfmt \"7 sacc stitle qframe evalue bitscore qstart qend qlen sstart send\" -evalue 0.1 -out "+outputFile1)
    # Load Protein list
    proteins = list()
    f = open(db,"r")
    line = f.readline()
    while line !="":
        proteins.append(line.split(".")[0].lower())
        line = f.readline()
    f.close()
    # Count Protein Occurence
    f1 = open(outputFile2, "w+")
    f2 = open(outputFile1, "r")
    line = f2.readline()
    temp = list()
    lastTitle = ""
    lastLine = ""
    while line!="":
        if(line.startswith("#")):
            if len(temp)!=0:              
                f1.write(lastTitle.replace("@","")+"@")
                for i in range(0,len(temp)):
                    f1.write(temp[i]+"@")
                f1.write("\n")
            temp = list()
            for i in range(0,len(proteins)):
                temp.append("")
            if(line.startswith("# BLAST ")):
                if(lastLine.startswith("# 0 hits")):
                    f1.write(lastTitle.replace("@","")+"@")
                    for i in range(0,len(temp)):
                        f1.write(temp[i]+"@")
                    f1.write("\n")
                break
            else:
                lastTitle = f2.readline()[9:].replace("\n","")
                lastLine = line
                line = f2.readline()
                while(line!="" and line.startswith("#") and not line.startswith("# RPSTBLASTN") and not line.startswith("# BLAST")):
                    lastLine = line
                    line = f2.readline()
                if(line.startswith("# BLAST ")):
                    if(lastLine.startswith("# 0 hits")):
                        f1.write(lastTitle.replace("@","")+"@")
                        for i in range(0,len(temp)):
                            f1.write(temp[i]+"@")
                        f1.write("\n")
                    break
        else:
            parts = line.replace("\n","").split("\t")
            for i in range(0,len(parts)):
                parts[i] = parts[i].lower()
            key = parts[1].split(",")[0]
            eVal = parts[3]
            for i in range(0,len(temp)):
                if key==proteins[i]:
                    temp[i] = temp[i] + "," + eVal
            lastLine = line
            line = f2.readline()
    f1.close()
    f2.close()