#!/usr/bin/env python
############################################################################
##### Transposon Classifier RFSB - part of Transposon Ultimate #############
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# @author Kevin Riehl 2020/21 <kevin.riehl.de@gmail.com>

# EXAMPLE 1: Predict using the model
# fastaFile           = "Transposon_Sequences_PRJEB28388.fa"
# modelFile           = "ergebnisTestModel.model"#"models/TE_Classifier_RFSB_models_All_Big.pickle"
# kmerConfigFile      = "config/kmers.txt"
# dbFile              = "config/RPSTBLASTN_LIB/db_large.pn"
# outputFile          = "ergebnisTestPythonPipeline3.txt"
# tempFileA = "testA3.txt"#"testA.txt"
# tempFileB = "testB3.txt"#"testB.txt"
# tempFileC = "testC3.txt"#"testC.txt"
# createFeatures_kmer(fastaFile, kmerConfigFile, tempFileA)
# createFeatures_prot(fastaFile, dbFile, tempFileB, tempFileC)
# classifyTransposons(tempFileA, tempFileC, dbFile, modelFile, fastaFile, outputFile)

# EXAMPLE 2: Train a model
# fastaFile           = "G:/CambridgeGenData/TE_DB/ALL_tiny.fasta"
## kmerConfigFile      = "config/kmers.txt"
# taxonomyConfigFile  = "config/taxonomy.txt"
# dbFile              = "config/RPSTBLASTN_LIB/db_large.pn"
# outputFile          = "ergebnisTestModel.model"
# eThreshold          = 5.0
# labelFile           = "G:/CambridgeGenData/TE_DB/ALL_tiny.label"
# tempFileA           = "testA2.txt"
# tempFileB           = "testB2.txt"
# tempFileC           = "testC2.txt"
# createFeatures_kmer(fastaFile, kmerConfigFile, tempFileA)
# createFeatures_prot(fastaFile, dbFile, tempFileB, tempFileC)
# classes, names, classes_parent, classes_sort, levels = loadTaxonomy(labelFile, configFile="")
# trainClassificationModel(classes, names, classes_sort, classes_parent, levels, tempFileA, tempFileC, labelFile, dbFile, eThreshold, outputFile)

## Imports
from RFSB_FeatureGenerator import createFeatures_kmer, createFeatures_prot
from RFSB_ModelManager import classifyTransposons, trainClassificationModel
from RFSB_Evaluator import algorithmEvaluationSummary
from RFSB_TaxonomyLoader import loadTaxonomy
from RFSB_Help import helpExplanations
from datetime import datetime
import os.path
import sys
import os

## Methods
def getArgument(args, title):
    for i in range(0, len(args)):
        if(args[i].startswith("-"+title)):
            if(i<len(args)-1):
                return args[i+1]
            else:
                return ""
    return ""

def checkEmptyFile(fastaFile):
    f = open(fastaFile, "r")
    content = f.readline()
    f.close()
    return content == ""

## MAIN CODE Examples
# Parameters
args = sys.argv[1:]
mode = getArgument(args, "mode")
script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in    

if(mode=="trainModel"):
    # get mandatory params
    fastaFile           = getArgument(args, "fastaFile")
    labelFile           = getArgument(args, "labelFile")
    outputFile          = getArgument(args, "outputModelFile")
    if(not os.path.isfile(fastaFile)):
        print("ERROR: Fasta File \"",fastaFile,"\" could not be found...EXIT")
    if(not os.path.isfile(labelFile)):
        print("ERROR: Label File \"",labelFile,"\" could not be found...EXIT")
    if(outputFile==""):
        print("ERROR: Please specify \"outputModelFile\"...EXIT")
    # get optional params
    taxonomyConfigFile  = getArgument(args, "taxonomyConfigFile")
    dbFile              = getArgument(args, "proteinDBFile")
    kmerConfigFile      = getArgument(args, "kmerConfigFile")
    eThreshold          = getArgument(args, "eThreshold")
    tempFolder          = getArgument(args, "tempFolder")
    delTemp             = getArgument(args, "deleteTempFiles")
    if(dbFile==""):
        dbFile = os.path.join(script_dir,"config","RPSTBLASTN_LIB","db_large.pn")
#        dbFile = "config/RPSTBLASTN_LIB/db_large.pn"
        if(not os.path.isfile(dbFile)):
            print("ERROR: ProteinDB \".pn\" File was not handed over. Standard file was used and could not be found. Please navigate into the folder of \"RFSB_TE_Classifier.py\" and run command again...EXIT")
    elif(not os.path.isfile(dbFile)):
        print("ERROR: ProteinDB \".pn\" File \"",dbFile,"\" could not be found...EXIT")
    if(kmerConfigFile==""):
        kmerConfigFile = os.path.join(script_dir,"config","kmers.txt")
#        kmerConfigFile = "config/kmers.txt"
        if(not os.path.isfile(kmerConfigFile)):
            print("ERROR: K-mer configuration file was not handed over. Standard file was used and could not be found. Please navigate into the folder of \"RFSB_TE_Classifier.py\" and run command again...EXIT")
    elif(not os.path.isfile(kmerConfigFile)):
        print("ERROR: K-mer configuration file \"",dbFile,"\" could not be found...EXIT")
    if(not eThreshold==""):
        eThreshold = float(eThreshold)
    else:
        print("[INFO] Parameter \"eThreshold\" not set. Set to default value: '5.0'.")
        eThreshold = 5.0
    if(tempFolder==""):
        if("/" in outputFile):
            parts = outputFile.split("/")[:-1]
            tempFolder = "/".join(parts)+"/"
        else:
            tempFolder = ""
    if(delTemp==""):
        print("[INFO] Parameter \"deleteTempFiles\" not set. Set to default value: 'True'.")
        delTemp = True
    elif(delTemp=="False" or delTemp=="FALSE"):
        delTemp = False
    else:
        delTemp = True
    # check if fasta file is empty
    if(checkEmptyFile(fastaFile)):
        print("ERORR: mode \"",mode,"\" FASTA File is emtpy...EXIT")
        sys.exit(0)
    # temporary files
    tempFileA           = tempFolder+outputFile.split("/")[-1]+".featuresA_"+str(datetime.now()).replace("-","").replace(":","").replace(" ","_").split(".")[0]+".temp"
    tempFileB           = tempFolder+outputFile.split("/")[-1]+".featuresB_"+str(datetime.now()).replace("-","").replace(":","").replace(" ","_").split(".")[0]+".temp"
    tempFileC           = tempFolder+outputFile.split("/")[-1]+".featuresC_"+str(datetime.now()).replace("-","").replace(":","").replace(" ","_").split(".")[0]+".temp"   
    # train model
    print("[INFO] Generate kMer features...")
    createFeatures_kmer(fastaFile, kmerConfigFile, tempFileA)
    print("[INFO] Generate protein features using RPSTBLASTN...")
    createFeatures_prot(fastaFile, dbFile, tempFileB, tempFileC)
    print("[INFO] Load / Generate taxonomy...")
    classes, names, classes_parent, classes_sort, levels = loadTaxonomy(labelFile, configFile=taxonomyConfigFile)
    print("[INFO] Train model...")
    trainClassificationModel(classes, names, classes_sort, classes_parent, levels, tempFileA, tempFileC, labelFile, dbFile, eThreshold, outputFile)
    print("[INFO] Finished...")
    # delete temporary files
    if(delTemp):
        os.remove(tempFileA) 
        os.remove(tempFileB) 
        os.remove(tempFileC) 
elif(mode=="classify"):
    # get mandatory params
    fastaFile           = getArgument(args, "fastaFile")
    outputFile          = getArgument(args, "outputPredictionFile")
    if(not os.path.isfile(fastaFile)):
        print("ERROR: Fasta File \"",fastaFile,"\" could not be found...EXIT")
    if(outputFile==""):
        print("ERROR: Please specify \"outputPredictionFile\"...EXIT")
    # get optional params
    modelFile           = getArgument(args, "modelFile")
    kmerConfigFile      = getArgument(args, "kmerConfigFile")
    dbFile              = getArgument(args, "proteinDBFile")
    tempFolder          = getArgument(args, "tempFolder")
    delTemp             = getArgument(args, "deleteTempFiles")
    if(not os.path.isfile(modelFile)):
        modelFile = os.path.join(script_dir,"models","TE_Classifier_RFSB_models_All_Big.pickle")
#        modelFile = "models/TE_Classifier_RFSB_models_All_Big.pickle"
        if(not os.path.isfile(modelFile)):
            print("ERROR: 'modelFile' was not handed over. Standard file was used and could not be found. Please navigate into the folder of \"RFSB_TE_Classifier.py\" and run command again...EXIT")
    if(dbFile==""):
        dbFile = os.path.join(script_dir,"config","RPSTBLASTN_LIB","db_large.pn")
#        dbFile = "config/RPSTBLASTN_LIB/db_large.pn"
        if(not os.path.isfile(dbFile)):
            print("ERROR: ProteinDB \".pn\" File was not handed over. Standard file was used and could not be found. Please navigate into the folder of \"RFSB_TE_Classifier.py\" and run command again...EXIT")
    elif(not os.path.isfile(dbFile)):
        print("ERROR: ProteinDB \".pn\" File \"",dbFile,"\" could not be found...EXIT")
    if(kmerConfigFile==""):
        kmerConfigFile = os.path.join(script_dir,"config","kmers.txt")
#        kmerConfigFile = "config/kmers.txt"
        if(not os.path.isfile(kmerConfigFile)):
            print("ERROR: K-mer configuration file was not handed over. Standard file was used and could not be found. Please navigate into the folder of \"RFSB_TE_Classifier.py\" and run command again...EXIT")
    elif(not os.path.isfile(kmerConfigFile)):
        print("ERROR: K-mer configuration file \"",kmerConfigFile,"\" could not be found...EXIT")
    if(tempFolder==""):
        if("/" in outputFile):
            parts = outputFile.split("/")[:-1]
            tempFolder = "/".join(parts)+"/"
        else:
            tempFolder = ""
    if(delTemp==""):
        print("[INFO] Parameter \"deleteTempFiles\" not set. Set to default value: 'True'.")
        delTemp = True
    elif(delTemp=="False" or delTemp=="FALSE"):
        delTemp = False
    else:
        delTemp = True
    # temporary files
    tempFileA           = tempFolder+outputFile.split("/")[-1]+".featuresA_"+str(datetime.now()).replace("-","").replace(":","").replace(" ","_").split(".")[0]+".temp"
    tempFileB           = tempFolder+outputFile.split("/")[-1]+".featuresB_"+str(datetime.now()).replace("-","").replace(":","").replace(" ","_").split(".")[0]+".temp"
    tempFileC           = tempFolder+outputFile.split("/")[-1]+".featuresC_"+str(datetime.now()).replace("-","").replace(":","").replace(" ","_").split(".")[0]+".temp"
    # check if fasta file is empty
    if(checkEmptyFile(fastaFile)):
        print("ERORR: mode \"",mode,"\" FASTA File is emtpy...EXIT")
        sys.exit(0)
    # classify TEs
    print("[INFO] Generate kMer features...")
    createFeatures_kmer(fastaFile, kmerConfigFile, tempFileA)
    print("[INFO] Generate Protein features using RPSTBLASTN...")
    createFeatures_prot(fastaFile, dbFile, tempFileB, tempFileC)
    print("[INFO] Load Model and classify...")
    classifyTransposons(tempFileA, tempFileC, dbFile, modelFile, fastaFile, outputFile)
    print("[INFO] Finished...")
    # delete temporary files
    if(delTemp):
        os.remove(tempFileA) 
        os.remove(tempFileB) 
        os.remove(tempFileC) 
elif mode=="evaluate":
    # get mandatory params
    predLabelFile = getArgument(args, "predLabelFile")
    trueLabelFile = getArgument(args, "trueLabelFile")
    if(not os.path.isfile(predLabelFile)):
        print("ERROR: Prediction Label File \"",predLabelFile,"\" could not be found...EXIT")
    if(not os.path.isfile(trueLabelFile)):
        print("ERROR: True Label File \"",trueLabelFile,"\" could not be found...EXIT")
    # get optional params
    printResults  = getArgument(args, "printResults")
    saveCSV       = getArgument(args, "saveResultsCSV")
    savePickle    = getArgument(args, "saveResultsPickle")
    csvFile       = getArgument(args, "outputCSVFile")
    pickleFile    = getArgument(args, "outputPickleFile")
    generatePNG   = getArgument(args, "generateDiagramPNG")
    generateSVG   = getArgument(args, "generateDiagramSVG")
    diagTitle     = getArgument(args, "diagramTitle")
    svgFile       = getArgument(args, "outputSVGFile")
    pngFile       = getArgument(args, "outputPNGFile")
    if(printResults=="False" or printResults=="FALSE"):
        printResults = False
    else:
        printResults = True
    if(saveCSV=="False" or saveCSV=="FALSE"):
        saveCSV = False
    else:
        saveCSV = True
    if(savePickle=="True" or savePickle=="true"):
        savePickle = True
    else:
        savePickle = False
    if(pickleFile==""):
        pickleFile = ".".join(predLabelFile.split(".")[:-1])+"_results.pickle"
    if(csvFile==""):
        csvFile = ".".join(predLabelFile.split(".")[:-1])+"_results.csv"
    if(generatePNG=="" or generatePNG=="True" or generatePNG=="true"):
        generatePNG = True
    else:
        generatePNG = False
    if(generateSVG=="" or generateSVG=="True" or generateSVG=="true"):
        generateSVG = True
    else:
        generateSVG = False
    if(svgFile==""):
        svgFile = ".".join(predLabelFile.split(".")[:-1])+"_diagram.svg"
    if(pngFile==""):
        pngFile = ".".join(predLabelFile.split(".")[:-1])+"_diagram.png"

    # calculate performance
    algorithmEvaluationSummary(predLabelFile, trueLabelFile, printResults, saveCSV, savePickle, generatePNG, generateSVG, diagTitle, csvFile, pickleFile, svgFile, pngFile)
elif mode=="help" or mode=="h" or "h" in args or "help" in args or "-h" in args or "-help" in args or "--h" in args or "--help" in args:
    helpExplanations()
else:
    print("ERORR: mode \"",mode,"\" unknown...EXIT")


