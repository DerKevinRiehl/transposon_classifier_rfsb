############################################################################
##### Transposon Classifier RFSB - part of Transposon Ultimate #############
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# classes = "1,1/1,1/1/1,1/1/2,1/1/3,1/2,1/2/1,1/2/2,2,2/1,2/1/1,2/1/2,2/1/3,2/1/4,2/1/5,2/1/6,2/2,2/3".split(",")

# classes_sort = "1,2,1/1,1/2,2/1,2/2,2/3,1/1/1,1/1/2,1/1/3,1/2/1,1/2/2,2/1/1,2/1/2,2/1/3,2/1/4,2/1/5,2/1/6".split(",")
# classes_parent = {"1":"no", "2":"no", "1/1":"1", "1/2":"1", "2/1":"2", "2/2":"2", "2/3":"2", "1/1/1":"1/1", "1/1/2":"1/1", "1/1/3":"1/1",
#                   "1/2/1":"1/2", "1/2/2":"1/2", "2/1/1":"2/1", "2/1/2":"2/1", "2/1/3":"2/1", "2/1/4":"2/1", "2/1/5":"2/1", "2/1/6":"2/1"}
# classes_neighb = {"1":["2"], "2":["1"], "1/1":["1/2"], "1/2":["1/1"], "2/1":["2/2", "2/3"], "2/2":["2/1", "2/3"], "2/3":["2/1", "2/2"]}
# names = {"1" :  "Class I, Retrotransposon" ,
# "1/1" :  "LTR, Retrotransposon" ,
# "1/1/1" :  "Copia, LTR, Retrotransposon" ,
# "1/1/2" :  "Gypsy, LTR, Retrotransposon" ,
# "1/1/3" :  "ERV, LTR, Retrotransposon" ,
# "1/2" :  "Non-LTR, Retrotransposon" ,
# "1/2/1" :  "LINE, Non-LTR, Retrotransposon" ,
# "1/2/2" :  "SINE, Non-LTR, Retrotransposon" ,
# "2" :  "DNA Transposon" ,
# "2/1" :  "TIR, DNA Transposon" ,
# "2/1/1" :  "Tc1-Mariner, TIR, DNA Transposon" ,
# "2/1/2" :  "hAT, TIR, DNA Transposon" ,
# "2/1/3" :  "CMC, TIR, DNA Transposon" ,
# "2/1/4" :  "Sola, TIR, DNA Transposon" ,
# "2/1/5" :  "Zator, TIR, DNA Transposon" ,
# "2/1/6" :  "Novosib, TIR, DNA Transposon" ,
# "2/2" :  "Helitron, DNA Transposon" ,
# "2/3" :  "MITE, DNA Transposon"}

## Imports
import sys

## Methods
def processClasses(classes):
    classes_parent = {}
    for c in classes:
        if(not "/" in c):
            classes_parent[c] = "no"
        else:
            assParent = ""
            parts = c.split("/")
            for i in range(0,len(parts)-1):
                assParent = assParent + parts[i] +"/"
            classes_parent[c] = assParent[:-1]
            
    levels = list()
    for c in classes:
        levels.append(len(c.split("/")))
    classes_sorted = [classes for _,classes in sorted(zip(levels,classes))]
    return classes_parent, classes_sorted, levels

def loadTaxonomy_from_config(configFile):
    f = open(configFile)
    lines=f.readlines()
    f.close()
    for i in range(0,len(lines)):
        lines[i] = lines[i].replace("\n","")
    
    classes = list()
    for line in lines:
        classes.append(line.split(":")[0].replace(".","/"))
        
    names = {}
    for line in lines:
        names[line.split(":")[0]] = line.split(":")[1]
  
    classes_parent, classes_sorted, levels = processClasses(classes)
    
    return classes, names, classes_parent, classes_sorted, levels

def loadTaxonomy(labelFile, configFile=""):
    # Load taxonomy from label file
    classesL = list()
    f = open(labelFile,"r")
    line = f.readline()
    while line!="":
        if(line.startswith(">")):
            line = f.readline().replace("\n","").replace(".","/")
            if not line in classesL:
                if(len(classesL)==0):
                    classesL.append(line)
                else:
                    inserted = False
                    for i in range(0,len(classesL)):
                        if(line < classesL[i]):
                            classesL.insert(i, line)
                            inserted = True
                            break
                    if not inserted:
                        classesL.append(line)
        line = f.readline()
    f.close()
    # If config exists,try to analyse it
    if(configFile!=""):
        classesC, namesC, classes_parentC, classes_sortedC, levelsC = loadTaxonomy_from_config(configFile)
        # double check if any class in labelFile does not occur in config file
        inconL = list()
        for c in classesL:
            if not c in classesC:
                inconL.append(c)
        if(len(inconL)!=0):
            print("ERROR: Labelfile and Configuration file have incongruent taxonomies")
            print("Labelfile:  ", labelFile)
            print("Configfile: ", configFile)
            print("List of incongruent classes:")
            for i in inconL:
                print(">>  \"",i,"\"")
            print("PROGRAMME SHUTDOWN")
            sys.exit(-1)
        # if everything okay, return config file based taxonomy
        return classesC, namesC, classes_parentC, classes_sortedC, levelsC
    # If no configuration available
    else:
        while True:
            namesL = {}
            for c in classesL:
                namesL[c] = "-"
            classes_parentL, classes_sortedL, levelsL = processClasses(classesL)  
            found = False
            for p in classesL:
                ex = classes_parentL[p]
                if ex!="no" and ex not in classesL:
                    inserted = False
                    for i in range(0,len(classesL)):
                        if(ex < classesL[i]):
                            classesL.insert(i, ex)
                            inserted = True
                            break
                    if not inserted:
                        classesL.append(ex)
                    found = True       
            if not found:
                break
        return classesL, namesL, classes_parentL, classes_sortedL, levelsL