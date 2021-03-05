############################################################################
##### Transposon Classifier RFSB - part of Transposon Ultimate #############
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.naive_bayes import GaussianNB

import warnings
import pickle
warnings.filterwarnings("ignore")

# Methods
def loadBinaryLabels(labelFile, targetLabel):
    f = open(labelFile, "r")
    line = f.readline()
    lL = list()
    while line!="":
        if(line.startswith(">")):
            lL.append(f.readline().replace("\n","").replace(".","/"))
        line = f.readline()
    f.close()
    labels = list()
    for i in range(0,len(lL)):
        if(str.startswith(lL[i], targetLabel)):
            labels.append(1)
        else:
            labels.append(0)
    labels = np.asarray(labels)
    return labels

def loadProteinFeatures(fileFeatures1, nfeatures, eThreshold):
    f = open(fileFeatures1, "r")
    line = f.readline()
    lF = list()
    while(line!=""):
        parts = line.replace("\n","").split("@")
        features = parts[-nfeatures:]
        lF.append(features)
        line = f.readline()
    f.close()
    featuresA = list()
    for i in range(0,len(lF)):
        r = list()
        for x in range(0,len(lF[i])):
            if(lF[i][x].replace("\n","")!=""):
                parts = lF[i][x].replace("\n","").split(",")
                parts = [z for z in parts if z]
                parts = [float(z) for z in parts]
                if(min(parts)<eThreshold):
                    r.append(1)
                else:
                    r.append(0)
            else:
                r.append(0)
        featuresA.append(r)
    featuresA = np.asarray(featuresA)
    return featuresA

def loadKmerFeatures(fileFeatures2):
    features2 = np.loadtxt(fileFeatures2, delimiter=',', skiprows=1, usecols=range(0,336))
    return features2

def loadFeatures(filePDB, fileFeatures1, fileFeatures2, eThreshold):
    NProteins = len(open(filePDB, "r").readlines())
    featuresProt = loadProteinFeatures(fileFeatures1, NProteins, eThreshold)
    featuresKmer = loadKmerFeatures(fileFeatures2)
    features = np.concatenate((featuresProt, featuresKmer), axis=1)
    return features

def generate(modelS):
    if(modelS=="RandomForest"):
        return RandomForestClassifier()
    elif(modelS=="LogisticRegression"):
        return LogisticRegression()
    elif(modelS=="GaussNB"):
        return GaussianNB()
    elif(modelS=="AdaBoost"):
        return AdaBoostClassifier()
    else:
        print("ERROR! UNKOWN MODEL!")

def trainModels_selStrat_binStruc(classes_sort, classes_parent, class_neighbors, labels_train, features_train, modelS, classes):
    # Generate Models
    models = {}
    for c in classes:
        models[c] = generate(modelS)

    # Train Models
    rulesD = {}
    for c in classes_sort:
        print("Training for class : "+c+"...")
        # Prepare Data from superior node
        parent = classes_parent[c]
        if(parent=="no"): # if highest node let pass all
            dataX = features_train
            print(features_train.shape)
            dataY = labels_train[c]
            rulesD[c] = [1 for x in range(0,len(labels_train[c]))]
        else: # if not, only use data passed by superior node
            res = list()
            res.append(models[parent].predict_proba(features_train)[:,1])
            parentIdx = classes.index(parent)
            for nei in class_neighbors[parentIdx]:
                predProba = models[nei].predict_proba(features_train)
                if(predProba.shape[1]==1):
                    res.append(predProba[:,0])
                else:
                    res.append(predProba[:,1])
            res = np.asarray(res)
            resArr = (res==np.max(res, axis=0))[0]
            ruleP = rulesD[parent]
            rule = [1 if (resArr[x]==True and ruleP[x]==1) else 0 for x in range(0,len(resArr)) ]
            rulesD[c] = rule
            
            dataX = [features_train[x] for x in range(0,len(rule)) if rule[x]==1]
            dataY = [labels_train[c][x] for x in range(0,len(rule)) if rule[x]==1]
        print("Class ",c,", Parent ",parent," original: ",len(labels_train[c]), " now: ", len(dataY))
        
        # Train Model
        if(modelS=="LogisticRegression"):
            if(len(list(np.unique(dataY)))==2):
                models[c].fit(dataX, dataY)
            else:
                dataY.pop()
                dataY.pop()
                dataY.append(0)
                dataY.append(1)
                models[c].fit(dataX, dataY)
        else:
            models[c].fit(dataX, dataY)
        print("...done")
    # Return Models
    return models

def predict_selStrat_binStruc(classes, features_test, models):
    predictions = {}
    for c in classes:
        predictions[c] = models[c].predict_proba(features_test)
    return predictions

def getArgument(args, title):
    for i in range(0, len(args)):
        if(args[i].startswith("-"+title)):
            if(i<len(args)-1):
                return args[i+1]
            else:
                return ""
    return ""

def getPredictionMatrix(Mprob, classes, levels, sNodes, nNodes, threshold):
    Mpred = list()
    for x in range(0,len(Mprob)):
        probs = Mprob[x]

        preds = list()
        for i in range(0,len(probs)):
            preds.append(0)
            
        maxLevel = max(levels)
        lastNode = ""
        for l in range(1, maxLevel+1): # for each level
            # determine candidates on that level
            candidates = list() 
            candidatesIdx = list()
            for i in range(0, len(classes)): 
                if(levels[i]==l):
                    if(lastNode==""):
                        candidates.append(classes[i])
                        candidatesIdx.append(i)
                    else:
                        if(lastNode in sNodes[i]):
                            candidates.append(classes[i])
                            candidatesIdx.append(i)
            
            # determine candidate with highest probability
            mx = -1
            ix = -1
            cx = -1
            for c in range(0, len(candidates)):
                if(probs[candidatesIdx[c]]>mx):
                    ix = candidatesIdx[c]
                    cx = c
                    mx = probs[candidatesIdx[c]]
                  
            # add prediction
            if(mx >= threshold):
                preds[ix] = 1
                lastNode = candidates[cx]
            else:
                break
           
        Mpred.append(preds)
    return Mpred

def getLevels(classes):
    levels = list()
    for c in classes:
        levels.append(len(c.split("/")))
    return levels

def getSuperiorNode(c):
    if(not "/" in c):
        return ""
    else:
        parts = c.split("/")
        newC = ""
        for p in range(0, len(parts)-1):
            newC = newC + parts[p] + "/"
        return newC[:-1]
    
def getSuperiorNodes(classes, levels):
    sNodes = list()
    for i in range(0,len(classes)):
        if(not "/" in classes[i]):
            sNodes.append(list())
        else:
            l = list()
            sup = getSuperiorNode(classes[i])
            while(sup!=""):
                l.append(sup)
                sup = getSuperiorNode(sup)
            sNodes.append(l)
    return sNodes

def getNeighboringNodes(classes, sClasses):
    nNodes = list()
    for i in range(0, len(classes)):
        l = list()
        for x in range(0, len(classes)):
            if(i!=x):
                if(sClasses[i]==sClasses[x]):
                    l.append(classes[x])
        nNodes.append(l)
    return nNodes

def getPredictionLabel(classes, names, Mpred):
    labels = list()
    for pred in Mpred:
        label = ""
        for i in range(0,len(pred)):
            if(pred[i]==1):
                lTemp = classes[i]
                if(label==""):
                    label = lTemp
                else:
                    if(len(label.split("/"))<len(lTemp.split("/"))):
                        label = lTemp
        labels.append(names[label].replace(" ","")+" "+label)
    return labels

def classifyTransposons(kmerFeatureFile, protFeatureFile, proteinDB, modelFile, fastaFile, outputFile):
    # Load model
    file   = open(modelFile,'rb')
    models =  pickle.load(file)#["model"]
    file.close()
    eThreshold  = models["eThreshold"]
    classes     = models["classes"]
    names       = models["names"]
    models      = models["model"]
    
     # Load data
    features = loadFeatures(proteinDB, protFeatureFile, kmerFeatureFile, eThreshold)
       
     # Parameters
    # eThreshold  = 5.0    
    levels = getLevels(classes)
    sNodes = getSuperiorNodes(classes,levels)
    nNodes = getNeighboringNodes(classes, sNodes)
       
    # Create predictions
    predictions = predict_selStrat_binStruc(classes, features, models)
    
    # Save PredictedLabel
    f = open(outputFile, "w+")
    f2 = open(fastaFile, "r")

    for i in range(0, len(features)):
        line = " "
        Mprob = list()
        for c in classes:
            if(len(predictions[c][i])==2):
                line = line + str(predictions[c][i][1])+" "
                Mprob.append(predictions[c][i][1])
            else:
                line = line + str(0)+" "
                Mprob.append(0)
        predLabel = getPredictionLabel(classes,names,getPredictionMatrix([Mprob], classes, levels, sNodes, nNodes, threshold=0.0))
        
        lineX = f2.readline()
        while(lineX!=""):
            if(lineX.startswith(">")):
                break
            else:
                lineX = f2.readline()
                
        f.write(lineX)
        f.write(predLabel[0]+" "+line+"\n")

    line = "\n\n#Explanation:\n#forecast forecast "
    for c in classes:
        line = line + c + " "
    f.write(line+"\n")
    
    f.close()
    f2.close()
    
def trainClassificationModel(classes, names, classes_sort, classes_parent, levels, featureFileKmer, featureFileProt, labelFile, filePDB, eThreshold, outputFile):
    classes_neighb = getNeighboringNodes(classes, getSuperiorNodes(classes, levels))
    # Load data
    features = loadFeatures(filePDB, featureFileProt, featureFileKmer, eThreshold)
    labels = {}
    for c in classes:
        labels[c] = loadBinaryLabels(labelFile, c)
    # Train Models
    features_train = features
    labels_train = {}
    for c in classes:
        labels_train[c] = labels[c]
    models = trainModels_selStrat_binStruc(classes_sort, classes_parent, classes_neighb, labels_train, features_train, "RandomForest", classes)
    # Save Models
    modelDic = {}
    modelDic["model"]      = models
    modelDic["eThreshold"] = eThreshold
    modelDic["modelType"]  = "RandomForest"
    modelDic["classes"]    = classes
    modelDic["names"]      = names
    with open(outputFile, "wb") as h:
        pickle.dump(modelDic, h, protocol=pickle.HIGHEST_PROTOCOL)

# from RFSB_TaxonomyLoader import loadTaxonomy
# labelFile           = "G:/CambridgeGenData/TE_DB/ALL_tiny.label"
# classes, names, classes_parent, classes_sort, levels = loadTaxonomy(labelFile, configFile="")    
# classes_neighb = getNeighboringNodes(classes, getSuperiorNodes(classes, levels))
