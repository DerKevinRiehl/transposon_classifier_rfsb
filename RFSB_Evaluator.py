############################################################################
##### Transposon Classifier RFSB - part of Transposon Ultimate #############
##### Kevin Riehl (kevin.riehl.de@gmail.com, 2021) #########################
############################################################################

# Imports
import matplotlib.pyplot as plt
import numpy as np
import pickle

# Methods
def loadLabelData(predLabel, trueLabel): 
    # Load classes
    predL = list()
    f = open(predLabel, "r")
    predictions = f.readlines()
    f.close()
    classes = predictions[-1].replace(" \n","").split(" ")[2:]
    # load predicted labels
    predictions = predictions[:-4]
    for i in range(0,len(predictions)):
        if(predictions[i].startswith(">")):
            pass
        else:
            predictions[i] = predictions[i].replace(" \n","").split(" ")[3:]
            # print(predictions[i])
            for x in range(0,len(predictions[i])):
                predictions[i][x] = float(predictions[i][x])
            predL.append(predictions[i])
    # load true labels
    trueL = list()
    f = open(trueLabel)
    truelabels = f.readlines()
    f.close()
    for i in range(0,len(truelabels)):
        if(truelabels[i].startswith(">")):
            pass
        else:
            label = truelabels[i].replace("\n","").replace(".","/")
            binL = list()
            for c in classes:
                if(label.startswith(c)):
                    binL.append(1)
                else:
                    binL.append(0)
            trueL.append(binL)
    return classes, predL, trueL
  
def getRandomData(path, randType="equalDist", folds=10):
    # Load Data
    classes = "1,1/1,1/1/1,1/1/2,1/1/3,1/2,1/2/1,1/2/2,2,2/1,2/1/1,2/1/2,2/1/3,2/1/4,2/1/5,2/1/6,2/2,2/3".split(",")
    
    predL = list()
    trueL = list()
    for z in range(1,folds+1):
        f = open(path+"/inference"+str(z)+"/truelabels.txt")#open("advanced/truelabels.txt")
        truelabels = f.readlines()
        predictions = list()
        f.close()
        for i in range(0,len(truelabels)):
            truelabels[i] = truelabels[i].split(" ")[:-1]
            pL = list()
            for x in range(0,len(truelabels[i])):
                truelabels[i][x] = int(truelabels[i][x])
                if(randType=="equalDist"):
                    pL.append(np.random.rand())
                else: # otherwise normal distribution mu=0.5, std=0.5
                    pL.append(np.random.normal(loc = 0.5, scale = 0.5))
            predictions.append(pL)
        trueL.append(truelabels)
        predL.append(predictions)
        
    return classes, predL, trueL
      
def getThresholdValues():
    thresholdValues = list()
#    thresholdValues.append(1.01)
    for z in range(100,-1,-1):
        thresholdValues.append(float(z)/100)
#    thresholdValues.append(-0.01)
    return thresholdValues

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
    
def getSuperiorNodesIncl(classes, levels):
    sNodes = getSuperiorNodes(classes, levels)
    for i in range(0,len(classes)):
        sNodes[i].insert(0,classes[i])
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

def getNeighboringNodesIncl(classes, levels):
    nNodes = getNeighboringNodes(classes, levels)
    for i in range(0,len(classes)):
        nNodes[i].insert(0,classes[i])
    return nNodes

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

def getPerformanceMeasures(confusionMat):
    [TP, FP, TN, FN, n] = confusionMat
    
    # numerical issues
#    if(TP==0):
#        TP += 0.1
#    if(FP==0):
#        FP += 0.1
#    if(TN==0):
#        TN += 0.1
#    if(FN==0):
#        FN += 0.1
    
    if(TP+TN+FP+FN!=0):
        ACC = (TP+TN)/(TP+TN+FP+FN)
    else:
        ACC = 0
    if(TP+FP!=0):
        PPV = (TP)/(TP+FP)
    else:
        PPV = 0
    if(TP+FN!=0):
        TPR = (TP)/(TP+FN)
    else:
        TPR = 0
    if(2*TP+FP+FN):
        F1 = (2*TP)/(2*TP+FP+FN)
    else:
        F1 = 0
    if((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)!=0):
        MCC = (TP*TN - FP*FN) / np.sqrt(TP+FP)  / np.sqrt(TP+FN)  / np.sqrt(TN+FP)  / np.sqrt(TN+FN)
    else:
        MCC = 0
    if(FP+TN!=0):
        FPR = (FP)/(FP+TN)
    else:
        FPR = 0
    return [ACC, PPV, TPR, F1, MCC, FPR]

def getLabels(M, classes, levels):
    labelL = list()
    levelL = list()
    clIdxL = list()
    for i in range(0,len(M)):
        idx = -1
        level = -1
        for c in range(0,len(classes)):
            if(M[i][c]==1):
                if(levels[c]>level):
                    level = levels[c]
                    idx = c
        labelL.append(classes[idx])
        levelL.append(level)
        clIdxL.append(idx)
    return labelL, levelL, clIdxL

def getConfusionMatrix_Simple(trueL, predL):
    TP = trueL*predL
    FP = (1-trueL)*predL
    TN = (1-trueL)*(1-predL)
    FN = trueL*(1-predL)
    return TP, FP, TN, FN
 
def getConfusionMatrix_Perspective1(Mtrue, Mpred, levels, nNodesIncl, Lbltrue, Lvltrue, c=1):    
    TP = 0 # counter variables
    FP = 0
    TN = 0
    FN = 0
    n  = 0
    j = c  # main loop 
    for i in range(0, len(Mtrue)):
        if(not (Lvltrue[i] >= levels[c])): # conditions that i in D_c set
            continue
        else:
            found = False
            for x in range(0,len(nNodesIncl[c])):
                if(Lbltrue[i].startswith(nNodesIncl[c][x])):
                    found = True
                    break
            if(not (found)):
                continue

        TPi, FPi, TNi, FNi = getConfusionMatrix_Simple(Mtrue[i][j], Mpred[i][j])
        TP += TPi
        FP += FPi
        TN += TNi
        FN += FNi
        n += 1
    return [TP, FP, TN, FN, n]
 
def getConfusionMatrix_Perspective2(Mtrue, Mpred, Lvltrue, levels, l=1):    
    TP = 0 # counter variables
    FP = 0
    TN = 0
    FN = 0
    n  = 0
    for i in range(0, len(Mtrue)): # main loop
        if(not (Lvltrue[i] >= l)): # conditions that i in D_z set
            continue
        for j in range(0, len(Mtrue[0])):
            if(not (levels[j] >= l)): # condtion that j in N_z
                continue
            TPi, FPi, TNi, FNi = getConfusionMatrix_Simple(Mtrue[i][j], Mpred[i][j])
            TP += TPi
            FP += FPi
            TN += TNi
            FN += FNi
            n += 1
    return [TP, FP, TN, FN, n]

def getConfusionMatrix_Perspective3(Mtrue, Mpred, Lbltrue, Lblpred, classes, levels, nNodes):
    TP = 0 # counter variables
    FP = 0
    TN = 0
    FN = 0
    sNodesIncl = getSuperiorNodesIncl(classes, levels)
    for i in range(0, len(Mtrue)): # main loop   
        pathTrue = sNodesIncl[classes.index(Lbltrue[i])]
        pathPred = sNodesIncl[classes.index(Lblpred[i])]
        if(len(pathTrue)<len(pathPred)):
            pathPred = pathPred[len(pathPred)-len(pathTrue):] # making  longer prediction Path to same length as pathTrue
        
        pathCommon = intersection(pathTrue, pathPred)
        childrenDeepest = cut(getChildren(getDeepestNode(pathCommon, classes, levels),classes),union(pathTrue,pathPred))
        neighborsOfCommonPath = list()
        for cl in pathCommon:
            neighborsOfCommonPath += nNodes[classes.index(cl)]       
        TP += power(pathCommon)
        FP += power(pathPred) - power(pathCommon)
        TN += power(neighborsOfCommonPath) + power(childrenDeepest)
        FN += power(pathTrue) - power(pathCommon)
    return [TP, FP, TN, FN, len(Mtrue)]

def union(A,B):
    return list(set(A+B))

def intersection(A,B):
    return list(set(A) & set(B)) 
    
def power(A):
    return len(A)

def getDeepestNode(A, classes, levels):
    lMax = -1
    node = ""
    for x in A:
        if(levels[classes.index(x)]>lMax):
            lMax = levels[classes.index(x)]
            node = x
    return node

def getChildren(node, classes):
    cList = list()
    for c in classes:
        if(c.startswith(node+"/")):
            cList.append(c)
    return cList

def cut(A,B):
    return list(set(A).difference(set(B)))

def calcStatistics(results, X, val):
    meaL = list()
    stdL = list()
    maxL = list()
    minL = list()
    for x in range(0,len(X)):
        dataL = list()
        for k in range(0,len(results)):
            dataL.append(results[k][val][x])
        meaL.append(np.average(dataL))
        stdL.append(np.std(dataL))
        maxL.append(np.max(dataL))
        minL.append(np.max(dataL))
    return meaL, stdL, maxL, minL
  
def calcFinalStatistics1(Mprob, Mtrue, classes, levels, threshold, c):
    Lbltrue, Lvltrue, Cidxtrue = getLabels(Mtrue, classes, levels)    
    sNodes = getSuperiorNodes(classes,levels)
    nNodes = getNeighboringNodes(classes, sNodes)
    nNodesIncl = getNeighboringNodesIncl(classes, sNodes)       
    Mpred = getPredictionMatrix(Mprob, classes, levels, sNodes, nNodes, threshold)
    Lblpred, Lvlpred, Cidxpred = getLabels(Mpred, classes,levels)
    confMat = getConfusionMatrix_Perspective1(Mtrue, Mpred, levels, nNodesIncl, Lbltrue, Lvltrue, c)
    measures = getPerformanceMeasures(confMat)
    return measures, confMat
     
def calcFinalStatistics2(Mprob, Mtrue, classes, levels, threshold, l):
    Lbltrue, Lvltrue, Cidxtrue = getLabels(Mtrue, classes, levels)    
    sNodes = getSuperiorNodes(classes,levels)
    nNodes = getNeighboringNodes(classes, sNodes)
    Mpred = getPredictionMatrix(Mprob, classes, levels, sNodes, nNodes, threshold)
    Lblpred, Lvlpred, Cidxpred = getLabels(Mpred, classes,levels)
    confMat = getConfusionMatrix_Perspective2(Mtrue, Mpred, Lvltrue, levels, l)
    measures = getPerformanceMeasures(confMat)
    return measures, confMat

def calcFinalStatistics3(Mprob, Mtrue, classes, levels, threshold):
    Lbltrue, Lvltrue, Cidxtrue = getLabels(Mtrue, classes, levels)    
    sNodes = getSuperiorNodes(classes,levels)
    nNodes = getNeighboringNodes(classes, sNodes)
    Mpred = getPredictionMatrix(Mprob, classes, levels, sNodes, nNodes, threshold)
    Lblpred, Lvlpred, Cidxpred = getLabels(Mpred, classes,levels)
    confMat = getConfusionMatrix_Perspective3(Mtrue, Mpred, Lbltrue, Lblpred, classes, levels, nNodes)
    measures = getPerformanceMeasures(confMat)
    return measures, confMat

def algorithmEvaluationSummary(predLabel, trueLabel, printResults, saveCSV, savePickle, genDiagPNG, genDiagSVG, diagTitle, csvFile, pickleFile, svgFile, pngFile):   
    classes, MprobL, MtrueL = loadLabelData(predLabel, trueLabel)
    levels = getLevels(classes)
    print("data loaded...")
    
    # Calculate Perspective 3 results
    thresholdValues = getThresholdValues()
    results = list()
    Mprob = MprobL
    Mtrue = MtrueL
    Lbltrue, Lvltrue, Cidxtrue = getLabels(Mtrue, classes, levels)    
    sNodes = getSuperiorNodes(classes,levels)
    nNodes = getNeighboringNodes(classes, sNodes)
    X = list() # threshold
    Y_ACC = list()
    Y_PPV = list()
    Y_TPR = list()
    Y_F1  = list()
    Y_MCC = list()
    Y_FPR = list()
    Y_conf = list()
    for threshold in thresholdValues:
        # print(threshold)
        Mpred = getPredictionMatrix(Mprob, classes, levels, sNodes, nNodes, threshold)
        Lblpred, Lvlpred, Cidxpred = getLabels(Mpred, classes,levels)
        confMat = getConfusionMatrix_Perspective3(Mtrue, Mpred, Lbltrue, Lblpred, classes, levels, nNodes)
        measures = getPerformanceMeasures(confMat)
        X.append(threshold)
        Y_ACC.append(measures[0])
        Y_PPV.append(measures[1])
        Y_TPR.append(measures[2])
        Y_F1.append (measures[3])
        Y_MCC.append(measures[4])
        Y_FPR.append(measures[5])
        Y_conf.append(confMat)
    dataSum = {}
    dataSum["X"] = X
    dataSum["ACC"] = Y_ACC    
    dataSum["PPV"] = Y_PPV
    dataSum["TPR"] = Y_TPR
    dataSum["F1"]  = Y_F1
    dataSum["MCC"] = Y_MCC
    dataSum["FPR"] = Y_FPR
    dataSum["CONF"] = Y_conf
    results.append(dataSum)
    print("done...start calcStatistics")
    
    # Calcualte Statistics
    X = results[0]["X"]
    accM, accS, accMx, accMi = calcStatistics(results, X, "ACC")
    ppvM, ppvS, ppvMx, ppvMi = calcStatistics(results, X, "PPV")
    tprM, tprS, tprMx, tprMi = calcStatistics(results, X, "TPR")
    f1sM, f1sS, f1sMx, f1sMi = calcStatistics(results, X, "F1")
    mccM, mccS, mccMx, mccMi = calcStatistics(results, X, "MCC")
    fprM, fprS, fprMx, fprMi = calcStatistics(results, X, "FPR")
    
    # Determine threshold as maximizing MCC
    threshold = X[mccM.index(max(mccM))]
            
    # Determine Final Results for given Threshold
    measureLabels = ["ACC", "PPV", "TPR", "F1", "MCC", "FPR"]
    finalResults = {}
    finalResults["threshold"] = threshold
    finalResults["PERFORMANCE"] = {}
    finalResults["CONFMAT"] = {}
    for i in range(0,len(classes)):
        finalResults["PERFORMANCE"]["Class "+classes[i]] = {}
        res, cres = calcFinalStatistics1(MprobL, MtrueL, classes, levels, threshold, i)
        for x in range(0,len(measureLabels)):
            finalResults["PERFORMANCE"]["Class "+classes[i]][measureLabels[x]] = res[x]         
        finalResults["CONFMAT"]["Class "+classes[i]] = cres
        
    for i in range(1,max(levels)+1):
        finalResults["PERFORMANCE"]["Level "+str(i)] = {}
        res, cres = calcFinalStatistics2(MprobL, MtrueL, classes, levels, threshold, i)
        for x in range(0,len(measureLabels)):
            finalResults["PERFORMANCE"]["Level "+str(i)][measureLabels[x]] = res[x]    
        finalResults["CONFMAT"]["Level "+str(i)] = cres
        
    finalResults["PERFORMANCE"]["Total"] = {}
    res, cres = calcFinalStatistics3(MprobL, MtrueL, classes, levels, threshold)
    for x in range(0,len(measureLabels)):
        finalResults["PERFORMANCE"]["Total"][measureLabels[x]] = res[x]         
    finalResults["CONFMAT"]["Total"] = cres
    
    finalResults["graphs-total"] = {}
    finalResults["graphs-total"]["X-thresholds"] = X
    finalResults["graphs-total"]["ACC"] = accM
    finalResults["graphs-total"]["PPV"] = ppvM
    finalResults["graphs-total"]["TPR"] = tprM
    finalResults["graphs-total"]["F1"] = f1sM
    finalResults["graphs-total"]["MCC"] = mccM
    finalResults["graphs-total"]["FPR"] = fprM
    
    # Save Final Results to Pickle File
    if(savePickle):
        finResultDic = {}
        finResultDic["finalResults"] = finalResults
        with open(pickleFile, 'wb') as handle:
            pickle.dump(finResultDic, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    # Print console results
    if(printResults):
        measures = ["ACC","PPV","TPR","F1","MCC","FPR"]
        print("Final results - performance measures")
        print("Prsp.\tIdx\tIdx\tACC\tPPV\tTPR\tF1\tMCC\tFPR\tTP\tFP\tTN\tFN\tN")
        # Perspective 1
        for i in range(0,len(classes)):
            print("1\tClass "+classes[i]+"\t", end="")
            for m in measures:
                print("{0:0.3f}".format(finalResults["PERFORMANCE"]["Class "+classes[i]][m]), "\t", end="")
            for x in range(0,5):
                print(str(finalResults["CONFMAT"]["Class "+classes[i]][x]), "\t", end="")
            print("\n")
        # Perspective 2
        for l in range(1,max(levels)+1):
            print("2\tLevel "+str(l)+"\t", end="")
            for m in measures:
                print("{0:0.3f}".format(finalResults["PERFORMANCE"]["Level "+str(l)][m]), "\t", end="")
            for x in range(0,5):
                print(str(finalResults["CONFMAT"]["Level "+str(l)][x]), "\t", end="")
            print("\n")
        # Perspective 3
        print("3\tTotal\tTotal\t", end="")
        for m in measures:
            print("{0:0.3f}".format(finalResults["PERFORMANCE"]["Total"][m]), "\t", end="")
        for x in range(0,5):
            print(str(finalResults["CONFMAT"]["Total"][x]), "\t", end="")
        print("\n")
    if(saveCSV):
        f = open(csvFile,"w+")
        measures = ["ACC","PPV","TPR","F1","MCC","FPR"]
        f.write("Final results - performance measures")
        f.write("\n")
        f.write("Prsp.\tIdx\tIdx\tACC\tPPV\tTPR\tF1\tMCC\tFPR\tTP\tFP\tTN\tFN\tN")
        f.write("\n")
        # Perspective 1
        for i in range(0,len(classes)):
            f.write("1\tClass "+classes[i]+"\t")
            for m in measures:
                f.write(str(finalResults["PERFORMANCE"]["Class "+classes[i]][m])+"\t")
            for x in range(0,5):
                f.write(str(finalResults["CONFMAT"]["Class "+classes[i]][x])+"\t")
            f.write("\n")
        # Perspective 2
        for l in range(1,max(levels)+1):
            f.write("2\tLevel "+str(l)+"\t")
            for m in measures:
                f.write(str(finalResults["PERFORMANCE"]["Level "+str(l)][m])+"\t")
            for x in range(0,5):
                f.write(str(finalResults["CONFMAT"]["Level "+str(l)][x])+"\t")
            f.write("\n")
        # Perspective 3
        f.write("3\tTotal\tTotal\t")
        for m in measures:
            f.write(str(finalResults["PERFORMANCE"]["Total"][m])+"\t")
        for x in range(0,5):
            f.write(str(finalResults["CONFMAT"]["Total"][x])+"\t")
        f.write("\n")
        f.close()
        
    # Plot graphics for Perspective 3 and save to PNG and SVG File
    plt.figure(figsize=(10,10))
    plt.plot(X,accM, label="Accuracy")
    plt.plot(X,ppvM, label="Precision")
    plt.plot(X,tprM, label="Recall")
    plt.plot(X,f1sM, label="F1")
    plt.plot(X,mccM, label="MCC")
    plt.plot(X,fprM, label="FPR")
    
    plt.plot([threshold, threshold], [0, 1], color="red", alpha=0.5)
    plt.fill_between(X, np.asarray(accM)+np.asarray(accS), np.asarray(accM)-np.asarray(accS), facecolor='blue', alpha=0.25)
    plt.fill_between(X, np.asarray(ppvM)+np.asarray(ppvS), np.asarray(ppvM)-np.asarray(ppvS), facecolor='blue', alpha=0.25)
    plt.fill_between(X, np.asarray(tprM)+np.asarray(tprS), np.asarray(tprM)-np.asarray(tprS), facecolor='blue', alpha=0.25)
    plt.fill_between(X, np.asarray(f1sM)+np.asarray(f1sS), np.asarray(f1sM)-np.asarray(f1sS), facecolor='blue', alpha=0.25)
    plt.fill_between(X, np.asarray(mccM)+np.asarray(mccS), np.asarray(mccM)-np.asarray(mccS), facecolor='blue', alpha=0.25)
    plt.fill_between(X, np.asarray(fprM)+np.asarray(fprS), np.asarray(fprM)-np.asarray(fprS), facecolor='blue', alpha=0.25)
    plt.legend()
    plt.xlim(-0.1,1.1)
    plt.ylim(-0.1,1.1)
    plt.gca().set_aspect("equal") # force square grid
    plt.grid()
    plt.title(diagTitle)
    plt.xlabel("Threshold")
    plt.ylabel("Classification Performance")
    #plt.show()
    if(genDiagPNG):
        plt.savefig(pngFile, format="png")
    if(genDiagSVG):
        plt.savefig(svgFile, format="svg")
    plt.close()
