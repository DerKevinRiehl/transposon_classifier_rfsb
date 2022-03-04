# Transposon Classifier "RFSB"
Transposon classification tool for nucleotide sequence classification, providing classification, model training and prediction evaluation. *RFSB* is part of [TransposonUltimate](https://github.com/DerKevinRiehl/TransposonUltimate).

- **Mode 1: Transposon classification**
  - **Description**: Mode 1 classifies given transposon DNA sequences in a fastaFile and stores the predicted results into an outputPredictionFile. (RFSB comes with a pretrained classification model ready to use)
  - **Input**:  Transposon nucleotide sequence(s) (FASTA File)
  - **Output**: Class predictions (and probabilities)
- **Mode 2: Prediction evaluation**
  - **Description**: Mode 2 evaluates given transposon classification results by a classification software predLabelFile onto the true class labels trueLabelFile and calculate results of the evaluation (accuracy, precision, recall, F1, MCC, false positive rate and the confusion matrix for all three perspectives). The results can be print to console, stored into CSV and PICKLE file, and in addition to that diagrams can be generated as well.
  - **Input**:  Prediction label file and true label file
  - **Output**: Classification performance evaluation
- **Mode 3: Model training**
  - **Description**: Mode 3 trains the proposed RFSB transposon classification model onto a given transposon DNA sequence database. The taxonomy can either be given or is automatically detected scanning the true label file corresponding to the database fasta file. Once trained, the model will store the taxonomic scheme and given e-threshold. The proteinDB necessary to run RPSTBLASTN however needs to be given, when running Mode 1.
  - **Input**:  Transposon database, true labels (optionally classification taxonomy)
  - **Output**: Trained classification model

## Installation
Installation as [CondaPackage](https://anaconda.org/DerKevinRiehl/transposon_classifier_rfsb):
```
 conda install -c derkevinriehl -c bioconda transposon_classifier_rfsb 
```
*Note: Otherwise you can find all source codes in this Github repository. Please extract the models out of the ZIP file in the "models" folder if you clone and use this Github.*

## Usage of Mode 1: Transposon classification
If you want to reproduce following examples, please download the "demoFiles" folder from this GitHub. Otherwise you can use your own FASTA files. (RFSB comes with a pretrained classification model ready to use)
```
transposon_classifier_RFSB -help
transposon_classifier_RFSB -mode classify -fastaFile demoFiles/demo1_seq.fasta -outputPredictionFile demoFiles/demo1_results.txt
```
Parameter | Mandatory | Description
------------ | ------------- | -------------
fastaFile | (mandatory) | Fasta file containing DNA sequences
outputPredictionFile | (mandatory) | Output prediction file containing classification results
modelFile | (optional) | File containing model for classification, default: *models/TE_Classifier_RFSB_models_ALL_Big.pickle*
kmerConfigFile | (optional) | Configuration file containing k-mers, default: *config/kmers.txt*
proteinDBFile | (optional) | Configuration file containing NCBI CDD PSSM model IDs, default: *config/RPSTBLASTN_LIB/db_large.pn*
tempFolder | (optional) | Folder in which temporary files are stored to, default: *(folder of the outputFile)*
deleteTempFiles | (optional) | Whether to delete temporary feature files, default: *True*
 
## Usage of Mode 2: Prediction evaluation
```
transposon_classifier_RFSB -mode evaluate -predLabelFile demoFiles/demo2_predLabel.txt -trueLabelFile demoFiles/demo2_trueLabel.txt -outputPickleFile True
```
Parameter | Mandatory | Description
------------ | ------------- | -------------
predLabelFile | (mandatory) | File containing predicted labels = output of MODE [1]
trueLabelFile | (mandatory) | File containing true labels
printResults | (optional) | Whether to print results into console, default: *True*
saveResultsCSV | (optional) | Whether to save results as CSV file, default: *True*
saveResultsPickle | (optional) | Whether to save output as .pickle file, default: *False*
outputPickleFile | (optional) | Desired Pickle output file containing statistics, default: *predLabelFile+"_results.pickle"*
outputCSVFile | (optional) | Desired CSV output file containing statistics, default: *predLabelFile+"_results.csv"*
generateDiagramPNG | (optional) | Whether to render diagram and save as PNG, default: *True*
generateDiagramSVG | (optional) | Whether to render diagram and save as SVG, default: *True*
diagTitle | (optional) | The rendered diagrams' title
outputPNGFile | (optional) | Desired PNG output file for rendered diagram, default: *predLabelFile+"_diagram.png"*
outputSVGFile | (optional) | Desired SVG output file for rendered diagram, default: *predLabelFile+"_diagram.svg"*
 
## Usage of Mode 3: Model training
```
transposon_classifier_RFSB -mode trainModel -fastaFile demoFiles/demo3_transposonDB.fasta -labelFile demoFiles/demo3_labels.txt -outputModelFile demoFiles/demo3_model.pickle -eThreshold 5.0
```
Parameter | Mandatory | Description
------------ | ------------- | -------------
fastaFile | (mandatory) | Fasta file containing DNA sequences
labelFile | (mandatory) | True class label file
outputModelFile | (mandatory) | Desired path to store output model
taxonomyConfigFile | (optional) | Configuration file containing the taxonomy scheme used, if not given, taxonomy will automatically be generated from labels in labelFile. (In config folder example available), default: *(none)*
proteinDBFile | (optional) | Configuration file containing NCBI CDD PSSM model IDs, default: *config/RPSTBLASTN_LIB/db_large.pn*
kmerConfigFile | (optional) | Configuration file containing k-mers, default: *config/kmers.txt*
eThreshold | (optional) | Default e-value threshold for RPSTBLASTN application of NCBI CDD protein models, default: *5.0*
tempFolder | (optional) | Folder in which temporary files are stored to, default: *(folder of the outputFile)*
deleteTempFiles | (optional) | Whether to delete temporary feature files, default: *True*

## Explanation of file structures

**Structure of *predLabelFile* resp. *outputPredictionFile***
After the header of each DNA sequence of the fasta file there are two fields separated by whitespace: transposon class name and position in hierarchy, followed by two whitespaces and then the whitespace-separated probabilities for each class of the hierarchy (between 0 and 1). After the last DNA sequence there needs to be two empty lines followed by the footer outlining the taxonomy as shown in the example.
```
>Transposon1
Gypsy,LTR,Retrotransposon 1/1/2  1.0 0.8 0.0 1.0 0.0 0.0 0.5 0.0 0.0 0.4 0.3 0.6 0.3 0.0 0.0 0.0 0.1 0.0 
>Transposon2
Gypsy,LTR,Retrotransposon 1/1/2  0.9 1.0 0.1 0.9 0.0 0.0 0.4 0.6 0.0 0.7 0.0 0.5 0.2 0.2 0.2 0.0 0.1 0.1 
...
>TransposonN
Copia,LTR,Retrotransposon 1/1/1  0.8 0.7 0.9 0.0 0.0 0.1 0.8 0.0 0.1 0.8 0.0 0.4 0.3 0.0 0.2 0.0 0.0 0.0 


#Explanation:
#forecast forecast 1 1/1 1/1/1 1/1/2 1/1/3 1/2 1/2/1 1/2/2 2 2/1 2/1/1 2/1/2 2/1/3 2/1/4 2/1/5 2/1/6 2/2 2/3
```

**Structure of *trueLabel* File**
After the header of each DNA sequence of the fasta file the transposons class position label (as defined in taxonomy) needs to be used (separated with slash).
```
>Transposon1
1/1/2
>Transposon2
2/1/2
>Transposon3
2/2
...
>TransposonN
2/1/3
```

**Structure of *taxonomyConfigFile***
For each class in taxonomy, create a single line with the class’ hierarchy position label, followed by colon, followed by a description. Please make sure to follow the order as shown here.
```
1:Class I, Retrotransposon
1/1:LTR, Retrotransposon
1/1/1:Copia, LTR, Retrotransposon
1/1/2:Gypsy, LTR, Retrotransposon
1/1/3:ERV, LTR, Retrotransposon
1/2:Non-LTR, Retrotransposon
1/2/1:LINE, Non-LTR, Retrotransposon
1/2/2:SINE, Non-LTR, Retrotransposon
2:Class II, DNA Transposon
2/1:TIR, DNA Transposon
2/1/1:Tc1-Mariner, TIR, DNA Transposon
2/1/2:hAT, TIR, DNA Transposon
2/1/3:CMC, TIR, DNA Transposon
2/1/4:Sola, TIR, DNA Transposon
2/1/5:Zator, TIR, DNA Transposon
2/1/6:Novosib, TIR, DNA Transposon
2/2:Helitron, DNA Transposon
2/3:MITE, DNA Transposon
```

**Structure of *kmerConfigFile***
In one line, separated by comma, list all kmer subsequences to scan for.
```
AA,AT,AC,...,TGGG,TGGT,TGGA,TTTC,TTTG,TTAT,TTCT,TTGT,TTGG,TTGA,TTGC,TTAA,TTAC,TTAG,TTCC,TTCA,TTCG
```

## Citations
Please cite our paper if you find TransposonUltimate useful:

Kevin Riehl, Cristian Riccio, Eric A Miska, Martin Hemberg, TransposonUltimate: software for transposon classification, annotation and detection, Nucleic Acids Research, 2022; gkac136, https://doi.org/10.1093/nar/gkac136

```
@article{riehl2022transposonultimate,
  title={TransposonUltimate: software for transposon classification, annotation and detection},
  author={Riehl, Kevin and Riccio, Cristian and Miska, Eric and Hemberg, Martin},
  journal={Nucleic Acids Research},
  year={2022}
}
```

## Acknowledgements
We would like to thank Sarah Buddle, Simone Procaccia, Fu Xiang Quah and Alexandra Dallaire for their assistance with testing and debugging the software.

## Benchmark Code
Please find the modified source codes and linux shell scripts that were used for benchmarking [here](https://github.com/DerKevinRiehl/transposon_classifier_rfsb/blob/main/benchmark/ClassifierCode.rar).
