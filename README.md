# Transposon Classifier "RFSB"
Transposon classification tool for nucleotide sequence classification, providing classification, model training and prediction evaluation.

- **Mode 1: Transposon classification**
  - **Description**: Mode 1 classifies given transposon DNA sequences in a fastaFile and stores the predicted results into an outputPredictionFile.
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
 conda install -c derkevinriehl transposon_classifier_rfsb 
```
*Note: Otherwise you can find all source codes in this Github repository.*

## Usage of Mode 1: Transposon classification
```
transposon_classifier_RFSB -help
transposon_classifier_RFSB –mode classify -fastaFile demo1_seq.fasta –outputPredictionFile demo1_results.txt
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
transposon_classifier_RFSB –mode evaluate -predLabelFile demo2_predLabel.txt –trueLabelFile demo2_trueLabel.txt –outputPickleFile True
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
transposon_classifier_RFSB –mode trainModel -fastaFile demo3_transposonDB.fasta –labelFile demo3_labels.txt –outputModelFile demo3_model.pickle –eThreshold 5.0
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

## Explanation of output files
SX3351_addisababa.SV.vcf.gff3.matches.gff3

For each filtered structural variant a set of potential transposon annotation candidates (IDs similar to transposon annotation file) is reported:
```
seq1	PBSV	duplication	2909118	2910241	.	+	.	['23769', '23770', '23771'];Sseq1TYPE=DUP;END=2910240;Sseq1LEN=1122
seq1	PBSV	deletion	163800	164962	.	+	.	['1', '11827 '];Sseq1TYPE=DEL;END=164961;Sseq1LEN=-1161
seq1	PBSV	deletion	290360	290514	.	+	.	['11843 '];Sseq1TYPE=DEL;END=290513;Sseq1LEN=-153
seq1	PBSV	insertion	343890	344420	.	+	.	['538 '];merged;Sseq1TYPE=seq4NS;END=344424;Sseq1LEN=533
...
```

SX3351_addisababa.transpositionEvents.gff3

For each final structural variant that is considered to be a transposition event, the given transposon annotation (IDs similar to transposon annotation file) and predicted class are reported:
```
seq1	PBSV	deletion	290360	290514	.	+	.	Transposon=11843;Class=2/1/2(hAT,TIR,DNATransposon);Sseq1TYPE=DEL;END=290513;Sseq1LEN=-153
seq1	PBSV	insertion	610241	614786	.	+	.	Transposon=545;Class=2/1/3(CMC,TIR,DNATransposon);merged;merged;Sseq1TYPE=seq4NS;END=611763;Sseq1LEN=1521
seq1	PBSV	deletion	879772	884345	.	+	.	Transposon=556;Class=1/1/2(Gypsy,LTR,Retrotransposon);Sseq1TYPE=DEL;END=884344;Sseq1LEN=-4572
seq1	PBSV	insertion	1126531	1126860	.	+	.	Transposon=23592;Class=2/1/1(Tc1-Mariner,TIR,DNATransposon);Sseq1TYPE=seq4NS;END=1126859;Sseq1LEN=327
...
```

## Citations
Please cite our paper if you find transposition event detector "deTEct" useful:
(in progress)
