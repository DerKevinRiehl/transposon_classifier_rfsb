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
 | | 
 | | 
 
## Usage of Mode 2: Prediction evaluation
```
transposon_classifier_RFSB –mode evaluate -predLabelFile demo2_predLabel.txt –trueLabelFile demo2_trueLabel.txt –outputPickleFile True
```
Parameter | Mandatory | Description
------------ | ------------- | -------------
 | | 
 | | 
 
## Usage of Mode 3: Model training
```
transposon_classifier_RFSB –mode trainModel -fastaFile demo3_transposonDB.fasta –labelFile demo3_labels.txt –outputModelFile demo3_model.pickle –eThreshold 5.0
```
Parameter | Mandatory | Description
------------ | ------------- | -------------
fastaFile | (mandatory) | Fasta file containing DNA sequences
labelFile | (mandatory) | 
outputModelFile | (mandatory) | 
taxonomyConfigFile | (optional) | 
proteinDBFile | (optional) | 
kmerConfigFile | (optional) | 
eThreshold | (optional) | 
tempFolder | (optional) | 
deleteTempFiles | (optional) | 

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
