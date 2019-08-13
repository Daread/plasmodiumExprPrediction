# README
## Datasets
The datasets used in "Predicting gene expression in the human malaria parasite Plasmodium falciparum" are available within the _modelData_ folder.

There are .csv files with features/labels for the three stages that were analyzed in the paper. The _Response_ column at the end of each .csv is the expression label ("Low" or "High"), which was the value we sought to predict throughout our work. The rest of the columns are feature values that were used for making model predictions.

In our work, motif features (putatively giving information about transcription factor binding sites around a given gene) seemed to be unhelpful during our preliminary analysis. Consequently, our final models didn't use datasets containing any motif features. However, both versions of .csv files are included for those interested in exploring the dataset including motifs (all "*withMotifs.csv" files) or without motifs (all "*noMotifs.csv" files, which are what we used for final analysis).

Also, all feature values in these files are not converted to z scores. Our analysis made that transformation, but the values presented here did not make that transformation in case the original values are of interest.

Finally, "modelData" includes a text file detailing the train/test splits. For each stage, it lists the chromosomes that were used to form each of the five folds of data, as well as the rows corresponding to those folds (there are some differences of 1 gene between stages, due to the need to discard genes where at least one feature was unavailable at a particular gene). The rows correspond to the csv files attached, starting counting at 0. The "test folds" are the two folds that were used for final predictions of model accuracy (and which we never used during model development, up until the final stage of model evaluation). The three "training folds" were sub-divided into training/validation for cross-validation during preliminary analysis.

## Data Processing
Additional files in this repository represent the code used to process data and carry out analysis in "Predicting gene expression in the human malaria parasite Plasmodium falciparum"

Shell scripts run within the "dataPreProcessing" subdirectory used the command ```echo "source ~/myenv.sh" >> $(a).job;```, or similar commands. The file ~/myenv.sh contained instructions for which software/versions to load before running the specified job file. These instructions are specific to your computing environment.

Anytime ```source ~/myenv.sh``` is run, you need to have the following software/versions loaded:
* python 2.7.3
* numpy 1.11.0
* pandas 0.18.1
* samtools 0.1.19
* pysam 0.8.1
* cython 0.25.2
* matplotlib 1.5.1
* bwa 0.7.3
* fastqc 0.11.3
* fastx-toolkit 0.0.14
* pysam 0.8.4
* curl 7.48.0

A number of the scripts contained here require reference datasets that we have not directly included in this repository. To run these scripts locally, you will need:
* An indexed Plasmodium genome for use in BWA alignment.
  * Use BWA version 0.7.3
  * Replace all instances of _"/net/noble/vol1/data/bwa-indices/pfal3D7-PlasmoDBv29/pfal3D7-PlasmoDBv29"_ with the appropriate path to your local indexed Plasmodium genome.
* The plasmodium genome in fasta format. 
  * Replace all instances of _"/net/noble/vol1/data/reference_genomes/pfal3D7-PlasmoDBv29/PlasmoDB-29_Pfalciparum3D7_Genome.fasta"_ with the path to your local copy of the P. falciparum genome, v29 on PlasmoDB