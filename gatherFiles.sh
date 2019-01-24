#1 /bin/bash
# This script, when run, will go through the directories where
#   analysis was carried out for Read et al, "Predicting gene expression in the human
#		 malaria parasite Plasmodium falciparum"
# 	Analysis was originally carried out in several locations across the Genome
# 	Sciences computing cluster directory tree. When run, this shell script
# 	will copy current versions (as they exist on the GS cluster) into a new
#	directory, to make a centralized location that will be shared as a public
#	github, accompnaying the Plasmodium manuscript

# Get environment
cp /net/noble/vol2/katecook/myenv.sh ./myenv.sh

# ChIP-Seq Data Pre-processing
# mkdir ./dataPreProcessing
# Bartfai et al data handling
mkdir ./dataPreProcessing/BartfaiEtAl
bartfaiDir="/net/trapnell/vol1/home/readdf/plasmodiumExprPred/results/katecook/20161115_bartfai2010_histone_data"
cp $bartfaiDir/Makefile ./dataPreProcessing/BartfaiEtAl/Makefile
cp $bartfaiDir/accession_list.txt ./dataPreProcessing/BartfaiEtAl/
cp $bartfaiDir/info.tab ./dataPreProcessing/BartfaiEtAl/
mkdir ./dataPreProcessing/BartfaiEtAl/Lib
cp $bartfaiDir/Lib/* ./dataPreProcessing/BartfaiEtAl/Lib/
# jiang et al data handling
mkdir ./dataPreProcessing/jiangEtAl
jiangDir="/net/trapnell/vol1/home/readdf/plasmodiumExprPred/results/katecook/20161117_jiang2013_histone_data"
cp $jiangDir/Makefile ./dataPreProcessing/jiangEtAl/Makefile
cp $jiangDir/accession_list.txt ./dataPreProcessing/jiangEtAl/
cp $jiangDir/info.tab ./dataPreProcessing/jiangEtAl/
mkdir ./dataPreProcessing/jiangEtAl/Lib
cp $jiangDir/Lib/* ./dataPreProcessing/jiangEtAl/Lib/
# bunnik et al data handling
mkdir ./dataPreProcessing/bunnikEtAl
bunnikDir="/net/trapnell/vol1/home/readdf/plasmodiumExprPred/results/katecook/20161119_bunnik2014_MNase_data"
cp $bunnikDir/Makefile ./dataPreProcessing/bunnikEtAl/Makefile
cp $bunnikDir/accession_list.txt ./dataPreProcessing/bunnikEtAl/
cp $bunnikDir/info.tab ./dataPreProcessing/bunnikEtAl/
mkdir ./dataPreProcessing/bunnikEtAl/Lib
cp $bunnikDir/Lib/* ./dataPreProcessing/bunnikEtAl/Lib/

# Get GC Assignments
mkdir ./dataPreProcessing/GCcontent
GCdir="/net/trapnell/vol1/home/readdf/plasmodiumExprPred/results/katecook/20170105_GC_content"
cp $GCdir/Makefile ./dataPreProcessing/GCcontent/Makefile
cp $GCdir/ids.lst ./dataPreProcessing/GCcontent/
mkdir ./dataPreProcessing/GCcontent/Lib
cp $GCdir/Lib/* ./dataPreProcessing/GCcontent/Lib/

# Get Hi-C Data
mkdir ./dataPreProcessing/HiCdata
# Download from website https://noble.gs.washington.edu/proj/plasmo3d/, where processed Hi C Data is held from Ay et al, 2014.
cd ./dataPreProcessing/HiCdata
wget https://noble.gs.washington.edu/proj/plasmo3d/ourData/Genes-Positions/RINGS.genes_distsToLandmarks.txt
wget https://noble.gs.washington.edu/proj/plasmo3d/ourData/Genes-Positions/TROPHOZOITES-XL.genes_distsToLandmarks.txt
wget https://noble.gs.washington.edu/proj/plasmo3d/ourData/Genes-Positions/SCHIZONTS.genes_distsToLandmarks.txt
cd ../..

# Motif scans
mkdir ./dataPreProcessing/motifData
motifDir="/net/trapnell/vol1/home/readdf/plasmodiumExprPred/results/katecook/20170424_scan_with_pfms"
cp $motifDir/Makefile ./dataPreProcessing/motifData/Makefile
cp $motifDir/ids.lst ./dataPreProcessing/motifData/
cp $motifDir/info.tab ./dataPreProcessing/motifData/
cp $motifDir/fimo.background ./dataPreProcessing/motifData/
mkdir ./dataPreProcessing/motifData/Lib
cp $motifDir/Lib/* ./dataPreProcessing/motifData/Lib/

# Gather files for loading genome data
mkdir ./dataPreProcessing/dataLoading
# Load the main set of ChIP-seq data
mkdir ./dataPreProcessing/dataLoading/mainChIPseq
mainChipDir="/net/trapnell/vol1/home/readdf/plasmodiumExprPred/results/katecook/20161115_genomedata"
cp $mainChipDir/Makefile ./dataPreProcessing/dataLoading/mainChIPseq
cp $mainChipDir/info.tab ./dataPreProcessing/dataLoading/mainChIPseq
# Loading the motif data
mkdir ./dataPreProcessing/dataLoading/motifData
motifLoadingDir="/net/trapnell/vol1/home/readdf/nobleLabProj/2018_readdf_predict-expression/results/2018_05_02_motif_scans"
cp $motifLoadingDir/Makefile ./dataPreProcessing/dataLoading/motifData
cp $motifLoadingDir/noRedundancyMotifStageToTrack.csv ./dataPreProcessing/dataLoading/motifData
# Re-load the MNase data
mkdir ./dataPreProcessing/dataLoading/mnaseData
mnaseLoadDir="/net/trapnell/vol1/home/readdf/nobleLabProj/2018_readdf_predict-expression/results/2018_04_16_MNaseReprocessing"
cp $mnaseLoadDir/Makefile ./dataPreProcessing/dataLoading/mnaseData
mkdir ./dataPreProcessing/dataLoading/mnaseData/Lib
cp $mnaseLoadDir/Lib/*	./dataPreProcessing/dataLoading/mnaseData/Lib

# Loading the data into a single csv file in preparation for model training
mkdir ./dataPreProcessing/loadIntoCSV
csvLoadDir="/net/trapnell/vol1/home/readdf/nobleLabProj/2018_readdf_predict-expression/results/2018_04_30_feature_revisit"
# Get all the variables from motifsIncludedFeatureToTrack.csv that aren't noNorm tracks. We ended up not using those, and
#	I think it'll be confusing to include them at all (even just to load and not use). So I'll cut those here
removeCols='6,7,8,15,16'
cat $csvLoadDir/motifsIncludedFeatureToTrack.csv | cut -d, -f$removeCols --complement > ./dataPreProcessing/loadIntoCSV/motifsIncludedFeatureToTrack.csv
cp $csvLoadDir/assignFeatures.py ./dataPreProcessing/loadIntoCSV
mkdir ./dataPreProcessing/loadIntoCSV/assignedCSVfiles




















