
# fix_first_col:
# 	sed 's/^\(GeneDB|\)*//' <Nucleosome_18h_gt_200.bedgraph >fixed18h_gt_200.bedgraph
# 	mv fixed18h_gt_200.bedgraph Nucleosome_18h_gt_200.bedgraph
# 	sed 's/^\(GeneDB|\)*//' <Nucleosome_36h_gt_200.bedgraph >fixed36h_gt_200.bedgraph
# 	mv fixed36h_gt_200.bedgraph Nucleosome_36h_gt_200.bedgraph
# 	sed 's/^\(GeneDB|\)*//' <Nucleosome_18h_lt_100.bedgraph >fixed18h_lt_100.bedgraph
# 	mv fixed18h_lt_100.bedgraph Nucleosome_18h_lt_100.bedgraph
# 	sed 's/^\(GeneDB|\)*//' <Nucleosome_36h_lt_100.bedgraph >fixed36h_lt_100.bedgraph
# 	mv fixed36h_lt_100.bedgraph Nucleosome_36h_lt_100.bedgraph
# 	sed 's/^\(GeneDB|\)*//' <Nucleosome_18h_100_to_200.bedgraph >fixed18h_100_to_200.bedgraph
# 	mv fixed18h_100_to_200.bedgraph Nucleosome_18h_100_to_200.bedgraph
# 	sed 's/^\(GeneDB|\)*//' <Nucleosome_36h_100_to_200.bedgraph >fixed36h_100_to_200.bedgraph
# 	mv fixed36h_100_to_200.bedgraph Nucleosome_36h_100_to_200.bedgraph

DATADIR = /net/noble/vol2/home/katecook/proj/2016predictExpression/data
GENOMEDATA_FILE = $(DATADIR)/pfal3D7.genomedata

GENOME = /net/noble/vol1/data/reference_genomes/pfal3D7-PlasmoDBv29/PlasmoDB-29_Pfalciparum3D7_Genome.fasta

BUNNIK_2014 = /net/gs/vol1/home/readdf/nobleLabProj/2018_readdf_predict-expression/results/2018_04_16_MNaseReprocessing

#all: load_bunnik2014
all: replace_bunnik2014

# load_bunnik2014:
# 	genomedata-open-data $(GENOMEDATA_FILE) --trackname bunnik2014_MNase_lt_100_18h bunnik2014_MNase_lt_100_36h bunnik2014_MNase_100_to_200_18h bunnik2014_MNase_100_to_200_36h bunnik2014_MNase_gt_200_18h bunnik2014_MNase_gt_200_36h
# 	genomedata-load-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_lt_100_18h < $(BUNNIK_2014)/Nucleosome_18h_lt_100.bedgraph
# 	genomedata-load-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_lt_100_36h < $(BUNNIK_2014)/Nucleosome_36h_lt_100.bedgraph
# 	genomedata-load-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_100_to_200_18h < $(BUNNIK_2014)/Nucleosome_18h_100_to_200.bedgraph
# 	genomedata-load-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_100_to_200_36h < $(BUNNIK_2014)/Nucleosome_36h_100_to_200.bedgraph
# 	genomedata-load-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_gt_200_18h < $(BUNNIK_2014)/Nucleosome_18h_gt_200.bedgraph
# 	genomedata-load-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_gt_200_36h < $(BUNNIK_2014)/Nucleosome_36h_gt_200.bedgraph
# 	genomedata-close-data $(GENOMEDATA_FILE)

replace_bunnik2014:
	genomedata-erase-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_lt_100_18h 
	genomedata-erase-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_lt_100_36h 
	genomedata-erase-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_100_to_200_18h 
	genomedata-erase-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_100_to_200_36h 
	genomedata-erase-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_gt_200_18h 
	genomedata-erase-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_gt_200_36h 
	genomedata-load-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_lt_100_18h < $(BUNNIK_2014)/Nucleosome_18h_lt_100.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_lt_100_36h < $(BUNNIK_2014)/Nucleosome_36h_lt_100.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_100_to_200_18h < $(BUNNIK_2014)/Nucleosome_18h_100_to_200.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_100_to_200_36h < $(BUNNIK_2014)/Nucleosome_36h_100_to_200.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_gt_200_18h < $(BUNNIK_2014)/Nucleosome_18h_gt_200.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) bunnik2014_MNase_gt_200_36h < $(BUNNIK_2014)/Nucleosome_36h_gt_200.bedgraph
	genomedata-close-data $(GENOMEDATA_FILE)

