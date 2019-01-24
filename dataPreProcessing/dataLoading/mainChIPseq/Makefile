
DATADIR = $(HOME)/proj/2016predictExpression/data
GENOMEDATA_FILE = $(DATADIR)/pfal3D7.genomedata

GENOME = /net/noble/vol1/data/reference_genomes/pfal3D7-PlasmoDBv29/PlasmoDB-29_Pfalciparum3D7_Genome.fasta

JIANG_2013 = $(HOME)/proj/2016predictExpression/results/katecook/20161117_jiang2013_histone_data
JIANG_MEAN = $(HOME)/proj/2016predictExpression/results/katecook/20170613_jiang2013_mean_bg
BUNNIK_2014 = $(HOME)/proj/2016predictExpression/results/katecook/20161119_bunnik2014_MNase_data
BARTFAI_2010 = $(HOME)/proj/2016predictExpression/results/katecook/20161115_bartfai2010_histone_data
GC =  $(HOME)/proj/2016predictExpression/results/katecook/20170105_GC_content

all: load_seq load_jiang2013 load_jiang_mean_bg load_bunnik2014 load_bartfai2010 load_GC

delete:
	rm -rf $(GENOMEDATA_FILE)


load_seq:
	genomedata-load-seq $(GENOMEDATA_FILE) $(GENOME)

load_jiang2013:
	genomedata-open-data $(GENOMEDATA_FILE) --tracknames jiang2013_H3K36me2_18h jiang2013_H3K36me2_42h jiang2013_H3K36me3_18h jiang2013_H3K4me3_18h jiang2013_H3K36me3_42h jiang2013_H3K4me3_42h jiang2013_H3K9me3_18h jiang2013_H3K9me3_42h jiang2013_H4K20me3_42h jiang2013_H4K20me3_18h
	genomedata-load-data $(GENOMEDATA_FILE) jiang2013_H3K36me2_18h < $(JIANG_2013)/SRR851666.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) jiang2013_H3K36me2_42h < $(JIANG_2013)/SRR851667.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) jiang2013_H3K36me3_18h < $(JIANG_2013)/SRR851668.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) jiang2013_H3K4me3_18h < $(JIANG_2013)/SRR851670.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) jiang2013_H3K36me3_42h < $(JIANG_2013)/SRR851669.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) jiang2013_H3K4me3_42h < $(JIANG_2013)/SRR851671.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) jiang2013_H3K9me3_18h < $(JIANG_2013)/SRR851672.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) jiang2013_H3K9me3_42h < $(JIANG_2013)/SRR851674.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) jiang2013_H4K20me3_42h < $(JIANG_2013)/SRR851676.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) jiang2013_H4K20me3_18h < $(JIANG_2013)/SRR851675.bedgraph
	genomedata-close-data $(GENOMEDATA_FILE)

load_jiang_mean_bg:
	genomedata-open-data $(GENOMEDATA_FILE) --tracknames jiang2013_mean
	genomedata-load-data $(GENOMEDATA_FILE) jiang2013_mean < $(JIANG_MEAN)/mean.bedgraph
	genomedata-close-data $(GENOMEDATA_FILE)

load_bunnik2014:
	genomedata-open-data $(GENOMEDATA_FILE) --trackname bunnik2014_MNase_lt_100_0h bunnik2014_MNase_lt_100_18h bunnik2014_MNase_lt_100_36h bunnik2014_MNase_100_to_200_0h bunnik2014_MNase_100_to_200_18h bunnik2014_MNase_100_to_200_36h bunnik2014_MNase_gt_200_0h bunnik2014_MNase_gt_200_18h bunnik2014_MNase_gt_200_36h
	genomedata-load-data $(GENOMEDATA_FILE) bunnik2014_MNase_lt_100_0h < $(BUNNIK_2014)/Nucleosome_0h_lt_100.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bunnik2014_MNase_lt_100_18h < $(BUNNIK_2014)/Nucleosome_18h_lt_100.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bunnik2014_MNase_lt_100_36h < $(BUNNIK_2014)/Nucleosome_36h_lt_100.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bunnik2014_MNase_100_to_200_0h < $(BUNNIK_2014)/Nucleosome_0h_100_to_200.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bunnik2014_MNase_100_to_200_18h < $(BUNNIK_2014)/Nucleosome_18h_100_to_200.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bunnik2014_MNase_100_to_200_36h < $(BUNNIK_2014)/Nucleosome_36h_100_to_200.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bunnik2014_MNase_gt_200_0h < $(BUNNIK_2014)/Nucleosome_0h_gt_200.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bunnik2014_MNase_gt_200_18h < $(BUNNIK_2014)/Nucleosome_18h_gt_200.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bunnik2014_MNase_gt_200_36h < $(BUNNIK_2014)/Nucleosome_36h_gt_200.bedgraph
	genomedata-close-data $(GENOMEDATA_FILE)

load_bartfai2010:
	genomedata-open-data $(GENOMEDATA_FILE) --tracknames bartfai2010_H2A.Z_10hpi bartfai2010_H2A.Z_20hpi bartfai2010_H2A.Z_30hpi bartfai2010_H2A.Z_40hpi bartfai2010_H2A_40hpi bartfai2010_H3K9ac_10hpi bartfai2010_H3K9ac_20hpi bartfai2010_H3K9ac_30hpi bartfai2010_H3K9ac_40hpi bartfai2010_H3K4me3_10hpi bartfai2010_H3K4me3_20hpi bartfai2010_H3K4me3_30hpi bartfai2010_H3K4me3_40hpi bartfai2010_input_10hpi bartfai2010_input_20hpi bartfai2010_input_30hpi bartfai2010_input_40hpi
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_H2A.Z_10hpi < $(BARTFAI_2010)/ChIPseq_H2A.Z_10hpi_log2ratio.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_H2A.Z_20hpi < $(BARTFAI_2010)/ChIPseq_H2A.Z_20hpi_log2ratio.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_H2A.Z_30hpi < $(BARTFAI_2010)/ChIPseq_H2A.Z_30hpi_log2ratio.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_H2A.Z_40hpi < $(BARTFAI_2010)/ChIPseq_H2A.Z_40hpi_log2ratio.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_H3K9ac_10hpi < $(BARTFAI_2010)/ChIPseq_H3K9ac_10hpi_log2ratio.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_H3K9ac_20hpi < $(BARTFAI_2010)/ChIPseq_H3K9ac_20hpi_log2ratio.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_H3K9ac_30hpi < $(BARTFAI_2010)/ChIPseq_H3K9ac_30hpi_log2ratio.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_H3K9ac_40hpi < $(BARTFAI_2010)/ChIPseq_H3K9ac_40hpi_log2ratio.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_H3K4me3_10hpi < $(BARTFAI_2010)/ChIPseq_H3K4me3_10hpi_log2ratio.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_H3K4me3_20hpi < $(BARTFAI_2010)/ChIPseq_H3K4me3_20hpi_log2ratio.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_H3K4me3_30hpi < $(BARTFAI_2010)/ChIPseq_H3K4me3_30hpi_log2ratio.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_H3K4me3_40hpi < $(BARTFAI_2010)/ChIPseq_H3K4me3_40hpi_log2ratio.bedgraph
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_input_10hpi < $(BARTFAI_2010)/SRR065660.bedgraph	
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_input_20hpi < $(BARTFAI_2010)/SRR065661.bedgraph	
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_input_30hpi < $(BARTFAI_2010)/SRR065662_3.bedgraph	
	genomedata-load-data $(GENOMEDATA_FILE) bartfai2010_input_40hpi < $(BARTFAI_2010)/SRR065664.bedgraph	
	genomedata-close-data $(GENOMEDATA_FILE)

load_GC:
	genomedata-open-data $(GENOMEDATA_FILE) --tracknames GC_window100
	genomedata-load-data $(GENOMEDATA_FILE) GC_window100 < $(GC)/all_GC.bedgraph
	genomedata-close-data $(GENOMEDATA_FILE)


earlyring_all:
	genomedata-load-seq $(DATADIR)/$@ $(GENOME)
	genomedata-open-data $(DATADIR)/$@ --tracknames bunnik2014_MNase_lt_100_0h bunnik2014_MNase_100_to_200_0h bunnik2014_MNase_gt_200_0h bartfai2010_H2A.Z_10hpi bartfai2010_H3K9ac_10hpi bartfai2010_H3K4me3_10hpi bartfai2010_input_10hpi GC_window100
	genomedata-load-data $(DATADIR)/$@ bunnik2014_MNase_lt_100_0h < $(BUNNIK_2014)/Nucleosome_0h_lt_100.bedgraph
	genomedata-load-data $(DATADIR)/$@ bunnik2014_MNase_100_to_200_0h < $(BUNNIK_2014)/Nucleosome_0h_100_to_200.bedgraph
	genomedata-load-data $(DATADIR)/$@ bunnik2014_MNase_gt_200_0h < $(BUNNIK_2014)/Nucleosome_0h_gt_200.bedgraph
	genomedata-load-data $(DATADIR)/$@ bartfai2010_H2A.Z_10hpi < $(BARTFAI_2010)/ChIPseq_H2A.Z_10hpi_log2ratio.bedgraph
	genomedata-load-data $(DATADIR)/$@ bartfai2010_H3K9ac_10hpi < $(BARTFAI_2010)/ChIPseq_H3K9ac_10hpi_log2ratio.bedgraph
	genomedata-load-data $(DATADIR)/$@ bartfai2010_H3K4me3_10hpi < $(BARTFAI_2010)/ChIPseq_H3K4me3_10hpi_log2ratio.bedgraph
	genomedata-load-data $(DATADIR)/$@ bartfai2010_input_10hpi < $(BARTFAI_2010)/SRR065660.bedgraph
	genomedata-load-data $(DATADIR)/$@ GC_window100 < $(GC)/all_GC.bedgraph
	genomedata-close-data $(DATADIR)/$@

troph_all:
	genomedata-load-seq $(DATADIR)/$@ $(GENOME)
	genomedata-open-data $(DATADIR)/$@ --tracknames jiang2013_H3K36me2_18h jiang2013_H3K36me3_18h jiang2013_H3K4me3_18h jiang2013_H3K9me3_18h jiang2013_H4K20me3_18h bunnik2014_MNase_lt_100_18h bunnik2014_MNase_100_to_200_18h bunnik2014_MNase_gt_200_18h bartfai2010_H2A.Z_20hpi bartfai2010_H3K9ac_20hpi bartfai2010_H3K4me3_20hpi bartfai2010_input_20hpi GC_window100
	genomedata-load-data $(DATADIR)/$@ jiang2013_H3K36me2_18h < $(JIANG_2013)/SRR851666.bedgraph
	genomedata-load-data $(DATADIR)/$@ jiang2013_H3K36me3_18h <  $(JIANG_2013)/SRR851668.bedgraph
	genomedata-load-data $(DATADIR)/$@ jiang2013_H3K4me3_18h < $(JIANG_2013)/SRR851670.bedgraph
	genomedata-load-data $(DATADIR)/$@ jiang2013_H3K9me3_18h < $(JIANG_2013)/SRR851672.bedgraph
	genomedata-load-data $(DATADIR)/$@ jiang2013_H4K20me3_18h < $(JIANG_2013)/SRR851675.bedgraph
	genomedata-load-data $(DATADIR)/$@ bunnik2014_MNase_lt_100_18h < $(BUNNIK_2014)/Nucleosome_18h_lt_100.bedgraph
	genomedata-load-data $(DATADIR)/$@ bunnik2014_MNase_100_to_200_18h < $(BUNNIK_2014)/Nucleosome_18h_100_to_200.bedgraph
	genomedata-load-data $(DATADIR)/$@ bunnik2014_MNase_gt_200_18h < $(BUNNIK_2014)/Nucleosome_18h_gt_200.bedgraph
	genomedata-load-data $(DATADIR)/$@ bartfai2010_H2A.Z_20hpi < $(BARTFAI_2010)/ChIPseq_H2A.Z_20hpi_log2ratio.bedgraph
	genomedata-load-data $(DATADIR)/$@ bartfai2010_H3K4me3_20hpi < $(BARTFAI_2010)/ChIPseq_H3K4me3_20hpi_log2ratio.bedgraph
	genomedata-load-data $(DATADIR)/$@ bartfai2010_input_20hpi < $(BARTFAI_2010)/SRR065661.bedgraph
	genomedata-load-data $(DATADIR)/$@ GC_window100 < $(GC)/all_GC.bedgraph
	genomedata-close-data $(DATADIR)/$@


schizont_all:
	genomedata-load-seq $(DATADIR)/$@ $(GENOME)
	genomedata-open-data $(DATADIR)/$@ --tracknames bartfai2010_input_40hpi bartfai2010_H2A.Z_40hpi bartfai2010_H2A_40hpi bartfai2010_H3K9ac_40hpi bartfai2010_H3K4me3_40hpi bunnik2014_MNase_lt_100_36h bunnik2014_MNase_100_to_200_36h bunnik2014_MNase_gt_200_36h jiang2013_H3K36me2_42h jiang2013_H3K36me3_42h jiang2013_H3K4me3_42h jiang2013_H3K9me3_42h jiang2013_H4K20me3_42h GC_window100
	genomedata-load-data $(DATADIR)/$@ bartfai2010_input_40hpi < $(BARTFAI_2010)/SRR065664.bedgraph
	genomedata-load-data $(DATADIR)/$@ bartfai2010_H2A.Z_40hpi < $(BARTFAI_2010)/ChIPseq_H2A.Z_40hpi_log2ratio.bedgraph
	genomedata-load-data $(DATADIR)/$@ bartfai2010_H3K9ac_40hpi < $(BARTFAI_2010)/ChIPseq_H3K9ac_40hpi_log2ratio.bedgraph
	genomedata-load-data $(DATADIR)/$@ bartfai2010_H3K4me3_40hpi < $(BARTFAI_2010)/ChIPseq_H3K4me3_40hpi_log2ratio.bedgraph
	genomedata-load-data $(DATADIR)/$@ bunnik2014_MNase_lt_100_36h < $(BUNNIK_2014)/Nucleosome_36h_lt_100.bedgraph
	genomedata-load-data $(DATADIR)/$@ bunnik2014_MNase_100_to_200_36h < $(BUNNIK_2014)/Nucleosome_36h_100_to_200.bedgraph
	genomedata-load-data $(DATADIR)/$@ bunnik2014_MNase_gt_200_36h < $(BUNNIK_2014)/Nucleosome_36h_gt_200.bedgraph
	genomedata-load-data $(DATADIR)/$@ jiang2013_H3K36me2_42h < $(JIANG_2013)/SRR851667.bedgraph
	genomedata-load-data $(DATADIR)/$@ jiang2013_H3K36me3_42h < $(JIANG_2013)/SRR851669.bedgraph
	genomedata-load-data $(DATADIR)/$@ jiang2013_H3K4me3_42h < $(JIANG_2013)/SRR851671.bedgraph
	genomedata-load-data $(DATADIR)/$@ jiang2013_H3K9me3_42h < $(JIANG_2013)/SRR851674.bedgraph
	genomedata-load-data $(DATADIR)/$@ jiang2013_H4K20me3_42h < $(JIANG_2013)/SRR851676.bedgraph
	genomedata-load-data $(DATADIR)/$@ GC_window100 < $(GC)/all_GC.bedgraph
	genomedata-close-data $(DATADIR)/$@





















