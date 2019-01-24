DATADIR = /net/noble/vol2/home/katecook/proj/2016predictExpression/data
GENOMEDATA_FILE = $(DATADIR)/pfal3D7.genomedata

GENOME = /net/noble/vol1/data/reference_genomes/pfal3D7-PlasmoDBv29/PlasmoDB-29_Pfalciparum3D7_Genome.fasta

motifDir = /net/gs/vol1/home/readdf/plasmodiumExprPred/results/katecook/20170424_scan_with_pfms

all: load_campbell2010 load_CISBP load_CISBP_RNA_DNAencoded

# Load each directory worth of bedgraphs
load_campbell2010:
	#genomedata-open-data $(GENOMEDATA_FILE) --tracknames campbell2010_PF07_0126_DLD campbell2010_PF10_0075_D1 campbell2010_PF10_0075_D2 campbell2010_PF10_0075_D3 campbell2010_PF11_0091 campbell2010_PF11_0404_D1 campbell2010_PF11_0442 campbell2010_PF13_0026 campbell2010_PF13_0097 campbell2010_PF13_0235_D1 campbell2010_PF13_0267 campbell2010_PF14_0079 campbell2010_PF14_0533 campbell2010_PFD0985w_D1 campbell2010_PFD0985w_D2 campbell2010_PFE0840c_D2 campbell2010_PFF0670w_D1 campbell2010_PFF0670w_D2 campbell2010_PFL1075w campbell2010_PFL1085w campbell2010_PFL1900w_D1 campbell2010_PFL1900w_DLD
#	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_MAL8P1 < $(motifDir)/campbell2010/campbell2010_MAL8P1.bedgraph
# RE-ADD campbell2010_MAL8P1  to tracknames when re-adding!
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PF07_0126_DLD < $(motifDir)/campbell2010/campbell2010_PF07_0126_DLD.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PF10_0075_D1 < $(motifDir)/campbell2010/campbell2010_PF10_0075_D1.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PF10_0075_D2 < $(motifDir)/campbell2010/campbell2010_PF10_0075_D2.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PF10_0075_D3 < $(motifDir)/campbell2010/campbell2010_PF10_0075_D3.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PF11_0091 < $(motifDir)/campbell2010/campbell2010_PF11_0091.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PF11_0404_D1 < $(motifDir)/campbell2010/campbell2010_PF11_0404_D1.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PF11_0442 < $(motifDir)/campbell2010/campbell2010_PF11_0442.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PF13_0026 < $(motifDir)/campbell2010/campbell2010_PF13_0026.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PF13_0097 < $(motifDir)/campbell2010/campbell2010_PF13_0097.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PF13_0235_D1 < $(motifDir)/campbell2010/campbell2010_PF13_0235_D1.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PF13_0267 < $(motifDir)/campbell2010/campbell2010_PF13_0267.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PF14_0079 < $(motifDir)/campbell2010/campbell2010_PF14_0079.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PF14_0533 < $(motifDir)/campbell2010/campbell2010_PF14_0533.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PFD0985w_D1 < $(motifDir)/campbell2010/campbell2010_PFD0985w_D1.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PFD0985w_D2 < $(motifDir)/campbell2010/campbell2010_PFD0985w_D2.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PFE0840c_D2 < $(motifDir)/campbell2010/campbell2010_PFE0840c_D2.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PFF0670w_D1 < $(motifDir)/campbell2010/campbell2010_PFF0670w_D1.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PFF0670w_D2 < $(motifDir)/campbell2010/campbell2010_PFF0670w_D2.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PFL1075w < $(motifDir)/campbell2010/campbell2010_PFL1075w.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PFL1085w < $(motifDir)/campbell2010/campbell2010_PFL1085w.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PFL1900w_D1 < $(motifDir)/campbell2010/campbell2010_PFL1900w_D1.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) campbell2010_PFL1900w_DLD < $(motifDir)/campbell2010/campbell2010_PFL1900w_DLD.bedgraph	
	genomedata-close-data $(GENOMEDATA_FILE)
         
load_CISBP:
	genomedata-open-data $(GENOMEDATA_FILE) --tracknames CIS-BP_M0046_1.02 CIS-BP_M0047_1.02 CIS-BP_M0048_1.02 CIS-BP_M0051_1.02 CIS-BP_M0052_1.02 CIS-BP_M0054_1.02 CIS-BP_M0055_1.02 CIS-BP_M0056_1.02 CIS-BP_M0057_1.02 CIS-BP_M0058_1.02 CIS-BP_M0059_1.02 CIS-BP_M0060_1.02 CIS-BP_M0061_1.02 CIS-BP_M0063_1.02 CIS-BP_M0064_1.02 CIS-BP_M0066_1.02 CIS-BP_M0067_1.02 CIS-BP_M0069_1.02 CIS-BP_M0071_1.02
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0046_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0046_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0047_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0047_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0048_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0048_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0051_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0051_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0052_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0052_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0054_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0054_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0055_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0055_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0056_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0056_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0057_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0057_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0058_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0058_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0059_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0059_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0060_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0060_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0061_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0061_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0063_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0063_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0064_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0064_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0066_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0066_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0067_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0067_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0069_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0069_1.02.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CIS-BP_M0071_1.02 < $(motifDir)/CIS-BP/CIS-BP_M0071_1.02.bedgraph
	genomedata-close-data $(GENOMEDATA_FILE)

load_CISBP_RNA_DNAencoded:
	genomedata-open-data $(GENOMEDATA_FILE) --tracknames CISBP-RNA_DNAencoded_M179_0.6 CISBP-RNA_DNAencoded_M180_0.6 CISBP-RNA_DNAencoded_M181_0.6 CISBP-RNA_DNAencoded_M182_0.6 CISBP-RNA_DNAencoded_M183_0.6 CISBP-RNA_DNAencoded_M202_0.6 CISBP-RNA_DNAencoded_M203_0.6 CISBP-RNA_DNAencoded_M208_0.6
	genomedata-load-data -v $(GENOMEDATA_FILE) CISBP-RNA_DNAencoded_M179_0.6 < $(motifDir)/CISBP-RNA_DNAencoded/CISBP-RNA_DNAencoded_M179_0.6.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CISBP-RNA_DNAencoded_M180_0.6 < $(motifDir)/CISBP-RNA_DNAencoded/CISBP-RNA_DNAencoded_M180_0.6.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CISBP-RNA_DNAencoded_M181_0.6 < $(motifDir)/CISBP-RNA_DNAencoded/CISBP-RNA_DNAencoded_M181_0.6.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CISBP-RNA_DNAencoded_M182_0.6 < $(motifDir)/CISBP-RNA_DNAencoded/CISBP-RNA_DNAencoded_M182_0.6.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CISBP-RNA_DNAencoded_M183_0.6 < $(motifDir)/CISBP-RNA_DNAencoded/CISBP-RNA_DNAencoded_M183_0.6.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CISBP-RNA_DNAencoded_M202_0.6 < $(motifDir)/CISBP-RNA_DNAencoded/CISBP-RNA_DNAencoded_M202_0.6.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CISBP-RNA_DNAencoded_M203_0.6 < $(motifDir)/CISBP-RNA_DNAencoded/CISBP-RNA_DNAencoded_M203_0.6.bedgraph
	genomedata-load-data -v $(GENOMEDATA_FILE) CISBP-RNA_DNAencoded_M208_0.6 < $(motifDir)/CISBP-RNA_DNAencoded/CISBP-RNA_DNAencoded_M208_0.6.bedgraph
	genomedata-close-data $(GENOMEDATA_FILE)


