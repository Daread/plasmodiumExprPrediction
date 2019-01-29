import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff() #http://matplotlib.org/faq/usage_faq.html (interactive mode)
import pandas as pd
import pdb
import numpy as np
import itertools
import copy
from math import log10
import sys

class TSSassigner(object):

	def __init__(self, geneToTSSdictIn):
		self.geneToTSS = geneToTSSdict

	def __call__(self, row):
		# In the dictionary to use?
		if (row["Description"] in self.geneToTSS):
			return (self.geneToTSS[row["Description"]])
		# No. Use CDS coordinates if needed
		else:
			print("No TSS for " + row["Description"])
			# Old method
			if row['Strand'] == "+":
				tssVal = row["Start"]
			else:
				tssVal = row["End"]
			return tssVal


# Functions

def checkStrandedness(row):
    if row['Strand'] == "+":
        tssVal = row["Start"]
    else:
        tssVal = row["End"]
    return tssVal

# Reads in a dataframe with TSS annotations for genes to classify
# Return a list where each entry is, for a given row, the distance to the nearest
#		upstream gene, wrt the orientation of the gene on that row.
#		the upstream gene could be a TSS, or the end of a transcribed region
def getUpstreamTSS(groSeqDF):
	defaultValue = 100000
	#pdb.set_trace()
	# Get the orientations and TSS's
	TSSlist = groSeqDF["TSS"].tolist()
	orientationList = groSeqDF["Strand"].tolist()
	chromosomeList = groSeqDF["ChromInfo"].tolist()
	startList = groSeqDF["Start"].tolist()
	endList = groSeqDF["End"].tolist()
	orientationList = [1 if (x == "+") else -1 for x in orientationList]
	# Now got through and find the orientation
	lastGenelist = []
	for eachIndex, eachTSS in enumerate(TSSlist):
		thisOrient = orientationList[eachIndex]
		# Forward facing gene
		if (thisOrient == 1):
			# Make sure this isn't the first gene in the list
			if (eachIndex == 0):
				# Start is in sense orientation. No previous
				lastGenelist.append(defaultValue)
			# See if this is a new chromosome
			elif (chromosomeList[eachIndex - 1] != chromosomeList[eachIndex]):
				# Yes. Infinite distance
				lastGenelist.append(defaultValue)
			# internal. Find the latest TSS/end of transcript
			else:
				thisVal =  min(abs(eachTSS - startList[eachIndex - 1]),
				 						abs(eachTSS -  endList[eachIndex - 1]))
				#abs(eachTSS - TSSlist[eachIndex - 1])
				# min(abs(eachTSS - startList[eachIndex - 1]),
				# 						abs(eachTSS -  endList[eachIndex - 1]))
				if (thisVal == 0):
					thisVal = lastGenelist[-1]
				lastGenelist.append(thisVal)
		# End + orientation
		# So this must be anti-sense
		else:
			if (eachIndex == (len(TSSlist) - 1)):
				# Lst gene. Nothing beyond
				lastGenelist.append(defaultValue)
			elif (chromosomeList[eachIndex + 1] != chromosomeList[eachIndex]):
				# New chrom
				lastGenelist.append(defaultValue)
			else:
				thisVal = min(abs(eachTSS - startList[eachIndex + 1]),
				 						abs(eachTSS -  endList[eachIndex + 1]))
				#abs(eachTSS - TSSlist[eachIndex + 1])
				# min(abs(eachTSS - startList[eachIndex + 1]),
				# 						abs(eachTSS -  endList[eachIndex + 1]))
				if (thisVal == 0):
					thisVal = lastGenelist[-1]
				lastGenelist.append(thisVal)

	# Make an output for the distances to previous genes
	#pdb.set_trace()
	logGeneList = [log10(x) for x in lastGenelist]
	plt.hist(logGeneList)
	plt.savefig("./outputFiles/distanceToGeneHist.pdf",
		bbox_inches='tight')
	plt.clf()
	# Now make a cdf
	y,binEdges = np.histogram(logGeneList, bins=80)
	binCenters = 0.5 * (binEdges[1:] + binEdges[:-1])
	plt.plot(binCenters, np.cumsum(y), "-")
	plt.title("CDF of Distance to Nearest Gene/TSS From a TSS")
	plt.xlabel("Log10(Distance to nearest upstream gene)")
	plt.ylabel("Number of genes")
	plt.xlim(0,6)
	plt.savefig("./outputFiles/distanceToGeneCDF.pdf",
		bbox_inches='tight')
	plt.clf()

	# Return this list
	#pdb.set_trace()
	return(lastGenelist)

# Hard-coded: List of lists for chromosommes to divide up into Xtiles
# Output from divideGenomeIntoFiveSets: (0, 1, 0, 2, 2, 3, 4, 3, 1, 3, 1, 2, 0, 4)
#                                        1  2  3  4  5  6  7  8  9 10 11 12 13 14
chromBins = [[1, 3, 13], [2, 9, 11], [4, 5, 12], [6, 8, 10], [7, 14]]
xTilesToMake = 3

# Read in the GRO-seq data
geneCountFile = "./LuEtAlFullyNormalizedCounts.csv"
myGeneDF = pd.read_csv(geneCountFile)


# Repeated code, gets the data formated and into one dataframe
# First, I will read the csv file holding all genes
print("Starting Script Now")
geneCountFile = "./LuEtAlFullyNormalizedCounts.csv"
myGeneDF = pd.read_csv(geneCountFile)

# Get the rows with no NA values
myGeneDF = myGeneDF.dropna(axis=0,thresh=3)

# Now get the gff file
if (len(sys.argv) > 1):
	gffFile = str(sys.argv[1])
else:
	gffFile = "../../data/2018_03_21_Dividing_P_Falciparum_Genome/PlasmoDB-29_Pfalciparum3D7.gff"
# No header, only need the 1st col (chrom) 3rd col (gene, exon, transcript) 
#	   4th/5th/6th (start/stop/orientation)
#     and 9th (ID and descript)
myGFFinput = pd.read_csv(gffFile, sep='\t', header=None, usecols=[0,2,3,4,6,8])
myGFFinput.columns = ["ChromInfo", "Type", "Start", "End", "Strand", "Description"]

#pdb.set_trace()


# Keep ones that are in the data file as protein-coding
protCoding = list(myGeneDF.iloc[:,0])
# Format of Description column is ID=PF3D7_0621350;description=protein...
#                 for genes. need after =, before ;
myGFFinput["Description"] = myGFFinput["Description"].str.split("=").str[1] # After =
myGFFinput["Description"] = myGFFinput["Description"].str.split(";").str[0] # Before ;
myGFFinput = myGFFinput.loc[myGFFinput["Description"].isin(protCoding)]

# Now get the genes out of the geneCountFile. These are the GRO-seq entries
myGFFinput["Chrom"] = myGFFinput["ChromInfo"].str.split("_").str[1] # format is Pf3D7_CHROMNUM_v3, where
							# CHROMNUM is 13, 02, etc.


# Read in file that has assigned TSS's for all genes based on Adjalley et al datasets
# See 2018_05_10_get_TSS_data/assignMainTSS.py for details
TSSassignmentFile = "../2018_05_10_get_TSS_data/mainTSSassignments.tsv"
geneToTSSdict = {}
with open(TSSassignmentFile, 'r') as inFile:
	for eachLine in inFile:
		splitLine = eachLine.split()
		geneToTSSdict[splitLine[0]] = int(splitLine[1])
myTSSassigner = TSSassigner(geneToTSSdict)
myGFFinput["TSS"] = myGFFinput.apply(myTSSassigner, axis=1)

# Also get the start codon position
myGFFinput["Start_Codon"] = myGFFinput.apply(checkStrandedness, axis=1)

# Discard mitochondrial genes
myGFFinput = myGFFinput.loc[myGFFinput["Chrom"] != "M76611"]
myGFFinput.sort_values(inplace = True, by=["Chrom", "TSS"])
myGFFinput["Next TSS Upstream"] = getUpstreamTSS(myGFFinput)
#pdb.set_trace()

keptIDlist = myGFFinput["Description"].tolist()

# Keep only the genes matching this
myGeneDF = myGeneDF.loc[myGeneDF["geneID"].isin(keptIDlist)]

# Get the chromosome info. Convert to integer format
myGeneDF = pd.merge(myGeneDF, myGFFinput, left_on = "geneID", right_on = "Description")
chromListInts = myGeneDF["Chrom"].tolist()
chromListInts = [int(x) for x in chromListInts]
myGeneDF["Chrom"] = chromListInts
myGeneDF.drop("Description", axis = 1, inplace=True)

# Get the TSS

# Main loop. Go through, make sub-lists for each set of chromosomes. Use this to make Xtiles
subdivideDir = "./dividedGenFiles/"
for chromList in chromBins:
	thisBinDF = myGeneDF.loc[myGeneDF["Chrom"].isin(chromList)]
	# Output
	thisSet = "_".join(str(x) for x in chromList) + "fullListGROseqGenes.csv"
	with open(subdivideDir + thisSet, 'w') as outFile:
		thisBinDF.to_csv(outFile, index=False)

# Now drop the early trophozite and early schizont to get subsets with feature coverage
myGeneDF.drop("ET", axis = 1, inplace=True)
myGeneDF.drop("ES", axis = 1, inplace=True)
myGeneDF.drop("G.III", axis = 1, inplace=True)
myGeneDF.drop("G.V", axis = 1, inplace=True)
# Write the data for each chromosome group. I want an output of the full data for those chroms,
# 		as well as top and bottom Xtile gene/stage pairs (two files)
modelList = [["LR"], ["LT"], ["LS"], ["LR", "LS"], ["LR", "LT", "LS"]]
singleModelDict = {}

#pdb.set_trace()

# setsToKeep = ["LR",  "LS"]
for setsToKeep in modelList:
	for chromList in chromBins:
		thisBinDF = copy.deepcopy(myGeneDF.loc[myGeneDF["Chrom"].isin(chromList)])
		# Output the full data for this chromosome group
		thisSet = ("chroms" +  "_".join(str(x) for x in chromList) + "_" +
		 	"_".join(setsToKeep)+ "_GROseqGenes.csv")
		with open(subdivideDir + thisSet, 'w') as outFile:
			thisBinDF.to_csv(outFile, index=False)

		# Break into Xtiles
		thisBinDF.drop("Gene Description", axis = 1, inplace=True)
		
		for eachSet in setsToKeep:
			# See if this ISN'T the first
			try:
				geneStagePairDF
				tempDF = thisBinDF.copy()
				tempDF["Stage"] = eachSet
				tempDF["GROcount"] = tempDF[eachSet]
				tempDF.drop(setsToKeep, axis = 1, inplace=True)
				geneStagePairDF = geneStagePairDF.append(tempDF)
			except:
				tempDF = thisBinDF.copy()
				tempDF["Stage"] = eachSet
				tempDF["GROcount"] = tempDF[eachSet]
				geneStagePairDF = tempDF.drop(setsToKeep, axis = 1)

		#geneStagePairDF.sort(inplace = True, columns=["GROcount"]) # Old 
		geneStagePairDF.sort_values(inplace = True, by=["GROcount"])
		rowCount = geneStagePairDF.shape[0]
		# Get the index cutoffs for the Xtiles desired
		rowPerXtile = int(rowCount / xTilesToMake)
		# Write outputs for each chromosome bin
		topTileDF = copy.deepcopy(geneStagePairDF.head(rowPerXtile))
		topTileDF["Expression Level"] = "Low"
		#pdb.set_trace()
		bottomTileDF = copy.deepcopy(geneStagePairDF.tail(rowPerXtile))
		bottomTileDF["Expression Level"] = "High"
		# Get a df with just top and bottom Xtiles
		highAndLowDF = topTileDF.append(bottomTileDF)
		# Add labels to the full DF as well
		geneStagePairDF["Expression Level"] = "Intermediate"
		geneStagePairDF.iloc[:rowPerXtile, geneStagePairDF.columns.get_loc('Expression Level')] = "Low"
		geneStagePairDF.iloc[-rowPerXtile:, geneStagePairDF.columns.get_loc('Expression Level')] = "High"
		#pdb.set_trace()
		# Output
		fullPairsFile = ("chroms" + "_".join(str(x) for x in chromList) + "_".join(setsToKeep) + 
					"stageCountPairFull.csv")
		with open(subdivideDir + fullPairsFile, 'w') as outFile:
			geneStagePairDF.to_csv(outFile, index=False)
		# Output Xtiles
		topAndBottomFile = ("chroms" + "_".join(str(x) for x in chromList) + ";" + "_".join(setsToKeep) + 
					";stageCountPairTopAndBottom" + str(xTilesToMake) + "tiles.csv")
		with open(subdivideDir + topAndBottomFile, 'w') as outFile:
			highAndLowDF.to_csv(outFile, index=False)
		print("Done Writing Sets for ", ", ".join(str(x) for x in chromList))
		# Clear the DF
		geneStagePairDF = None

		# If this is a single stage model, save it
		if (len(setsToKeep) == 1):
			singleModelDict[setsToKeep[0] + str(chromList[0])] = highAndLowDF
		else:
			# First, add the first
			#pdb.set_trace()
			combinedDF = singleModelDict[setsToKeep[0] + str(chromList[0])]
			#pdb.set_trace()
			# Append the rest
			for eachInd in range(1,len(setsToKeep)):
				#pdb.set_trace()
				combinedDF = combinedDF.append(singleModelDict[setsToKeep[eachInd] + str(chromList[0])])
			# Output
			evenByStageFile = ("chroms" + "_".join(str(x) for x in chromList) + ";" + "EvenByStage" + "_".join(setsToKeep) + 
					";stageCountPairTopAndBottom" + str(xTilesToMake) + "tiles.csv")
			#pdb.set_trace()
			with open(subdivideDir + evenByStageFile, 'w') as outFile:
				combinedDF.to_csv(outFile, index=False)



print("All Done")