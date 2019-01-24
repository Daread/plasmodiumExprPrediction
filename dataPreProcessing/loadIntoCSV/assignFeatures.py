import pandas as pd
import pdb
import numpy as np
import itertools
from genomedata import Genome
import sys

# Classes
class featureAssinger(object):

	def __init__(self, trackNameDF, featureToAdd, winBackIn, winForIn, valueToUse = "Not Specified"):
		self.trackNames = trackNameDF
		self.thisFeature = featureToAdd
		self.windowSizeBack = winBackIn
		self.windowSizeForward = winForIn
		self.defaultVal = 0 # Default value for no signal in an existing track
		self.missingDataVal = np.nan # Value if no track available
		self.summaryStat = valueToUse

	# Make this callable to use in apply
	def __call__(self, row):
		# Get the TSS and stage
		thisTSS = row[tssOrStartCodon]
		thisStage = row["Stage"]
		thisChrom = row["ChromInfo"]
		thisOrientation = row["Strand"]

		# Track to use
		thisTrack = self.trackNames.get_value(thisStage, self.thisFeature)

		# # Debug the NAs for LT and LS
		# if((thisStage == "LT") and (self.thisFeature == "Nucleosome Occupancy")):
		# 	pdb.set_trace()
		#pdb.set_trace()
		# Make sure this is real
		if (thisTrack == "NOTRACK"):
			# No data here
			return (self.missingDataVal)

		# Get the track value
		thisChrom = myGenome[thisChrom]

		# Check strand
		if (thisOrientation == "+"):
			# Default. Use -windowback thru +windowforward
			trackVals = thisChrom[(thisTSS - self.windowSizeBack):(thisTSS + self.windowSizeForward), 
						thisTrack]
			#pdb.set_trace()
		# Reverse
		else:
			trackVals = thisChrom[(thisTSS - self.windowSizeForward):(thisTSS + self.windowSizeBack), 
			thisTrack][::-1]


		if (self.summaryStat == "Average"):
			trackAverage = np.nanmean(trackVals)

			# if (self.thisFeature == "Nucleosome Occupancy"):
			# 	print("Nucleosome")
			# 	pdb.set_trace()
			if (np.isnan(trackAverage)):
				return(self.defaultVal)
			else:
				return (trackAverage)
		# Max?
		elif (self.summaryStat == "Max"):
			trackMax = max(trackVals)
			if (np.isnan(trackMax)):
				#pdb.set_trace()
				return(self.defaultVal)
			else:
				return (trackMax)
		else:
			print("Error: No summary statistic specified")
			pdb.set_trace()
# End featureAssigner

class landmarkAssigner(object): 

	def __init__(self, ringHiCfileIn, trophHiCfileIn, schizHiCfileIn, columnToAssignIn):
		self.ringFile = ringHiCfileIn
		self.trophFile = trophHiCfileIn
		self.schizFile = schizHiCfileIn
		self.colToAssign = columnToAssignIn
		# Now get data frames
		self.ringDF = pd.read_csv(self.ringFile, sep='\t')
		self.trophDF = pd.read_csv(self.trophFile, sep='\t')
		self.schizDF = pd.read_csv(self.schizFile, sep='\t')
		# Make into dictionaries for easier look-up
		IDlist = self.ringDF["new_ID"].tolist()
		ringList = self.ringDF[self.colToAssign].tolist()
		ringMin = min(ringList)
		ringMax = max(ringList)
		trophList = self.trophDF[self.colToAssign].tolist()
		trophMin = min(trophList)
		trophMax = max(trophList)
		schizList = self.schizDF[self.colToAssign].tolist()
		schizMin = min(schizList)
		schizMax = max(schizList)
		# Assign
		self.stageSpecHiCDict = {}
		self.stageSpecHiCDict["ER"] = {}
		self.stageSpecHiCDict["LR"] = {}
		self.stageSpecHiCDict["LT"] = {}
		self.stageSpecHiCDict["LS"] = {}

		for eachInd, eachID in enumerate(IDlist):
			self.stageSpecHiCDict["ER"][eachID] = (ringList[eachInd] - ringMin) / (ringMax - ringMin)
			self.stageSpecHiCDict["LR"][eachID] = (ringList[eachInd] - ringMin) / (ringMax - ringMin)
			self.stageSpecHiCDict["LT"][eachID] = (trophList[eachInd] - trophMin) / (trophMax - trophMin)
			self.stageSpecHiCDict["LS"][eachID] = (schizList[eachInd] - schizMin) / (schizMax - schizMin)

	# Set up how each row will be assigned
	def __call__(self, row):
		# Get the TSS and stage
		thisTSS = row[tssOrStartCodon]
		thisStage = row["Stage"]
		thisGene = row["geneID"]
		# Assign based on the stage-specific dictionary assinging positions
		try:
			if (self.stageSpecHiCDict[thisStage][thisGene]):
				return (self.stageSpecHiCDict[thisStage][thisGene])
		except:
			return (np.nan)
			print("Error finding HiC data for gene " + str(thisGene))
			pdb.set_trace()

# Main

print("Starting Now")

# Set the size of windows to use
# windowBack = 1000
# windowFor = 1000

# Dictionary for having multiple windows per feature
windowDict = {}
#windowDict["TSS_Centered"] = [1000, 1000]

# Promoter vs. internal gene features
windowDict["Promoter"] = [500, 0]
windowDict["Gene_Body"] = [0, 1000]

# Alternate windows for maximums
maxDict = {}
maxDict["Promoter"] = [500, 0]
maxDict["Gene_Body"] = [0, 500]

# Define start point to use.
if (len(sys.argv) > 1):
	tssOrStartCodon = str(sys.argv[1])
	# Make sure this is valid
	if ((tssOrStartCodon != "TSS") and (tssOrStartCodon != "Start_Codon")):
		sys.exit("Error: Invalid Start Input (Must be TSS or Start_Codon")
else:
	tssOrStartCodon =  "Start_Codon" #"TSS" # "Start_Codon"

# First, read in the csv file containing stage/expression pairs
if (len(sys.argv) > 2):
	stageGenePairFile = str(sys.argv[2])
else:
	stageGenePairFile = "../2018_03_26_morePreProcessing/dividedGenFiles/"

#stageGenePairName = "chroms1_3_13ER_LR_LT_LSstageCountPairTopAndBottom3tiles.csv"
nameList = [
			"chroms1_3_13;LR;stageCountPairTopAndBottom3tiles.csv",
			"chroms1_3_13;LS;stageCountPairTopAndBottom3tiles.csv",
			"chroms1_3_13;LT;stageCountPairTopAndBottom3tiles.csv",

			"chroms7_14;LR;stageCountPairTopAndBottom3tiles.csv",
			"chroms7_14;LS;stageCountPairTopAndBottom3tiles.csv",
			"chroms7_14;LT;stageCountPairTopAndBottom3tiles.csv",

			"chroms2_9_11;LR;stageCountPairTopAndBottom3tiles.csv",
			"chroms2_9_11;LS;stageCountPairTopAndBottom3tiles.csv",
			"chroms2_9_11;LT;stageCountPairTopAndBottom3tiles.csv",

			"chroms4_5_12;LR;stageCountPairTopAndBottom3tiles.csv",
			"chroms4_5_12;LS;stageCountPairTopAndBottom3tiles.csv",
			"chroms4_5_12;LT;stageCountPairTopAndBottom3tiles.csv",

			"chroms6_8_10;LR;stageCountPairTopAndBottom3tiles.csv",
			"chroms6_8_10;LS;stageCountPairTopAndBottom3tiles.csv",
			"chroms6_8_10;LT;stageCountPairTopAndBottom3tiles.csv"


		# "chroms7_14;LR;stageCountPairTopAndBottom3tiles.csv",

					
		# 			"chroms2_9_11;LR;stageCountPairTopAndBottom3tiles.csv",

		# 	"chroms7_14;EvenByStageLR_LS;stageCountPairTopAndBottom3tiles.csv",
			
		# 	"chroms2_9_11;EvenByStageLR_LS;stageCountPairTopAndBottom3tiles.csv",

		# 	"chroms7_14;EvenByStageLR_LT_LS;stageCountPairTopAndBottom3tiles.csv",
			
		# 	"chroms2_9_11;EvenByStageLR_LT_LS;stageCountPairTopAndBottom3tiles.csv",

		# 	"chroms1_3_13;EvenByStageLR_LT_LS;stageCountPairTopAndBottom3tiles.csv",
		# 	"chroms1_3_13;EvenByStageLR_LS;stageCountPairTopAndBottom3tiles.csv",

		# 	# "chroms4_5_12;EvenByStageLR_LT_LS;stageCountPairTopAndBottom3tiles.csv",
		# 	# "chroms4_5_12;EvenByStageLR_LS;stageCountPairTopAndBottom3tiles.csv",
		# 	# "chroms4_5_12;LR;stageCountPairTopAndBottom3tiles.csv",
		# 	# "chroms4_5_12;LS;stageCountPairTopAndBottom3tiles.csv",
		# 	# "chroms4_5_12;LT;stageCountPairTopAndBottom3tiles.csv",
		# 	# "chroms4_5_12;LR_LS;stageCountPairTopAndBottom3tiles.csv",
		# 	# "chroms4_5_12;LR_LT_LS;stageCountPairTopAndBottom3tiles.csv",

		# 	"chroms7_14;LS;stageCountPairTopAndBottom3tiles.csv",
		# 	"chroms7_14;LT;stageCountPairTopAndBottom3tiles.csv",
		# 	"chroms7_14;LR_LS;stageCountPairTopAndBottom3tiles.csv",
		# 	"chroms7_14;LR_LT_LS;stageCountPairTopAndBottom3tiles.csv",
			
		# 	"chroms2_9_11;LS;stageCountPairTopAndBottom3tiles.csv",
		# 	"chroms2_9_11;LT;stageCountPairTopAndBottom3tiles.csv",
		# 	"chroms2_9_11;LR_LS;stageCountPairTopAndBottom3tiles.csv",
		# 	"chroms2_9_11;LR_LT_LS;stageCountPairTopAndBottom3tiles.csv",
			

		# 	"chroms1_3_13;LR;stageCountPairTopAndBottom3tiles.csv",
		# 	"chroms1_3_13;LS;stageCountPairTopAndBottom3tiles.csv",
		# 	"chroms1_3_13;LT;stageCountPairTopAndBottom3tiles.csv",
		# 	"chroms1_3_13;LR_LS;stageCountPairTopAndBottom3tiles.csv",
		# 	"chroms1_3_13;LR_LT_LS;stageCountPairTopAndBottom3tiles.csv"



			]

# Loop for each name
for stageGenePairName in nameList:
	print("Working on " + stageGenePairName)
	#stageGenePairFile += stageGenePairName
	groSeqData = pd.read_csv(stageGenePairFile + stageGenePairName)
	#pdb.set_trace()

	# Next, I'll open a file containing my pairings of stage-feature file. For example, 
	#	I'll use bunnik2014_MNase_100_to_200_0h for defining 
	#	nucleosome occupancy in ER/gene count pairs
	# stageToFileList = "logScaleStageToFeatureFile.csv"
	stageToFileList = "motifsIncludedFeatureToTrack.csv"
	stageFileMatches = pd.read_csv(stageToFileList)

	# Get a list of the stages
	stageList = stageFileMatches["Stage"].tolist()
	# Set the indices to be stage names
	stageFileMatches = stageFileMatches.set_index("Stage")
	# Get the list of features
	featureList = list(stageFileMatches.columns.values)
	# Since "Stage" is the index, don't need to drop it specifically.

	stageFileDict = {}
	for eachStage in stageList:
		stageFileDict[eachStage] = {}
		# Now each stage is a key to an inner dictionary. Fill this inner dictionary
		for eachFeature in featureList:
			stageFileDict[eachStage][eachFeature] = stageFileMatches.get_value(eachStage, eachFeature)


	# Next, I need to read in a genomedata file to acccess tracks for defining our features
	genomeDataDir = "/net/noble/vol2/home/katecook/proj/2016predictExpression/data/pfal3D7.genomedata"
	with Genome(genomeDataDir) as myGenome:
		# Add a functor for applying operations to each row
		for eachFeature in featureList:
			print("Assinging " + str(eachFeature))
			# fileLookupFunctor = featureAssinger(stageFileMatches, eachFeature, windowBack, windowFor)
			# # Add the row
			# groSeqData[eachFeature] = groSeqData.apply(fileLookupFunctor, axis = 1)
			# Loop for each section where averages should be taken
			for eachSection in windowDict:
				fileLookupFunctorAverage = featureAssinger(stageFileMatches, eachFeature,
					 windowDict[eachSection][0], windowDict[eachSection][1], valueToUse="Average")
				fileLookupFunctorMax = featureAssinger(stageFileMatches, eachFeature,
					 maxDict[eachSection][0], maxDict[eachSection][1], valueToUse="Max")
				# Add the row
				groSeqData["Mean " + eachFeature + ": " + eachSection] = groSeqData.apply(
					fileLookupFunctorAverage, axis = 1)
				groSeqData["Max " + eachFeature + ": " + eachSection] = groSeqData.apply(
					fileLookupFunctorMax, axis = 1)

	# 4-13-18 add
	# Add look-up for assinging distance to genomic landmarks
	ringHiCfile = "./ringDistToLandmark.txt"
	trophHiCfile = "./trophDistToLandmark.txt"
	schizHiCfile = "./schizDistToLandmark.txt"
	hiCFeaturesToAdd = ["distToTelomeres", "distToCentromeres", "distToCenter"]
	for eachFeature in hiCFeaturesToAdd:
		landmarkFunctor = landmarkAssigner(ringHiCfile, trophHiCfile, schizHiCfile,
										eachFeature)
		# Apply to make a new row
		groSeqData[eachFeature] = groSeqData.apply(landmarkFunctor, axis = 1)

	#pdb.set_trace()

	# Now output the modified file
	outputDir = "./assignedCSVfiles/"
	# Same name as before. Just has feature data added
	#windowLabel = str(float(windowBack) / 1000) + "kb" + str(float(windowFor) / 1000) + "kbFeatWindows"
	windowLabel = "_".join([x for x in windowDict])
	outputFile = outputDir + tssOrStartCodon + str(windowDict["Gene_Body"][1] / 1000) + "kbgenic" + windowLabel + stageGenePairName
	print("Outputting to " + outputFile)
	with open(outputFile, 'w') as outFile:
		groSeqData.to_csv(outFile)

#pdb.set_trace()
print("All Done")