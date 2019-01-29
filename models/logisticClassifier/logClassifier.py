import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff() #http://matplotlib.org/faq/usage_faq.html (interactive mode)
import pandas as pd
import pdb
import numpy as np
import math
from sklearn import linear_model
from sklearn import metrics
from sklearn.grid_search import GridSearchCV
import random
import copy
import pickle
import sys
myRandomSeed = 77
random.seed(myRandomSeed)

# Returns: the inputDf, with all columns standardized to mean 0, sd = 1
def getZdf(inputDF, rowsToSkip = []):
	#Combine
	comboDF = copy.deepcopy(inputDF)

	# Z normalize
	for eachIndex, eachCol in enumerate(comboDF.columns):
		# Skip this?
		if (comboDF.columns.tolist()[eachIndex] in rowsToSkip):
			continue

		# Get the array
		colAsArray = np.array(comboDF[eachCol].tolist())
		# Normalize and reset. Normalize to mean 0, var = 1
		comboDF[eachCol] = (colAsArray - colAsArray.mean()) / math.sqrt(colAsArray.var())

	# Return
	return (comboDF)

#
# Functions
# Requires: dataframeList is a list of files holdingg saved pandas dataframes
#			featureColNames are columns in all df's in the dataframeList, and
#			are to be used for maxing a matrix to return
#			responseColName is a valid column in all dataframes
# Modifies: Nothing
# Effects: Returns a matrix of row-bound values from all df's, matching all columns in
#			featureColNames.
# 			Also returns a numpy array holding all values of responseColName, in order that
# 			the df's were listed in dataframeList
def getFeatureMatrixAndLabels(dataframeList, featureColNames, responseColName):
	# Need to loop through for eacch dataframe
	for eachFilename in dataframeList:
		eachDF = pd.read_csv(eachFilename)
		# Get the file set
		chromosomeSet = eachFilename.split(";")[0].split("chroms")[1]
		eachDF['chromSet'] = chromosomeSet
		# Get rid of any NA-containing rows?
		#eachDF = eachDF.dropna(axis=0, how="any")
		#pdb.set_trace()

		# Get the columns of interest
		try:
			thisFeatureDF = eachDF[featureColNames + [responseColName]]
		except:
			pdb.set_trace()
		# Get rid of NAs
		thisFeatureDF = thisFeatureDF.dropna(axis=0, how="any")

		# Get the labels from the retained rows
		featureLabels = thisFeatureDF[responseColName].tolist()
		# Set to only the real features
		thisFeatureDF = thisFeatureDF[featureColNames]
		# Try to get the labels and matrix, otherwise make the combo matrix
		try:
			combinedDF = combinedDF.append(thisFeatureDF)
			responseList += featureLabels
		except:
			# First to add
			combinedDF = thisFeatureDF
			responseList = featureLabels
	#pdb.set_trace()

	# Check if there are any NA's. Fill them anyway, but if there is a full column, raise
	#	a warning/pause
	if ( len(combinedDF.dropna(axis=1, how="all").columns) != 
		len(combinedDF.columns)):
		print("Warning: All NA column is being dropped")
		pdb.set_trace()

	# Fill missing values with zero (may want to revisit and use mean)
	combinedDF = combinedDF.fillna(0)

	# return the vlaues
	#combinedMatrix = combinedDF.as_matrix()
	return combinedDF, responseList

class LogisticClassifier(object):

	def __init__(self):
		self.AUROCdf = pd.DataFrame()
		self.trainAUROCdf = pd.DataFrame()
		self.classifAccuracyDF = pd.DataFrame()
		self.coefDF = pd.DataFrame()
		self.cvDF = pd.DataFrame()
		self.testProbsDict = {}
		#self.testProbabilsDF = pd.DataFrame()#; self.testProbabilsDF.columns = ["LR", "LT", "LS"]
		self.usedInternalCV = True
		self.classificationScores = []
		self.classificationLabels = []

	# Takes in two matrices (training and test features) and returns two matrices, corresponding
	#      to a training and test matrix (respectively) normalized so that all values are
	#      column-normalized (aka feature-normalized) z scores, normalized with a combined
	#	   training and test matrix
	def getZmatrices(self, trainDF, testDF):
		# Get rows
		trainRows = trainDF.shape[0]
		testRows = testDF.shape[0]
		#Combine
		comboDF = trainDF.append(testDF)

		# Z normalize
		for eachCol in comboDF.columns:
			# Get the array
			colAsArray = np.array(comboDF[eachCol].tolist())
			# Normalize and reset. Normalize to mean 0, var = 1
			comboDF[eachCol] = (colAsArray - colAsArray.mean()) / math.sqrt(colAsArray.var())

		# Break up again for use as separate data in training/testing
		normalizedTrainDF = comboDF[:trainRows]
		normalizedTestDF =  comboDF[trainRows:]

		#pdb.set_trace()

		# Return the values as matrices for use in the classifier
		return (normalizedTrainDF.as_matrix(), normalizedTestDF.as_matrix())


	def runLogisticModel(self, featuresToUse, crossValFileSet, runName, 
						 featuresLabel = "FullFeatures", makeOutput=True, useZscores = True,
						 useCrossValForHyper = False,
						  skipChromSets = []):

		plt.figure()
		# Get lists to track the test error and training error in each N-fold cross-validation run
		testErrorCVList = []
		trainErrorCVList = []
		testAUCcrossValList = []
		trainAUCcrossValList = []
		# alphaToUse = 0.01# .01
		# l1ToUse = 0.7
		itersToUse = 10000 #00
		# Parameters for highRegularization run
		# alphaToUse = .1
		# l1ToUse = .8

		# Each loop through, select one file to be held out
		for testFile in crossValFileSet:
			thisChromSet = testFile.split("Gene_Body")[1].split("stage")[0]
			chromOnlyList = thisChromSet.split("chroms")[1].split(";")[0]
			# Make sure this isn't one to skip
			if (chromOnlyList in skipChromSets):
				continue

			# Get the training files. These are just all the files that aren't the test file
			trainFiles = [x for x in crossValFileSet if (x != testFile)]
			trainFiles = [stagePairDir + x for x in trainFiles]
			#trainChroms = []
			testChroms = testFile.split("stage")[0]

			# Get the training and test data formatted
			trainDF, trainLabels = getFeatureMatrixAndLabels(trainFiles, featuresToUse,
										columnToPredict)
			trainMatrix = trainDF.as_matrix()
			

			testDF, testLabels = getFeatureMatrixAndLabels([stagePairDir + testFile], 
										featuresToUse, columnToPredict)
			testMatrix = testDF.as_matrix()
			testRows = trainDF.shape[0]

			# Z normalize to help in coefficient interpretation
			if (useZscores):
				trainMatrix, testMatrix = self.getZmatrices(trainDF, testDF)

			if (useCrossValForHyper):
				# Use internal CV on the training data to pick hyperparameters
				print("Beginning Internal CV")
				myClassifier = linear_model.SGDClassifier(loss='log', penalty="elasticnet",
							#alpha=0.0001, l1_ratio=0.15, 
							random_state=myRandomSeed,
							max_iter = itersToUse)

				# Cross-val to get optimal parameters
				paramSearchGrid = {"l1_ratio":[ .4, .5, .6, .7, .8, .9, .95], 
									"alpha":[1e-1, 1e-2, 1e-3, 1e-4]}
				cvSearch = GridSearchCV(myClassifier, paramSearchGrid)
				#pdb.set_trace()
				cvSearch.fit(trainMatrix, trainLabels)
				# Get the best scores
				print("Working on files:")
				print(crossValFileSet)
				print("Best parameters:")
				print(cvSearch.best_params_)
				# Set alpha and l1
				alphaToUse = cvSearch.best_params_["alpha"]
				self.cvDF.loc["alpha", thisChromSet] = alphaToUse
				l1ToUse =    cvSearch.best_params_["l1_ratio"]
				self.cvDF.loc["l1_ratio", thisChromSet] = l1ToUse
				print(" ")
			# Use from dictionary
			else:
				print("Using hard-coded values for hyper-parameters")
				self.usedInternalCV = False
				# Get for a specific stage
				thisStage = testFile.split(";")[1]
				# Get for the correct CV parameters (from CV with motifs or without)
				if ("noMotifs" in runName):
					if (thisStage == "LR"):
						alphaToUse = .0001
						l1ToUse = .95
					elif (thisStage == "LT"):
						alphaToUse = .0001
						l1ToUse = .8
					elif (thisStage == "LS"):
						alphaToUse = .0001
						l1ToUse = .9	
					else:
						print("Error on stage selection")
						sys.exit(1)  
				elif ("withMotifs" in runName):
					if (thisStage == "LR"):
						alphaToUse = .01
						l1ToUse = .95
					elif (thisStage == "LT"):
						alphaToUse = .0001
						l1ToUse = .9
					elif (thisStage == "LS"):
						alphaToUse = .001
						l1ToUse = .95
					else:
						print("Error on stage selection")
						sys.exit(1) 
				else:
					print("Needs a motif usage designation to hard-set hyperparameters")
					sys.exit(1)


			# Fit a model
			print("Fitting Model")
			myClassifier = linear_model.SGDClassifier(loss='log', penalty= "elasticnet",
						alpha = alphaToUse, l1_ratio = l1ToUse, 
						random_state=myRandomSeed,
						max_iter = itersToUse)
			#print("Classifying Now")
			# Train on the training data
			myClassifier.fit(trainMatrix, trainLabels)

			# Get the dataframe for retrieving column names
			featureDF = pd.read_csv(stagePairDir + testFile)[featuresToUse]

			# Match up the coefficients with the feature labels
			coefDict = {}
			for eachIndex, eachCoeff in enumerate(myClassifier.coef_[0]):
				coefDict[featureDF.columns[eachIndex]] = eachCoeff
			# Output these coefficients
			print("Coefficients from model:")
			
			for eachFeature in coefDict:
				self.coefDF.loc[eachFeature, thisChromSet] = coefDict[eachFeature]
				#print(eachFeature + ": " + str(coefDict[eachFeature]))

			# Now get the training and test predictions, to assess accuracy
			trainPredictions = myClassifier.predict(trainMatrix)
			trainMatches = [np.array([True for eachInd in range(len(trainPredictions)) if 
							trainPredictions[eachInd] == trainLabels[eachInd]])]
			trainPredCorrect = (1.0 * np.sum(trainMatches)) / len(trainPredictions)
			#justScoreCorrect = myClassifier.score(trainMatrix, trainLabels)

			# Test as well
			testPredictions = myClassifier.predict(testMatrix)
			#pdb.set_trace()
			testMatches = [np.array([True for eachInd in range(len(testPredictions)) if 
							testPredictions[eachInd] == testLabels[eachInd]])]
			testPredCorrect = (1.0 * np.sum(testMatches)) / len(testPredictions)

			# Add these to the lists to average for a final CV score
			testErrorCVList.append(testPredCorrect)
			trainErrorCVList.append(trainPredCorrect)

			# Get the training AUROC
			#testScores = myClassifier.decision_function(testMatrix)
			trainScores = myClassifier.decision_function(trainMatrix)
			#pdb.set_trace()
			falsePos, truePos, _ = metrics.roc_curve(y_true = trainLabels, 
						y_score=trainScores, pos_label="Low")
			# Get the curve plotted
			thisChromSet = testFile.split("Gene_Body")[1].split("stage")[0]
			binarytrainLabels = [1 if (x == "Low") else 0 for x in trainLabels]
			trainAUC = metrics.roc_auc_score(y_true = binarytrainLabels,
							y_score=trainScores)
			trainAUCcrossValList.append(trainAUC)

			# Re-scale the labels as 0/1 for use in the auc_score function
			testScores = myClassifier.decision_function(testMatrix)
			falsePos, truePos, _ = metrics.roc_curve(y_true = testLabels, 
						y_score=testScores, pos_label="Low")
			# Get the curve plotted
			thisChromSet = testFile.split("Gene_Body")[1].split("stage")[0]
			binaryTestLabels = [1 if (x == "Low") else 0 for x in testLabels]
			thisAUC = metrics.roc_auc_score(y_true = binaryTestLabels,
							y_score=testScores)
			testAUCcrossValList.append(thisAUC)
			# Add the score to the main DF for output later
			self.AUROCdf.loc[featuresLabel, thisChromSet] = thisAUC
			self.trainAUROCdf.loc[featuresLabel, thisChromSet] = trainAUC
			self.classifAccuracyDF.loc[featuresLabel, thisChromSet] = testPredCorrect

			# Get AUCs for individual genes
			if (thisStage in self.testProbsDict):
				tempDF = pd.DataFrame()
				tempDF["Probabilities"] = testScores
				tempDF["Labels"] = binaryTestLabels
				# Add to previous df
				self.testProbsDict[thisStage] = self.testProbsDict[thisStage].append(tempDF, 
							ignore_index = True)
			else:
				# Make new DF
				self.testProbsDict[thisStage] = pd.DataFrame()
				self.testProbsDict[thisStage]["Probabilities"] = testScores
				self.testProbsDict[thisStage]["Labels"] = binaryTestLabels


			# 8-22-18
			print(" ")
			print("At breakpoint")
			#pdb.set_trace()
			# Add scores:
			if (len(self.classificationScores) == 0):
				self.classificationScores = testScores
			else:
				self.classificationScores = np.concatenate((self.classificationScores, testScores),
															axis=None )
			# Labels:
			self.classificationLabels = self.classificationLabels + testLabels

			# Add a ROC curve for this x-fold validation sets
			if (makeOutput):
				plt.plot(falsePos, truePos, label=thisChromSet)
				# Save true and false pos
				print("Writing pickled output")
				#pdb.set_trace()
				pickleOut = "./pickledOutput/" + runName + "TruPos" + testFile.split(";")[1] + testFile.split(";")[0].split("chroms")[1]
				with open(pickleOut, 'wb') as outFile:
					pickle.dump(truePos, outFile)
				pickleOut = "./pickledOutput/" + runName + "FalsePos" + testFile.split(";")[1] + testFile.split(";")[0].split("chroms")[1]
				with open(pickleOut, 'wb') as outFile:
					pickle.dump(falsePos, outFile)
			
		# Output the total CV results
		print(" ")
		print("Fitting a model for:")
		print(crossValFileSet)
		print("Accuracy: ")
		print("Test accuracy:")
		print(testErrorCVList)
		testErrorAverage = sum(testErrorCVList) * 1.0 / len(testErrorCVList)
		print("Average: " + str(testErrorAverage))
		print("AUC:")
		print(testAUCcrossValList)
		print("Average: " + str(sum(testAUCcrossValList) * 1.0 / len(testAUCcrossValList) ) )
		print(" ")
		print("train accuracy:")
		print(trainErrorCVList)
		trainErrorAverage = sum(trainErrorCVList) * 1.0 / len(trainErrorCVList)
		print("Average: " + str(trainErrorAverage))
		print("AUC:")
		print(trainAUCcrossValList)
		print(" ")
		print(" " )

		# Calculate a combined ROC across the folds
		#pdb.set_trace()
		combinedFalsePos, combinedTruePos, _ = metrics.roc_curve(y_true = self.classificationLabels, 
						y_score=self.classificationScores, pos_label="Low")
		# Add a ROC curve for this x-fold validation sets
		#pdb.set_trace()
		if (makeOutput):
			plt.plot(falsePos, truePos, label=thisChromSet)
			# Save true and false pos
			print("Writing pickled combined ROC curve")
			#pdb.set_trace()
			pickleOut = "./pickledOutput/" + runName + "TruPos" + testFile.split(";")[1] + foldSet + motifUse
			with open(pickleOut, 'wb') as outFile:
				pickle.dump(combinedTruePos, outFile)
			pickleOut = "./pickledOutput/" + runName + "FalsePos" + testFile.split(";")[1] + foldSet + motifUse
			with open(pickleOut, 'wb') as outFile:
				pickle.dump(combinedFalsePos, outFile)
		

		# Output the ROC curve
		if (makeOutput):
			outputLabel = testFile.split("chroms")[1].split("stage")[0]
			outputLabel = "".join([x for x in outputLabel if not x.isdigit()])
			plt.title("ROC for " + str(len(crossValFileSet)) + "-fold Test Error, " + outputLabel)
			plt.xlabel("False Positive Rate")
			plt.ylabel("True Positive Rate")
			leg=plt.legend(bbox_to_anchor=(1.105, 1), loc=2, borderaxespad=0., title="Test Set")
			plt.savefig("./ROCoutputs/" + runName + outputLabel + featuresLabel + "_ROC.pdf", 
				#bbox_extra_artists=(leg), 
				bbox_inches='tight')
		# Always clear
		plt.clf()
	# End runLogisticModel


#######################################################################3
#Main
print("Starting Now")

# Usage:
# python logClassifier.py <file set, threeDevFolds or fullFiveFolds> <motifs, noMotifs or withMotifs>
#		


# Start by reading in the assigned value dataframe
if (len(sys.argv) > 3):
	stagePairDir = str(sys.argv[3])
else:
	stagePairDir = "../2018_04_30_feature_revisit/assignedCSVfiles/"
# Fold given?
if (len(sys.argv) > 1):
	foldSet = str(sys.argv[1])
else:
	foldSet = "threeDevFolds" # fullFiveFolds

useTSS = False

#
if (foldSet == "fullFiveFolds"):
	fileListsToUse =    [
			  
			["startCodon1kbgenicPromoter_Gene_Bodychroms6_8_10;LR;stageCountPairTopAndBottom3tiles.csv",
			"startCodon1kbgenicPromoter_Gene_Bodychroms4_5_12;LR;stageCountPairTopAndBottom3tiles.csv",

			"startCodon1kbGenicPromoter_Gene_Bodychroms1_3_13;LR;stageCountPairTopAndBottom3tiles.csv",
			"startCodon1kbGenicPromoter_Gene_Bodychroms2_9_11;LR;stageCountPairTopAndBottom3tiles.csv",
			"startCodon1kbGenicPromoter_Gene_Bodychroms7_14;LR;stageCountPairTopAndBottom3tiles.csv"],

			["startCodon1kbgenicPromoter_Gene_Bodychroms6_8_10;LT;stageCountPairTopAndBottom3tiles.csv",
			"startCodon1kbgenicPromoter_Gene_Bodychroms4_5_12;LT;stageCountPairTopAndBottom3tiles.csv",

			"startCodon1kbGenicPromoter_Gene_Bodychroms1_3_13;LT;stageCountPairTopAndBottom3tiles.csv",
			"startCodon1kbGenicPromoter_Gene_Bodychroms2_9_11;LT;stageCountPairTopAndBottom3tiles.csv",
			"startCodon1kbGenicPromoter_Gene_Bodychroms7_14;LT;stageCountPairTopAndBottom3tiles.csv"],

			["startCodon1kbgenicPromoter_Gene_Bodychroms6_8_10;LS;stageCountPairTopAndBottom3tiles.csv",
			"startCodon1kbgenicPromoter_Gene_Bodychroms4_5_12;LS;stageCountPairTopAndBottom3tiles.csv",

			"startCodon1kbGenicPromoter_Gene_Bodychroms1_3_13;LS;stageCountPairTopAndBottom3tiles.csv",
			"startCodon1kbGenicPromoter_Gene_Bodychroms2_9_11;LS;stageCountPairTopAndBottom3tiles.csv",
			"startCodon1kbGenicPromoter_Gene_Bodychroms7_14;LS;stageCountPairTopAndBottom3tiles.csv"]

			]
	setsToSkipThisFold = ["1_3_13", "2_9_11", "7_14"]
	cvForHyper = False

# 3 fold:
elif (foldSet == "threeDevFolds"):
	fileListsToUse = [
			["startCodon1kbGenicPromoter_Gene_Bodychroms1_3_13;LR;stageCountPairTopAndBottom3tiles.csv",
			"startCodon1kbGenicPromoter_Gene_Bodychroms2_9_11;LR;stageCountPairTopAndBottom3tiles.csv",
			"startCodon1kbGenicPromoter_Gene_Bodychroms7_14;LR;stageCountPairTopAndBottom3tiles.csv"],

			["startCodon1kbGenicPromoter_Gene_Bodychroms1_3_13;LT;stageCountPairTopAndBottom3tiles.csv",
			"startCodon1kbGenicPromoter_Gene_Bodychroms2_9_11;LT;stageCountPairTopAndBottom3tiles.csv",
			"startCodon1kbGenicPromoter_Gene_Bodychroms7_14;LT;stageCountPairTopAndBottom3tiles.csv"],

			["startCodon1kbGenicPromoter_Gene_Bodychroms1_3_13;LS;stageCountPairTopAndBottom3tiles.csv",
			"startCodon1kbGenicPromoter_Gene_Bodychroms2_9_11;LS;stageCountPairTopAndBottom3tiles.csv",
			"startCodon1kbGenicPromoter_Gene_Bodychroms7_14;LS;stageCountPairTopAndBottom3tiles.csv"]
			]

	# Use TSS instead?
	if (useTSS):
		print("Using TSS, not start codon")
		pdb.set_trace()
		fileListsToUse = [
					["TSS1kbgenicPromoter_Gene_Bodychroms1_3_13;LR;stageCountPairTopAndBottom3tiles.csv",
			"TSS1kbgenicPromoter_Gene_Bodychroms2_9_11;LR;stageCountPairTopAndBottom3tiles.csv",
			"TSS1kbgenicPromoter_Gene_Bodychroms7_14;LR;stageCountPairTopAndBottom3tiles.csv"],

			["TSS1kbgenicPromoter_Gene_Bodychroms1_3_13;LT;stageCountPairTopAndBottom3tiles.csv",
			"TSS1kbgenicPromoter_Gene_Bodychroms2_9_11;LT;stageCountPairTopAndBottom3tiles.csv",
			"TSS1kbgenicPromoter_Gene_Bodychroms7_14;LT;stageCountPairTopAndBottom3tiles.csv"],

			["TSS1kbgenicPromoter_Gene_Bodychroms1_3_13;LS;stageCountPairTopAndBottom3tiles.csv",
			"TSS1kbgenicPromoter_Gene_Bodychroms2_9_11;LS;stageCountPairTopAndBottom3tiles.csv",
			"TSS1kbgenicPromoter_Gene_Bodychroms7_14;LS;stageCountPairTopAndBottom3tiles.csv"]
			]


				
	setsToSkipThisFold = []
	cvForHyper = True
else:
	print("Invalid fold set")
	sys.exit(1)


columnToPredict = "Expression Level"

# Motifs given?
if (len(sys.argv) > 2):
	motifUse = sys.argv[2]
else:
	#motifUse = "noMotifs"
	motifUse = "withMotifs"

if (cvForHyper):
	cvName = "internalCV"
else:
	cvName = "hardSetHyperparams"

thisRun = "logistic" + motifUse + foldSet + cvName
if (useTSS):
	thisRun = "TSSbased" + thisRun
trainClassifier = True
outputFullDFs = False
outputStandardizedDFs = False

if (motifUse == "withMotifs"):
	fullFeatureFile = "./averageAndMaxFeatureList.txt"
	limitedFeatueFile = "./threeStageAverageAndMaxFeatureList.txt"
elif (motifUse == "noMotifs"):
	fullFeatureFile = "./noMotifAverageAndMax.txt"
	limitedFeatueFile = "./threeStageNoMotifAverageAndMax.txt"
elif (motifUse == "motifOnly"):
	fullFeatureFile = "./motifOnlyList.txt"
	limitedFeatueFile = "./motifOnlyList.txt"
elif (motifUse == "motifAndGC"):
	fullFeatureFile = "./motifAndGC.txt"
	limitedFeatueFile = "./motifAndGC.txt"
elif (motifUse == "GConly"):
	fullFeatureFile = "./GConlyList.txt"
	limitedFeatueFile = "./GConlyList.txt"	
else:
	# Into 
	print("Invalid motif selection")
	sys.exit(1)

print("Selections:")
print("Fold set: " + foldSet)
print("Motif usage: " + motifUse)

#stagePairFile = "TSS_Centered_Promoter_Gene_Bodychroms1_3_13LRstageCountPairTopAndBottom3tiles.csv"
myLogModel = LogisticClassifier()


for crossValFileSet in fileListsToUse:
	# Get the appropriate features
	featuresToUse = []
	thisModel = crossValFileSet[0].split(";")[1]
	# Full or limited feature list? Full list = available at only 2 timpepoints
	if ("LT" in crossValFileSet[0]):
		featureInputFile = limitedFeatueFile
	else:
		featureInputFile = fullFeatureFile
	# Open, read in features
	with open(featureInputFile, 'r') as inFile:
		for eachLine in inFile:
			# Append this to the list
			featuresToUse.append(eachLine.strip())
			#print(eachLine.strip())
		# Store the master list
		fullFeatureList = featuresToUse

	if (trainClassifier):
		myLogModel.runLogisticModel(fullFeatureList, crossValFileSet, thisRun,
								 makeOutput=True, useCrossValForHyper = cvForHyper,
								  skipChromSets = setsToSkipThisFold)

	# Get the full data output
	if (outputFullDFs):
		# Also get the raw dataframe from all data
		crossValFileSet = [stagePairDir + x for x in crossValFileSet]
		fullDataFrame, fullLabels = getFeatureMatrixAndLabels(crossValFileSet, 
										 ["geneID", "Stage",  "chromSet"] + 
										 fullFeatureList, 
										columnToPredict)

		fullDataFrame["Expression Label"] = fullLabels
		fullDataFrame.reset_index(inplace=True)
		set6_8_10 = fullDataFrame.index[fullDataFrame["chromSet"] == "6_8_10"].tolist()
		set4_5_12 = fullDataFrame.index[fullDataFrame["chromSet"] == "4_5_12"].tolist()
		set1_3_13 = fullDataFrame.index[fullDataFrame["chromSet"] == "1_3_13"].tolist()
		set2_9_11 = fullDataFrame.index[fullDataFrame["chromSet"] == "2_9_11"].tolist()
		set7_14 = fullDataFrame.index[fullDataFrame["chromSet"] == "7_14"].tolist()

		#pdb.set_trace()
		# Get the start and end of these
		#fullDataFrame.drop(["chromSet"], axis=1, inplace=True)
		print(" ")
		print(motifUse)
		print("6_8_10 from " + str(set6_8_10[0]) + " to " + str(set6_8_10[-1]))
		print("4_5_12 from " + str(set4_5_12[0]) + " to " + str(set4_5_12[-1] )) 
		print("1_3_13 from " + str(set1_3_13[0]) + " to " + str(set1_3_13[-1] )) 
		print("2_9_11 from " + str(set2_9_11[0]) + " to " + str(set2_9_11[-1] )) 
		print("7_14 from " + str(set7_14[0]) + " to " + str(set7_14[-1] )) 
		#pdb.set_trace()
		# output
		
		#pdb.set_trace()
		outputFile = "./featureAndResponseFiles/" + foldSet + "/fullData" + thisModel  + ".csv"
		with open(outputFile, 'w') as outFile:
			print("Writing " + str(outputFile))
			fullDataFrame.to_csv(outFile)

	# Get standardized DF for use in DeepPink
	if (outputStandardizedDFs):
		# Also get the raw dataframe from all data
		crossValFileSet = [stagePairDir + x for x in crossValFileSet]
		fullDataFrame, fullLabels = getFeatureMatrixAndLabels(crossValFileSet, 
										fullFeatureList + ["chromSet"], 
										columnToPredict)
		# Standardize the DF
		fullDataFrame = getZdf(fullDataFrame, rowsToSkip = ["chromSet"])

		fullDataFrame["Response"] = fullLabels
		fullDataFrame.reset_index(inplace=True)
		set6_8_10 = fullDataFrame.index[fullDataFrame["chromSet"] == "6_8_10"].tolist()
		set4_5_12 = fullDataFrame.index[fullDataFrame["chromSet"] == "4_5_12"].tolist()

		#pdb.set_trace()
		# Get the start and end of these
		fullDataFrame.drop(["chromSet"], axis=1, inplace=True)
		print(" ")
		print(motifUse)
		print("6_8_10 from " + str(set6_8_10[0]) + " to " + str(set6_8_10[-1]))
		print("4_5_12 from " + str(set4_5_12[0]) + " to " + str(set4_5_12[-1] )) 
		#pdb.set_trace()
		# output
		
		#pdb.set_trace()
		outputFile = "./featureAndResponseFiles/" + foldSet + "/standardized" + thisModel + motifUse + ".csv"
		with open(outputFile, 'w') as outFile:
			print("Writing " + str(outputFile))
			fullDataFrame.to_csv(outFile)

	
	# # Single-feature runs
	# for miniFeatureList in fullFeatureList:
	# 	# Run for this
	# 	myLogModel.runLogisticModel([miniFeatureList], crossValFileSet,
	# 	 featuresLabel = miniFeatureList, makeOutput = False)
	# #pdb.set_trace()

# Output
print("About to save outputs")

#pdb.set_trace()
if (trainClassifier):
	myLogModel.AUROCdf.to_csv("./AUROCcomparisons/" + thisRun + "AUROCvals.csv")
	myLogModel.trainAUROCdf.to_csv("./AUROCcomparisons/" + thisRun + "trainAUROCvals.csv")
	myLogModel.classifAccuracyDF.to_csv("./AUROCcomparisons/" + thisRun + "classificationVals.csv")
	myLogModel.coefDF.to_csv("./AUROCcomparisons/" + thisRun + "Coefficients.csv")
	# Save hyperparams?
	if (myLogModel.usedInternalCV):
		myLogModel.cvDF.to_csv("./hyperparams/" + thisRun + "hyperparamsByCV.csv")

	# Get the raw probabilistic scores
	#pdb.set_trace()
	for eachStage in myLogModel.testProbsDict:
		# Note: Deep Pink output clipped the 1396th entry in LT data. Clip
		#		here for consistency
		print("Writing for " + eachStage)
		myLogModel.testProbsDict[eachStage].iloc[:1396].to_csv("./probabilisticOutputs/" +
						eachStage + "probabilies" + thisRun + ".csv")
	# pdb.set_trace()
	# myLogModel.testProbabilsDF.to_csv("./probabilisticOutputs/" + 
	# 			thisRun + "classificationScores.csv")
#pdb.set_trace()

print("All Done")