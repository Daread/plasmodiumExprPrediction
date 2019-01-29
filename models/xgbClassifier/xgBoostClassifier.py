import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff() #http://matplotlib.org/faq/usage_faq.html (interactive mode)
import xgboost as xgb
import pandas as pd
import pdb
import numpy as np
import math
import copy
from sklearn import linear_model
from sklearn import metrics
from sklearn.grid_search import GridSearchCV
from scipy.stats import spearmanr
import random
import shap
import pickle
import sys
myRandomSeed = 77
random.seed(myRandomSeed)

# Return the hard-set hyper-parametesr for XGboost, as found previously using CV
#    within the dev set
def getStageParamDict(runName):
	# Motifs?
	if ("withMotifs" in runName):
		print("Error: Haven't gotten parameters with motifs")
		sys.exit(1)

	# No motifs
	if ("noMotifs" in runName):
		stageParamDict = {}
		# Ring
		LRdict = {}
		LRdict["max_depth"] = 4
		LRdict["min_child_weight"] = 5
		LRdict["subsample"] = .7
		LRdict["colsample_bylevel"] = .3
		LRdict["n_estimators"] = 100
		stageParamDict["LR"] = LRdict
		# Troph
		LTdict = {}
		LTdict["max_depth"] = 6
		LTdict["min_child_weight"] = 5
		LTdict["subsample"] = .4
		LTdict["colsample_bylevel"] = .5
		LTdict["n_estimators"] = 60
		stageParamDict["LT"] = LTdict
		# Schizont
		LSdict = {}
		LSdict["max_depth"] = 4
		LSdict["min_child_weight"] = 6
		LSdict["subsample"] = .4
		LSdict["colsample_bylevel"] = .3
		LSdict["n_estimators"] = 80
		stageParamDict["LS"] = LSdict
	# Return
	return(stageParamDict)


#pdb.set_trace()

#
# Functions
# Requires: dataframeList is a list of pandas dataframes
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
		#pdb.set_trace()
		# Get rid of any NA-containing rows?
		#eachDF = eachDF.dropna(axis=0, how="any")

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

class XGClassifierClass(object):

	def __init__(self):
		self.AUROCdf = pd.DataFrame()
		self.trainAUROCdf = pd.DataFrame()
		self.classifAccuracyDF = pd.DataFrame()
		self.coefDF = pd.DataFrame()
		self.cvDF = pd.DataFrame()
		self.testProbsDict = {}

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


	def runXGBModel(self, featuresToUse, crossValFileSet, runName,
						 featuresLabel = "FullFeatures", 
						 makeROCOutput=True, 
						 useZscores = True,
						 useCrossVal = True, 
						 outputImportance = True,
						 skipSets = []):

		plt.figure()
		# Get lists to track the test error and training error in each N-fold cross-validation run
		testErrorCVList = []
		trainErrorCVList = []
		testAUCcrossValList = []
		trainAUCcrossValList = []

		# Each loop through, select one file to be held out
		for testFile in crossValFileSet:

			thisChromSet = testFile.split("Gene_Body")[1].split("stage")[0]
			chromOnlyList = thisChromSet.split("chroms")[1].split(";")[0]
			# Skip?
			#pdb.set_trace()
			if (chromOnlyList in skipSets):
				# One of the sets where we aren't interested in the AUROC val
				#    (like the dev sets w/in the 5-fold runs)
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
			#pdb.set_trace()
			if (useZscores):
				trainMatrix, testMatrix = self.getZmatrices(trainDF, testDF)

			# Fit a model
			if (useCrossVal):
				print ("Beginning Internal CV")
				myClassifier = xgb.XGBClassifier()
				#pdb.set_trace()
				#print("Classifying Now")
				# Train on the training data
				#myClassifier.fit(trainMatrix, trainLabels)

				# CV for hyperparameter tuning
				print("Beginning CV for hyperparameters")
				paramDict = {}
				paramDict["max_depth"] = range(4, 9, 1)
				paramDict["min_child_weight"] = range(2, 7, 1)
				paramDict["subsample"] = [ .3, .4, .7, .8, 1.0]
				paramDict["colsample_bylevel"] = [.3, .5, .7, 1.]
				paramDict["n_estimators"] = [40, 60, 80, 100]

				# Smaller version, testing only
				# print("Using subset of cv options, for optimization!")
				# paramDict = {}
				# paramDict["max_depth"] = range(4, 5, 1)
				# paramDict["min_child_weight"] = range(2, 3, 1)
				# paramDict["subsample"] = [ .3, .8]
				# paramDict["colsample_bylevel"] = [.3,  1]
				# paramDict["n_estimators"] = [40,  100]

				myGridSearch = GridSearchCV(myClassifier, param_grid = paramDict,
									scoring = "roc_auc")
				binarytrainLabels = [1 if (x == "Low") else 0 for x in trainLabels]
				myGridSearch.fit(trainMatrix, binarytrainLabels)
				print("Best parameters:")
				print(myGridSearch.best_params_)
				paramToUseDict = myGridSearch.best_params_
				self.cvDF.loc["max_depth", thisChromSet] = paramToUseDict["max_depth"]
				self.cvDF.loc["min_child_weight", thisChromSet] = paramToUseDict["min_child_weight"]
				self.cvDF.loc["subsample", thisChromSet] = paramToUseDict["subsample"]
				self.cvDF.loc["colsample_bylevel", thisChromSet] = paramToUseDict["colsample_bylevel"]
				self.cvDF.loc["n_estimators", thisChromSet] = paramToUseDict["n_estimators"]
				#pdb.set_trace()
			else:
				print("Not using interanl CV!")
				
				stageParamDict = getStageParamDict(runName)

				thisStage = testFile.split(';')[1]
				try:
					paramToUseDict = stageParamDict[thisStage]
				except:
					print("Error getting hard-coded hyperparameters")
					sys.exit("Error with hyperparameter setting")

			print("Beginning to fit model")
			myClassifier = xgb.XGBClassifier(max_depth = paramToUseDict["max_depth"],
						min_child_weight = paramToUseDict["min_child_weight"],
						subsample = paramToUseDict["subsample"],
						colsample_bylevel = paramToUseDict["colsample_bylevel"],
						n_estimators = paramToUseDict["n_estimators"])


			myClassifier.fit(trainMatrix, trainLabels)

			print("Model Trained")
			#pdb.set_trace()

			# Get SHAP explanations of features
			shapValues = shap.TreeExplainer(myClassifier).shap_values(trainMatrix)
			# The output contains a last column which is the explanation model output
			# 		for this data point. We don't want this for our analysis of SHAP
			#       value vs. feature value. Drop it
			fullShapValues = copy.deepcopy(shapValues)
			shapValues = shapValues[:,:-1]
			# Now we need to get the correlation direction for each feature
			#shapCorrelArray = np.empty(shape=shapValues.shape[1])
			for eachColNum in range(shapValues.shape[1]):
				# Get this column
				thisFeatShapVals = shapValues[:,eachColNum]
				thisFeatValues = trainMatrix[:,eachColNum]
				# Get the correlation
				#pdb.set_trace()
				#shapCorrelArray[eachColNum] 
				shapFeatInteractionList = [a / b for a,b in zip(thisFeatShapVals, thisFeatValues)]
				thisAverage = sum(shapFeatInteractionList) / (1.0 * len(shapFeatInteractionList))
				thisCorrel = spearmanr(thisFeatShapVals, thisFeatValues)[0] #0th index is the rho
				shapSum = sum([abs(x) for x in thisFeatShapVals]) / len(thisFeatShapVals)
				if (thisCorrel > 0.0):
					thisCorrel = 1
				else:
					thisCorrel = -1
				valToRecord = thisCorrel * shapSum
				# Add to the DF
				
				eachFeature = trainDF.columns[eachColNum]
				self.coefDF.loc[eachFeature, thisChromSet] = valToRecord

			# Also make outputs for SHAP scores on a gene-by-gene basis for some 
			# 	gene families of interst, like var genes
			#outputShapByGeneType = True
			shapValDF = pd.DataFrame(fullShapValues)
			nullMatrix, trainIDs = getFeatureMatrixAndLabels(trainFiles, 
												featuresToUse, 
													"geneID")
			shapValDF["geneID"] = trainIDs
			trainDFwithID = copy.deepcopy(trainDF)
			trainDFwithID["geneID"] = trainIDs
			# Can either do SHAP or ROC, otherwise over-writes
			# if (makeROCOutput == False):
			# 	geneFamIDDir = "./shapData/"
			# 	geneFamIDfiles = ["varGeneIDs.txt", "exportedProteinIDs.txt", "AP2proteinIDs.txt",
			# 						"fullGFFIDs.txt"]
			# 	geneFamIDfiles = [geneFamIDDir + x for x in geneFamIDfiles]
			# 	# Open and read each
			# 	for eachFile in geneFamIDfiles:
			# 		# Read in, append each
			# 		thisIDlist = []
			# 		with open(eachFile, 'r') as inFile:
			# 			for eachLine in inFile:
			# 				thisIDlist.append(eachLine.rstrip())
			# 		# Now get this sub-DF
			# 		#pdb.set_trace()
			# 		thisShapValDF = shapValDF[shapValDF["geneID"].isin(thisIDlist)]
			# 		#pdb.set_trace()
			# 		# Get rid of the IDs
			# 		thisShapValDF.drop(["geneID"], axis=1, inplace=True)
			# 		thisTrainDF = trainDFwithID[trainDFwithID["geneID"].isin(thisIDlist)]

			# 		thisTrainDF.drop(["geneID"], axis=1, inplace=True)
			# 		# Now make the plot
			# 		#pdb.set_trace()
			# 		shap.summary_plot(thisShapValDF.as_matrix(), thisTrainDF)
			# 		# Save the figure
			# 		outputPath = geneFamIDDir + eachFile.split(".txt")[0].split("/")[-1] + "SHAPvalues.pdf"
			# 		plt.savefig(outputPath)
			# 		plt.clf()


			#pdb.set_trace()

			# Get the dataframe for retrieving column names
			featureDF = pd.read_csv(stagePairDir + testFile)[featuresToUse]

			thisChromSet = testFile.split("Gene_Body")[1].split("stage")[0]

			# Now get the training and test predictions, to assess accuracy
			#pdb.set_trace()
			trainPredictions = myClassifier.predict(trainMatrix)
			trainMatches = [np.array([True for eachInd in range(len(trainPredictions)) if 
							trainPredictions[eachInd] == trainLabels[eachInd]])]
			trainPredCorrect = (1.0 * np.sum(trainMatches)) / len(trainPredictions)

			# Test as well
			#pdb.set_trace()
			testPredictions = myClassifier.predict(testMatrix)
			testMatches = [np.array([True for eachInd in range(len(testPredictions)) if 
							testPredictions[eachInd] == testLabels[eachInd]])]
			testPredCorrect = (1.0 * np.sum(testMatches)) / len(testPredictions)

			# Add these to the lists to average for a final CV score
			testErrorCVList.append(testPredCorrect)
			trainErrorCVList.append(trainPredCorrect)

			# Get the training AUROC
			trainScores = myClassifier.predict_proba(trainMatrix)[:,1]
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
			#testScores = myClassifier.decision_function(testMatrix)
			testScores = myClassifier.predict_proba(testMatrix)[:,1]
			#pdb.set_trace()
			falsePos, truePos, _ = metrics.roc_curve(y_true = testLabels, 
						y_score=testScores, pos_label="Low")
			# Get the curve plotted
			thisChromSet = testFile.split("Gene_Body")[1].split("stage")[0]
			binaryTestLabels = [1 if (x == "Low") else 0 for x in testLabels]
			thisAUC = metrics.roc_auc_score(y_true = binaryTestLabels,
							y_score=testScores)
			testAUCcrossValList.append(thisAUC)

			# Add the test scores to the main DF for output later
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

				# 		# Add a ROC curve for this x-fold validation sets
			if (makeROCOutput):
				plt.plot(falsePos, truePos, label=thisChromSet)

				# Save the models to make shap gene-by-gene plots later
				pickledOutput = "./pickledOutput/" + runName + testFile.split(";")[1] + testFile.split(";")[0].split("chroms")[1]
				with open(pickledOutput, 'wb') as outFile:
					pickle.dump(myClassifier, outFile)
				
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

		# Output the values of the importances of features
		#if (outputImportance):

			# plt.bar(range(len(featureArray)), featureArray)#, tick_label = featuresToUse)
			# theseStages = crossValFileSet[0].split(";")[1]
			# plt.title("Feature Importances for " + theseStages)
			# plt.xticks(range(len(featuresToUse)), featuresToUse, rotation='vertical')
			# thisFile = "./featureImportances/" + theseStages + "XGboost"
			# plt.savefig(thisFile + ".pdf", bbox_inches = "tight")

		# Output the ROC curve
		if (makeROCOutput):
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


	# End 


#######################################################################3
#Main
print("Starting Now")

# Start by reading in the assigned value dataframe

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

thisRun = "xgboost" + motifUse + foldSet + cvName

#thisRun = "xgboostInternalCV" + motifUse + foldSet
trainClassifier = True
outputFullDFs = False

if (motifUse == "withMotifs"):
	fullFeatureFile = "../2018_04_19_logistic_classifier/averageAndMaxFeatureList.txt"
	limitedFeatueFile = "../2018_04_19_logistic_classifier/threeStageAverageAndMaxFeatureList.txt"
elif (motifUse == "noMotifs"):
	fullFeatureFile = "../2018_04_19_logistic_classifier/noMotifAverageAndMax.txt"
	limitedFeatueFile = "../2018_04_19_logistic_classifier/threeStageNoMotifAverageAndMax.txt"
else:
	print("Invalid motif selection")
	sys.exit(1)

print("Selections:")
print("Fold set: " + foldSet)
print("Motif usage: " + motifUse)

#stagePairFile = "TSS_Centered_Promoter_Gene_Bodychroms1_3_13LRstageCountPairTopAndBottom3tiles.csv"
myXGBmodel = XGClassifierClass()
for crossValFileSet in fileListsToUse:
	# Get the appropriate features
	featuresToUse = []
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

	myXGBmodel.runXGBModel(fullFeatureList, crossValFileSet, thisRun,
							useCrossVal = cvForHyper, skipSets = setsToSkipThisFold) #, makeOutput=True)

	# # Single-feature runs
	# for miniFeatureList in fullFeatureList:
	# 	# Run for this
	# 	myXGBmodel.runXGBModel([miniFeatureList], crossValFileSet,
	# 	 featuresLabel = miniFeatureList, makeOutput = False)
	# #pdb.set_trace()

# Output
pdb.set_trace()

myXGBmodel.AUROCdf.to_csv("./AUROCcomparisons/" + thisRun + "AUROCvals.csv")
myXGBmodel.trainAUROCdf.to_csv("./AUROCcomparisons/" + thisRun + "trainAUROCvals.csv")
myXGBmodel.classifAccuracyDF.to_csv("./AUROCcomparisons/" + thisRun + "lassificationVals.csv")
myXGBmodel.coefDF.to_csv("./AUROCcomparisons/" + thisRun + "Coefficients.csv")
myXGBmodel.cvDF.to_csv("./hyperparams/" + thisRun + "hyperparamsByCV.csv")

# Get the raw probabilistic scores
#pdb.set_trace()
for eachStage in myXGBmodel.testProbsDict:
	# Note: Deep Pink output clipped the 1396th entry in LT data. Clip
	#		here for consistency
	print("Writing for " + eachStage)
	myXGBmodel.testProbsDict[eachStage].iloc[:1396].to_csv("./probabilisticOutputs/" +
					eachStage + "probabilies" + thisRun + ".csv")

pdb.set_trace()
	# Now individual features


print("All Done")