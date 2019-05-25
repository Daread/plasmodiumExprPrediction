import pandas as pd
import pdb
import numpy as np
import itertools

# First, I will read the csv file holding all genes
print("Starting Now")
geneCountFile = "./LuEtAlFullyNormalizedCounts.csv"
myGeneDF = pd.read_csv(geneCountFile)

# Get the rows with no NA values
myGeneDF = myGeneDF.dropna(axis=0,thresh=3)

# Now get the gff file
gffFile = "../../data/2018_03_21_Dividing_P_Falciparum_Genome/PlasmoDB-29_Pfalciparum3D7.gff"
# No header, only need the 1st col (chrom) 3rd col (gene, exon, transcript) 
#     and 9th (ID and descript)
myGFFinput = pd.read_csv(gffFile, sep='\t', header=None, usecols=[0,2,8])
myGFFinput.columns = ["ChromInfo", "Type", "Description"]
# Keep only genes
#myGFFinput = myGFFinput.loc[myGFFinput["Type"] == "gene"]

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

# Discard mitochondrial genes
myGFFinput = myGFFinput.loc[myGFFinput["Chrom"] != "M76611"]

# Now get the counts for each chrom
chromCount = myGFFinput.groupby('Chrom').count()
chromCountList = chromCount["ChromInfo"].tolist()
print("Chromosome gene counts are:")
print(chromCountList)

# Store these for dict lookup
chromCountDict = {}
for eachInd in range(14):
	chromCountDict[eachInd + 1] = chromCountList[eachInd]

# Now time to run the permutations to get close to even. First, calc true even
subGroupsToForm = 5
mostEvenDist = float(sum(chromCountList)) / subGroupsToForm
print("Most even distribution is " + str(subGroupsToForm) + " groups of " + str(mostEvenDist))

groupAssignments = [[],[],[],[],[]] # Fill with chrom values
groupTotals = [0,0,0,0,0]
pickCombList = [[0], [0,1], [0,1,2], [0,1,2,3],
				[0,1,2,3,4], [0,1,2,3,4], [0,1,2,3,4], [0,1,2,3,4],
				[0,1,2,3,4], [0,1,2,3,4], [0,1,2,3,4], [0,1,2,3,4],
				[0,1,2,3,4], [0,1,2,3,4]  ]

pdb.set_trace()
# Need to track best MSE
bestMSEerror = float("inf")
comboCount = 0
# Loop and assign
for eachCombo in itertools.product(*pickCombList):
	# This will make a tuple for eachCombo for each possible allocation of chroms
	for eachIndex in range(len(eachCombo)):
		# Add to the total list
		groupTotals[eachCombo[eachIndex]] += chromCountDict[eachIndex + 1]
	# Get the MSE
	thisMSE = sum([(mostEvenDist - x)**2 for x in groupTotals])
	# Better?
	if (thisMSE < bestMSEerror):
		# yes
		bestMSEerror = thisMSE
		# Save this combo
		bestCombination = eachCombo
	# Reset the totals
	groupTotals = [0,0,0,0,0]

	# Track loop count, print updates
	comboCount += 1
	if (comboCount % 100000 == 0):
		print ("Combo " + str(comboCount))


pdb.set_trace()

print("Best combination is ")
print(bestCombination)
print("MSE: " + str(bestMSEerror))

# What are the non-gene hits here?
# nonGeneHits = myGFFinput.loc[myGFFinput["Type"] != "gene"]
# print("Getting nongenic hits")
# pdb.set_trace()

print("All Done")

