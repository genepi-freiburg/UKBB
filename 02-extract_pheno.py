#!/usr/bin/python

import sys
import getopt
import argparse
import csv
import collections
from operator import itemgetter
import os
from os import listdir
import re
import glob




######################################
######## predefined variables ########
######################################

ukbbPhenoFilePath1 = "/dsk/data1/programs/scripts/UKBB/test_extract_pheno/ukb8974_new_head.csv"
ukbbPhenoFilePath2 = "/dsk/data1/programs/scripts/UKBB/test_extract_pheno/ukb28633_new_head.csv"

ukbbPhenoFilePath1 = "/data/studies/06_UKBB/UKBB_150k/01_Raw_Data/ukb8974.csv"
ukbbPhenoFilePath2 = "/data/studies/06_UKBB/Biomarkers/ukb28633.csv"


exclusionsDefaultPath = "/data/studies/06_UKBB/EXCLUSIONS/"






####################################
######## read in parameters ########
####################################

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help = "Path of input file", action = "store", required = True)
parser.add_argument("-o", "--out", help = "File where the output should be written", action = "store", required = True)
parser.add_argument("-e", "--exclude", help = "Path to exclude list. One PatId per line", action = "store", required = False)

args = parser.parse_args()

inputFilePath = args.input
outFileName = args.out
excludeListPath = args.exclude

#### check if excludeListPath is set else use default file
if not excludeListPath:
	# print (glob.glob(exclusionsDefaultPath + "*"))
	excludeListPath = exclusionsDefaultPath + "/" + sorted(listdir(exclusionsDefaultPath))[-1]
	args.exclude = excludeListPath



#### print out chosen options
print ("Chosen options")
for arg in vars(args):
    print ("\t--"+arg, getattr(args, arg))
print ("\n")

outDir = os.path.dirname(outFileName)
if not outDir:
	outDir = os.getcwd()
if (not os.path.exists(outDir)):
	os.makedirs(outDir)



#########################
######## methods ########
#########################

###########################
######## read in input file


def readInInputFile():
	print("Reading in input file:", inputFilePath)
	fileIn = open(inputFilePath, "r")
	input = fileIn.readlines()
	fileIn.close()


	ukbbVarMappingDict = dict()
	wildCardMappingDict = collections.defaultdict(list)

	phenoPathVarNameMappingDict = collections.defaultdict(list)
	for curLine in input:

		if (curLine.startswith("#") or curLine.isspace()):
			continue
		splitLine = curLine.replace("\n","").split("\t")
		ukbbVar = splitLine[0]
		outVar = splitLine[1]
		ukbbVarMappingDict[ukbbVar] = outVar

	#### check for corresponding UKBB phenotype file
	print ("Checking and mapping variants in phenotpe files")
	ukbbPhenoHeader1 = open(ukbbPhenoFilePath1, "r").readline().replace("\"", "").replace("\n","").split(",")
	ukbbPhenoHeader2 = open(ukbbPhenoFilePath2, "r").readline().replace("\"", "").replace("\n","").split(",")
	nonMappedVar = list()

	for curUkbbVar in ukbbVarMappingDict.keys():

		#### check for whildcards
		if (curUkbbVar.find("*") != -1):
			mainString = curUkbbVar.replace(".*", "")
			if (not re.search("\d+-\d+.*", mainString)):
				print("\nERROR: Not allowed wildcard search. Check your input file.")
				sys.exit(2)
			if (any(s for s in ukbbPhenoHeader1 if mainString in s)):
				wildCardMappingDict[curUkbbVar].append(list(s for s in ukbbPhenoHeader1 if mainString in s))
				phenoPathVarNameMappingDict[ukbbPhenoFilePath1].append(curUkbbVar)
			elif (any(s for s in ukbbPhenoHeader2 if mainString in s)):
				wildCardMappingDict[curUkbbVar].append(s for s in ukbbPhenoHeader1 if mainString in s)
				phenoPathVarNameMappingDict[ukbbPhenoFilePath2].append(curUkbbVar)
			else:
				nonMappedVar.append(curUkbbVar)

		#### sort normal variables
		else:
			if curUkbbVar in ukbbPhenoHeader1:
				phenoPathVarNameMappingDict[ukbbPhenoFilePath1].append(curUkbbVar)
			elif curUkbbVar in ukbbPhenoHeader2:
				phenoPathVarNameMappingDict[ukbbPhenoFilePath2].append(curUkbbVar)
			else:
				nonMappedVar.append(curUkbbVar)

	if len(nonMappedVar) > 0:
		nonMatchedOutFileName = outDir + "/non_matched_variants.txt"
		nonMatchedOutFile = open(nonMatchedOutFileName, "w")
		print ("\n".join(nonMappedVar), file = nonMatchedOutFile)
		nonMatchedOutFile.close()


	return ukbbVarMappingDict, phenoPathVarNameMappingDict, wildCardMappingDict


################################
######## read in exclude Id list

def readInExcludeList(excludeListPath):

	excludeList = list()
	if not excludeListPath is None and os.path.isfile(excludeListPath):
		print ("Reading in exclusion list:", excludeListPath)
		with open(excludeListPath, "r") as input:
			for row in input:
				excludeList.append(row.replace("\n",""))
	else:
		excludeList.append("")
	return excludeList



##################################################
######## create map of header entries with indeces

def createheaderIndexMap(row):
	headerIndexMap = dict()
	index = 0;
	for curEntry in row:
		headerIndexMap[curEntry] = index
		index += 1
	return headerIndexMap



def createHeaderIndexList (headerToKeep, headerIndexMap):
	headerIndexList = list()
	wildcardheaderIndexMap = collections.defaultdict(list)
	for curHeader in headerToKeep:
		if curHeader.find("*" == -1):
			headerIndexList.append(headerIndexMap.get(curHeader))
		else:
			wildcardheaderIndexMap

	return headerIndexList



###############################
######## extract ukbb variables
def aggregateVariables(row, headerList, headerIndexMap):
	varList = list()

	for curHeader in headerList[0]:
		varList.append(row[headerIndexMap.get(curHeader)])
	finalVarList = list (filter(None, varList))
	return ",".join(finalVarList)




def extractDataFromRow (row, headerList, headerIndexMap, wildCardMappingDict):
	extractedDataList = list()


	for curHeader in headerList:
		if curHeader.find("*") == -1:
			extractedDataList.append(row[headerIndexMap.get(curHeader)])
		else:
			extractedDataList.append(aggregateVariables(row, wildCardMappingDict.get(curHeader), headerIndexMap))

	return extractedDataList




################################################
######## read in phenotype file and extract data


def extractData(phenoPath, phenoPathVarNameMappingDict, ukbbVarMappingDict, outFileDict, wildCardMappingDict):
	print ("Reading in and extracting data from: ", phenoPath)

	with open(phenoPath) as csvFile:
		phenoFile = csv.reader(csvFile)
		isHeader = True
		headerIndexMap = dict()
		index = -1
		for row in phenoFile:

			if isHeader == True:
				isHeader = False
				headerIndexMap = createheaderIndexMap(row)

				newHeader = list()
				for curHeader in phenoPathVarNameMappingDict.get(phenoPath):
						newHeader.append(ukbbVarMappingDict.get(curHeader))
				outFileDict["Individual_ID"].append(newHeader)
				continue

			outFileDict[row[0]].append(extractDataFromRow (row, phenoPathVarNameMappingDict[phenoPath], headerIndexMap, wildCardMappingDict))

			# index +=1
			# if index == 10:
			# 	break


	return outFileDict




#################################
######## write out phenotype file

def writeOutPhenoFile(outFileDict, excludeList):
	print ("Writing out to:", outFileName)
	fileOut = open(outFileName, "w")

	for patId in outFileDict:
		if patId in excludeList:
			continue
		outLine = patId
		for curList in outFileDict[patId]:
			outLine += "\t" + "\t".join(curList)


		# if (len(outLine.split("\t")) != 54):
		# 	print (patId)
			# print (outLine.split("\t"))

		fileOut.write(outLine + "\n")



#############################
######## main Script ########
#############################


ukbbVarMappingDict, phenoPathVarNameMappingDict, wildCardMappingDict = readInInputFile()

outFileDict = collections.defaultdict(list)

excludeList = readInExcludeList(excludeListPath)

for phenoPath in phenoPathVarNameMappingDict:
	outFileDict = extractData(phenoPath, phenoPathVarNameMappingDict, ukbbVarMappingDict, outFileDict, wildCardMappingDict)

writeOutPhenoFile(outFileDict, excludeList)
