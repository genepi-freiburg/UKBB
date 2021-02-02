#!/usr/bin/python

import sys
import getopt
import argparse
import glob
import collections
import os
from os import listdir
import subprocess
import shutil
from pathlib import Path

#################################################
######## fix variables, pathes and tools ########
#################################################

#### program pathes
bgenixPath = "/data/programs/bin/gwas/bgenix/bgenix"
qcToolPath = "/data/programs/bin/gwas/qctool/qctool_v2.0.1/qctool_v2.0.1"


#### imputed data
#impPath = "/data/studies/06_UKBB/UKBB_500k_V3/01_Raw_Data/EGAD00010001474/"
impPath = "/data/studies/06_UKBB/01_Data/01_Raw_Data/Appl_20272_Genetic_2018/EGAD00010001474/"
#impSampleFile = "/data/studies/06_UKBB/UKBB_500k_V3/01_Raw_Data/UKBB_500k_v3_clean.sample"
impSampleFile = "/data/studies/06_UKBB/01_Data/01_Raw_Data/Appl_20272_Genetic_2018/UKBB_500k_v3_clean.sample"

#### genotyped data
#genoPath = "/data/studies/06_UKBB/UKBB_500k_V3/01_Raw_Data/EGAD00010001497/"
genoPath = "/data/studies/06_UKBB/01_Data/01_Raw_Data/Appl_20272_Genetic_2018/EGAD00010001497/"
#genoFamPath = "/data/studies/06_UKBB/UKBB_500k/01_Raw_Data/genetic/ready_PLINK/ukb2027_cal_chr1_v2_s488366.fam"
genoFamPath = "/data/studies/06_UKBB/01_Data/01_Raw_Data/Appl_20272_Genetic_2017/genetic/ready_PLINK/ukb2027_cal_chr1_v2_s488366.fam"
plinkPath = "/data/programs/bin/gwas/plink/plink-2.0.0_alpha_20190724/plink2"




#### WES data
#wesDataPath = "/data/studies/06_UKBB/Exome_50k/01_Raw_Data/EFE/"
# wesDataPath = "/data/studies/06_UKBB/Exome_200k/02_Genotypes/BGEN_Format/output/"
wesDataPath = "/data/studies/06_UKBB/01_Data/01_Raw_Data/Appl_20272_Exome_50k/EFE/"
wesBedFile = wesDataPath + "ukb_efe_chr1_v1.bed"
wesBimFile = wesDataPath + "ukb_fe_exm_chrall_v1.bim"
wesFamFile = wesDataPath + "ukb20272_efe_chr1_v1_s49959.fam"

#wesPath = "/data/studies/06_UKBB/Exome_200k/02_Genotypes/BGEN_Format/output"
wesPath = "/data/studies/06_UKBB/01_Data/02_Genetic_Data/Exome_200k/BGEN_Format/output"


#### Default exclusions
exclusionsDefaultPath = "/data/studies/06_UKBB/01_Data/00_EXCLUSIONS/"





#########################################
######## prepare and get options ########
#########################################

parser = argparse.ArgumentParser()
parser.add_argument("-v", "--variantFile", help = "Path of input file containing chr and variant identifyer", action="store", required=True)
parser.add_argument("-m", "--mode", help = "possible values: \"geno\", \"imp\" or \"wes\"", action = "store", required = True)
parser.add_argument("-o", "--output", help = "name of output file", action = "store", required = True)
parser.add_argument("-e", "--exclude", help = "Path to exclude list. One PatId per line", action = "store", required = False)


parser.add_argument("-c", "--cmdOnly", help = "prints out execution commands instead of running them.", action = "store_true", required = False)



args = parser.parse_args()


rsFileInPath = args.variantFile
outFileName = args.output
mode = args.mode
excludeListPath = args.exclude
cmdOnly = args.cmdOnly



#### check if excludeListPath is set else use default file
if not excludeListPath:
	# print (glob.glob(exclusionsDefaultPath + "/*"))
	excludeListPath = exclusionsDefaultPath + "/" + sorted(listdir(exclusionsDefaultPath))[-1]
	args.exclude = excludeListPath





########################
#### prepare folder ####
########################

tmpDir = "tmp"
outDir = os.path.dirname(outFileName)
if not outDir:
	outDir = os.getcwd()

if not cmdOnly:
	if os.path.exists(tmpDir):
		shutil.rmtree(tmpDir)
	os.mkdir(tmpDir)
	if (not os.path.exists(outDir)):#
		os.makedirs(outDir)



################################
######## define methods ########
################################


##################################
######## prepare out file ########
def writeOutFromGen():
	genFiles = glob.glob(tmpDir + "/" + "*.gen")
	colNameList = list()
	rowList = list()
	headerList = list()
	headerList.append("Individual_ID")

	for genFile in genFiles:

		sampleFile = genFile.replace(".gen", ".sample")

		#### read in gen file and sample file and create matrix
		genFileIn = open(genFile, "r")
		genFile = genFileIn.readlines()
		genFileIn.close()

		sampleFileIn = open(sampleFile, "r")
		sampleFile = sampleFileIn.readlines()[2:]
		sampleFileIn.close();



		## prepare rowHeader (indiv id)
		if (len(colNameList) == 0):
			for curLine in sampleFile:
				colNameList.append(curLine.split()[0])


		## cut chr rsid pos ref alt from file and store as header
		for curLine in genFile:

			## prepare colHeader
			splitLine = curLine.split()
			splitLine.pop(1)
			splitLine[0] = "chr" + splitLine[0]
			tmpList = splitLine[:5]
			del(splitLine[:5])
			header = "_".join(tmpList)
			headerList.append(header)

			## calculate dosage entries
			dosageList = list()
			n = 3 # skip interval
			for index in range(len(splitLine)):
				if index % n == 0:
					dos = float(splitLine[index+1]) + 2 * float(splitLine[index+2])
					if (dos.is_integer()):
						dosageList.append(round(dos))
					else:
						dosageList.append(dos)
			rowList.append(dosageList)





	# write to outFile
	excludeList = readInExcludeList()
	outFile = open(outFileName, "w")

	print("\t".join(headerList), file = outFile)
	for i in range(len(colNameList)):
		patId = colNameList[i]
		if patId in excludeList:
			continue
		# print (colNameList[i])
		outFile.write(colNameList[i] + "\t")
		curOutRow = list()
		for j in range(len(rowList)):
			curOutRow.append(str(rowList[j][i]))
		print("\t".join(curOutRow), file = outFile)

	outFile.close()






###############################
######## read in exclusion file

def readInExcludeList():
	excludeList = list()
	if not excludeListPath is None and os.path.isfile(excludeListPath):
		with open(excludeListPath, "r") as input:
			for row in input:
				excludeList.append(row.rstrip())
	else:
		excludeList.append("")
	# print ("\n".join(excludeList))
	# sys.exit()
	return excludeList








#########################
######## read in var file

def readInVarFile():

	fileIn = open(rsFileInPath, "r")
	input = fileIn.readlines()
	fileIn.close()


	## prepare map sorting rsids to chr
	rsidChrMappingDict = collections.defaultdict(list)
	for curLine in input:

		if curLine.startswith("#") or curLine.isspace():
			continue

		splitLine = curLine.split()
		chr = splitLine[0]
		rsId = splitLine[1]
		rsidChrMappingDict[chr].append(rsId)
		rsIdOutFile = tmpDir + "/chr" + chr +".rsid"

	return rsidChrMappingDict






##############################
######## transform bgen to gen

def bgenToGen(bgenInFile, sampleInFile, genOutFile, sampleOutFile):
	qcToolCmd = []
	qcToolCmd.append(qcToolPath)
	qcToolCmd.append("-g " + bgenInFile)
	qcToolCmd.append("-s " + sampleInFile)
	qcToolCmd.append("-og " + genOutFile)
	qcToolCmd.append("-os " + sampleOutFile)

	if not cmdOnly:
		returnValue = subprocess.call(" ".join(qcToolCmd), shell = True)
		print(" ".join(qcToolCmd))
		# if returnValue != 0:
		# 	sys.exit("\n\n######## ERROR ########\nQcTool exited with non 0 status (" + str(returnValue) + "). Check your parameters.\n\n")
	if cmdOnly:
		print(" ".join(qcToolCmd))
		print ("\n")







#################################
######## extract SNPs using plink

def extractUsingPlink(bedFile, bimFile, famFile, rsidChrMappingDict, bgenTmpFile):

	plinkCmd = list()
	plinkCmd.append(plinkPath)
	plinkCmd.append("--bed " + bedFile)
	plinkCmd.append("--bim " + bimFile)
	plinkCmd.append("--fam " + famFile)

	# rsIdList = list()
	# for rsId in rsidChrMappingDict[chr]:
	# 	rsIdList.append(rsId)
	plinkCmd.append("--snps " + ",".join(rsidChrMappingDict[chr]))
	plinkCmd.append("--export bgen-1.1")
	plinkCmd.append("--out " + bgenTmpFile)

	if not cmdOnly:
		returnValue = subprocess.call(" ".join(plinkCmd), shell = True)
		if returnValue != 0:
			sys.exit("\n\n######## ERROR ########\nPlink exited with non 0 status (" + str(returnValue) + "). Check your parameters.\n\n")

	if cmdOnly:
		print(" ".join(plinkCmd))
		print ("\n")



######################################
######## check if varIds in input file

def  checkForVarIdsInInputfile (rsidChrMappingDict, isGeno):


	notMatchedVarFileName = outDir + "/excluded_var.txt"
	notMatchedVarFile = open(notMatchedVarFileName, "w")
	falseKeyList = list()
	for chr in rsidChrMappingDict:

		if isGeno:
			#### get correstponding bim files
			bimFile = glob.glob(genoPath + "/*chr" + chr + "_*.bim")[0]
		else:
			bimFile = wesBimFile

		## check if rs file exists in bimFile
		inputRsIdList = rsidChrMappingDict[chr]
		bimFileVarSet = set()
		excludeIdList = list()
		with open(bimFile, "r") as bimFileIn:
			for line in bimFileIn:
				bimFileVarSet.add(line.split("\t")[1])

		for curVar in inputRsIdList:
			if curVar not in bimFileVarSet:
				excludeIdList.append(curVar)

		if len(excludeIdList) > 0:
			for curExcludeVar in excludeIdList:
				inputRsIdList.remove(curExcludeVar)
		if len(inputRsIdList) > 0:
			rsidChrMappingDict[chr] = inputRsIdList
		else:
			falseKeyList.append(chr)

		print("\n".join(excludeIdList), file = notMatchedVarFile)
	notMatchedVarFile.close()

	for chr in falseKeyList:
		del rsidChrMappingDict[chr]
	return rsidChrMappingDict








def checkForVarIdsInOutputFile (rsidChrMappingDict):

	notExtractedVar = list()

	outFileIn = open(outFileName, "r")
	splitHeaderLine = outFileIn.readline().split("\t")
	isFirst = True
	outVarSet = set()
	for curHeader in splitHeaderLine:
		if isFirst:
			isFirst = False
			continue

		outVarSet = curHeader.split("_")[1]


	for chr in rsidChrMappingDict:
		for curInputVar in rsidChrMappingDict[chr]:
			if not curInputVar in outVarSet:
				notExtractedVar.append(curInputVar)


	notMatchedVarFileName = outDir + "/excluded_var.txt"
	notMatchedVarFile = open(notMatchedVarFileName, "w")
	print ("\n".join(notExtractedVar), file = notMatchedVarFile)
	notMatchedVarFile.close()







###########################################
#### performe action depending on mode ####
###########################################


################################
######## running on imputed data

if mode == "imp":

	#### read in rsid file
	rsidChrMappingDict = readInVarFile()

	for chr in rsidChrMappingDict:

		#### extract snps using bgenix
		bgenixCmd = []
		bgenFile =  glob.glob(impPath + "*chr" + chr + "_*.bgen")[0]
		bgenIndexFile = glob.glob(impPath + "*chr" + chr + "_*.bgi")[0]
		outFilePrefix = "chr" + chr
		bgenOutFile = tmpDir + "/" + outFilePrefix + ".bgen"

		bgenixCmd.append(bgenixPath)
		bgenixCmd.append("-g " + bgenFile)
		bgenixCmd.append("-i " + bgenIndexFile)
		for rsId in rsidChrMappingDict[chr]:
			bgenixCmd.append("-incl-rsids " + rsId)
		bgenixCmd.append("> " + bgenOutFile)
		separator = " "
		if not cmdOnly:
			returnValue = subprocess.call(" ".join(bgenixCmd), shell=True)
			if returnValue != 0:
				sys.exit("\n\n######## ERROR ########\nBGenix exited with non 0 status (" + str(returnValue) + "). Check your parameters.\n\n")
		if cmdOnly:
			print (" ".join(bgenixCmd))
			print ("\n")


		#### transform bgen to gen and sample using gcTool
		genOutFile = tmpDir + "/" + outFilePrefix + ".gen"
		sampleOutFile = tmpDir + "/" + outFilePrefix + ".sample"

		bgenToGen(bgenOutFile, impSampleFile, genOutFile, sampleOutFile)
		writeOutFromGen()
		checkForVarIdsInOutputFile(rsidChrMappingDict)





##############################
######## run on genotyped data

if mode == "geno":

	rsidChrMappingDict = readInVarFile()

	rsidChrMappingDict = checkForVarIdsInInputfile(rsidChrMappingDict, True)

	for chr in rsidChrMappingDict:

		#### get correstponding bed/bim files
		bedFile = glob.glob(genoPath + "/*chr" + chr + "_*.bed")[0]
		bimFile = glob.glob(genoPath + "/*chr" + chr + "_*.bim")[0]

		bgenTmpFile = tmpDir + "/chr" + chr
		sampleOutFile = tmpDir + "/chr" + chr + ".sample"
		bgenOutFile = tmpDir + "/chr" + chr + ".bgen"
		genOutFile = tmpDir + "/chr" + chr + ".gen"

		extractUsingPlink(bedFile, bimFile, genoFamPath, rsidChrMappingDict, bgenTmpFile)

		bgenToGen(bgenOutFile, sampleOutFile, genOutFile, sampleOutFile)

		writeOutFromGen()




if mode == "wes":

	rsidChrMappingDict = readInVarFile()
	# rsidChrMappingDict = checkForVarIdsInInputfile(rsidChrMappingDict, False)

	for chr in rsidChrMappingDict:
		bgenixCmd = []
		bgenFile =  glob.glob(wesPath + "/*chr"+ chr + ".bgen")[0]
		bgenIndexFile = glob.glob(bgenFile + ".bgi")[0]
		outFilePrefix = "chr" + chr
		bgenOutFile = tmpDir + "/" + outFilePrefix + ".bgen"
		sampleInFile = glob.glob(wesPath + "/*chr"+ chr + ".sample")[0]

		bgenixCmd.append(bgenixPath)
		bgenixCmd.append("-g " + bgenFile)
		bgenixCmd.append("-i " + bgenIndexFile)
		for rsId in rsidChrMappingDict[chr]:
			bgenixCmd.append("-incl-rsids " + rsId)
		bgenixCmd.append("> " + bgenOutFile)
		separator = " "
		if not cmdOnly:
			returnValue = subprocess.call(" ".join(bgenixCmd), shell=True)
			if returnValue != 0:
				sys.exit("\n\n######## ERROR ########\nBGenix exited with non 0 status (" + str(returnValue) + "). Check your parameters.\n\n")
		if cmdOnly:
			print (" ".join(bgenixCmd))
			print ("\n")

		#### transform bgen to gen and sample using gcTool
		genOutFile = tmpDir + "/" + outFilePrefix + ".gen"
		sampleOutFile = tmpDir + "/" + outFilePrefix + ".sample"

		bgenToGen(bgenOutFile, sampleInFile, genOutFile, sampleOutFile)
		writeOutFromGen()
		checkForVarIdsInOutputFile(rsidChrMappingDict)


		######## old version using bed bim fam
		# print (rsidChrMappingDict)
		# bgenTmpFile = tmpDir + "/chr" + chr
		# sampleOutFile = tmpDir + "/chr" + chr + ".sample"
		# bgenOutFile = tmpDir + "/chr" + chr + ".bgen"
		# genOutFile = tmpDir + "/chr" + chr + ".gen"
		#
		# extractUsingPlink(wesBedFile, wesBimFile, wesFamFile, rsidChrMappingDict, bgenTmpFile)
		#
		# print (wesPath)

		# bgenToGen(bgenInFile, sampleInFile, genOutFile, sampleOutFile)
		#
		# writeOutFromGen()



print ("\n#### Extraction done! ####")
