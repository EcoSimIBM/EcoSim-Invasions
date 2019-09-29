import csv
import sys
import string

#Made by Ryan Scott (scotto@uwindsor.ca)
#Any questions or ideas just let me know
#Simple program that takes a set of CSV files, extracts only the relevent column from them, and then creates another two columns for average and sdev of the given values
#This is a python program so it works in Windows too, as long as you have enough RAM to load the results files one-by-one
#Only problem is I wrote the other script in Bash, so to automate everything you'll have to write your own stuff
	
	
#sys.argv[1] is the structure of the input file name, using '?' as a wildcard
#sys.argv[2] is the number of files to iterate through
#sys.argv[3] is the directory in which to put the result file
#sys.argv[4] is the lowest number with which to iterate
#sys.argv[5] is the incrementation
#sys.argv[6] is appended to the output file name

#example: python AverageRuns.py "/work/nap_9022/Ryan/Brian/FluxResourcesFair/Original/Original?/Results_Prey_Female.csv" 15 /work/nap_9022/Ryan/Brian/FluxResourcesFair/Original/ 1 1 FemalePrey

def translateFileName(originalFileName, num):
	return string.replace(originalFileName, "?", str(num))

outputTemp = []
upperNumber = int(sys.argv[2])
new_file = sys.argv[3] + "ResultExtracted_" + sys.argv[6] + ".csv"
x = int(sys.argv[4])
incrementation = int(sys.argv[5])
maxLength = 0
while x <= upperNumber: #loop across all files
	filenameStructure = sys.argv[1] #input file for this iteration
	fileName = translateFileName(filenameStructure, x) #insert the number where the ? is
	currentFile = open (fileName)
	csv_reader = csv.reader(currentFile)
	transpose = zip(*csv_reader)
	transpose = list(transpose)
	maxLength = max(maxLength, len(transpose[0]))
	outputTemp.append(transpose)
	
	currentFile.close()
	x+=incrementation

output = []


for j in range(len(outputTemp[0])): #iterate across columns
	outputCol = []
	for k in range(maxLength): #iterate through all rows
		sum = float(0)
		elements = int(0)
		for i in range(len(outputTemp)): #iterate across files
			try:
				if (k == 0 and i == 0):
					outputCol.append(outputTemp[i][j][k])
				elif (k == 0):
					print "..."
				elif (outputTemp[i][j][k] != None):
					sum += float(outputTemp[i][j][k]) 
					elements = elements + 1
			except:
				sum += 0
				elements += 0
		if (elements > 0):
			outputCol.append(float(float(sum)/float(elements)))
	output.append(outputCol)

	
output = list(zip(*output))
with open(new_file, 'w') as csvfile2:
	spamwriter = csv.writer(csvfile2, lineterminator='\n')
	for i in range(len(output)):
		spamwriter.writerow(output[i])
		
		
		