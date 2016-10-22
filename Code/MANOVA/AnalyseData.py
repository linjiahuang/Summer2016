import math
from random import shuffle
import matplotlib.pyplot as plt

#gene_data is a 1 by n list
#Isoform is a n by K list
#Tested on Oct 22 and worked
def getFileData(filename, n, K):
	gene_data = []
	Isoform = []
	with open(filename) as f:
		for line in f:
			raw_data = map(float, line.split())
			gene_data.append(raw_data[0])
			Isoform.append(raw_data[1:K+1])
	return [gene_data, Isoform]

"""
#For Testing
filename = '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/New Simulation Result/file_n=' + str(10) + '-maf=' + str(int(0.05 * 100)) + '-K=' + str(3) + '/'
filepath = filename + 'simulation' + str(8) + '_n=' + str(10) + '-maf=' + str(int(5)) + '-K=' + str(3) + '.txt'
print(getFileData(filepath, 10, 3))
"""

#Each point is a list of doubles. For each point, the doubles sum to 1.
#Note that this is the Hellinger distance SQUARED
#Tested on Oct 22 and worked
def getHellingerSquared(point_a, point_b):
	result = 0
	for a, b in list(zip(point_a, point_b)):
		result = result + (math.sqrt(a) - math.sqrt(b))**2
	return result

def getSS(Isoform):
	n = len(Isoform)
	result = 0
	for i in range(0, n-1):
		for j in range(i+1, n):
			result = result + getHellingerSquared(Isoform[i], Isoform[j])
	return float(result)/n

#x=0 for minor minor, x=1 or minor major, x=2 for major major
def getSS_W(aggregated_data):
	minorminor_aggregated_data = []
	minormajor_aggregated_data = []
	majormajor_aggregated_data = []

	n = len(aggregated_data[1])
	for i in range(0, n):
		if aggregated_data[0][i] == 0:
			minorminor_aggregated_data.append(aggregated_data[1][i])
		if aggregated_data[0][i] == 1:
			minormajor_aggregated_data.append(aggregated_data[1][i])
		if aggregated_data[0][i] == 2:
			majormajor_aggregated_data.append(aggregated_data[1][i])
	return getSS(minorminor_aggregated_data) + getSS(minormajor_aggregated_data) + getSS(majormajor_aggregated_data)

#Given one file, calculate the F score
#Tested on Oct 22 and worked
def getFScore(filename, n, K):
	aggregated_data = getFileData(filename, n, K)
	SS_T = getSS(aggregated_data[1])
	SS_W = getSS_W(aggregated_data)
	F_original = ((SS_T - SS_W)/(K - 1.0))/(SS_W/(n - K))

	shuffle(aggregated_data[0])

	SS_T = getSS(aggregated_data[1])
	SS_W = getSS_W(aggregated_data)
	F_permuted = ((SS_T - SS_W)/(K - 1.0))/(SS_W/(n - K))

	return [F_original, F_permuted]

#Tested on Oct 22 and worked
def printFScore(n, K, maf):
	simulation_result = []

	for trial_number in range(1, 501):
		filepath = '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/New Simulation Result/file_n=' + str(n) + '-maf=' + str(int(maf * 100)) + '-K=' + str(K) + '/'
		#open file for simulation
		maf_formated = maf * 100
		filename_simulation = 'simulation' + str(trial_number) + '_n=' + str(n) + '-maf=' + str(int(maf_formated)) + '-K=' + str(K) + '.txt'
		FScore = getFScore(filepath + filename_simulation, n, K)
		simulation_result.append(FScore)

	simulation_result.sort(key = lambda x: x[0], reverse=True)

	total_trial = len(simulation_result)
	with open(filepath + 'simulation_result.txt', "a") as myfile:
		for i in range(0, total_trial):
			myfile.write(str(i+1))
			myfile.write('\t')
			myfile.write(str(simulation_result[i][0]))
			myfile.write('\t')
			myfile.write(str(simulation_result[i][1]))
			myfile.write('\n')

def precisionRecallSetup(n, K, maf):
	filepath = '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/New Simulation Result/file_n=' + str(n) + '-maf=' + str(int(maf * 100)) + '-K=' + str(K) + '/'
	#open file for simulation
	precisionRecall_data = []

	with open(filepath + 'simulation_result.txt') as myfile:
		for line in myfile:
			raw_data = map(float, line.split())
			if(raw_data[1] > raw_data[2]):
				precisionRecall_data.append(1)
			else:
				precisionRecall_data.append(0)

	return precisionRecall_data

#cut_off starts from 1
def getPrecisionRecallValues(precisionRecall_data, cut_off):
	Tp = 0
	Tn = 0
	Fp = 0
	Fn = 0

	for i in range(0, cut_off):
		if(precisionRecall_data[i] == 1):
			Tp = Tp + 1.0
		else:
			Fp = Fp + 1.0

	for i in range(cut_off, len(precisionRecall_data)):
		if(precisionRecall_data[i] == 1):
			Fn = Fn + 1.0
		else:
			Tn = Tn + 1.0

	return [Tp/(Tp + Fp), Tp/(Tp + Fn)]

def plotPrecisionRecall(n, K, maf):
	precisionRecall_data = precisionRecallSetup(n, K, maf)
	precision_list = []
	recall_list = []
	for i in range(1, len(precisionRecall_data)):
		PrecisionRecallValues = getPrecisionRecallValues(precisionRecall_data, i)
		precision_list.append(PrecisionRecallValues[0])
		recall_list.append(PrecisionRecallValues[1])
	plt.title(str(n) + ' individuals, ' + str(K) +  ' isoforms, ' + str(maf) + ' maf')
	plt.plot(recall_list, precision_list, 'ro')
	plt.ylabel('Precision')
	plt.xlabel('Recall')
	plt.axis([0, 1.0, 0, 1.0])
	plt.savefig(str(n) + ' individuals, ' + str(K) +  ' isoforms, ' + str(maf) + ' maf' + '.png')
	plt.close()
	#plt.show()

n = [10, 100]
maf = [0.05, 0.1, 0.2]
K = [3, 5]

"""
for n_item in n:
	for maf_item in maf:
		for K_item in K:
			printFScore(n_item, K_item, maf_item)
"""
"""
for n_item in n:
	for maf_item in maf:
		for K_item in K:
			plotPrecisionRecall(n_item, K_item, maf_item)
"""
plotPrecisionRecall(100, 5, 0.2)