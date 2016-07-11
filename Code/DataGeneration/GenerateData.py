import numpy as np
import math

#hyper-parameters
sigma_0 = 1
sigma_1 = 1.5

def generate_data(n, maf, K):
	"""Incomplete
	"""
	gene_data = []
	return gene_data

def generate_gene_data(n, maf):
	"""Generates the gene data for n people. 

	   E.g.
	   Suppose n=100 and maf = 0.1.
	   Then total number of alleles = 200 while total number of minor alleles is about 20.

       Given n, create a list of length n. For each entry in list, generate 2 random numbers independently and 
       uniformly in [0, 1]. The zygosity of this entry is equal to (random1 < maf) + (random2 < maf).
	"""
	gene_data = []
	for i in range(0, n):
		random1 = np.random.uniform()
		random2 = np.random.uniform()
		gene_data.append((random1 < maf) + (random2 < maf)) #code detail: note that in python (True) + (True) = 2
		total = total + (random1 < maf) + (random2 < maf)
	return gene_data


"""
sigma_0 = 1
sigma_1 = 1
position = 2 # from 0 to 3

def generateInputX(lociCount, n):
	Generates genotype data for n people for a particular gene
	   lociCount is the number of loci for the gene
	   Returns a n by K matrix (as a list)
	
	inputX = []
	for j in range(0, n):
		rawData = np.random.randint(3, size=lociCount)
		inputX.append(rawData.tolist())
	return inputX

def generateEta(inputX, K):
	Generates eta based on inputX. 
	   Each person has K eta.
	
	Eta = [] # eta for every one
	n = len(inputX)

	for j in range(0, n):
		individualEta = [] # eta for a single person
		for k in range(0, K):
			exponent = inputX[j][position] * np.random.normal(0, sigma_0) + np.random.normal(0, sigma_1)
			individualEta.append(math.exp(exponent))
		Eta.append(individualEta)
	return Eta

def generateY(inputEta):
	n = len(inputEta)
	K = len(inputEta[0])

	inputY = []
	for j in range(0, n):
		arrayEta = np.asarray(inputEta[j])
		inputY.append((np.random.dirichlet(arrayEta).tolist()))
	return inputY


resultX = generateInputX(4, 30)
print(resultX)
print("\n")

resultEta = generateEta(resultX, 5)
print(resultEta)
print("\n")
resultY = generateY(resultEta)
print(resultY)
"""