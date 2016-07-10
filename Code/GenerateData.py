import numpy as np
import math

sigma_0 = 1
sigma_1 = 1
position = 2 # from 0 to 3

def generateInputX(lociCount, n):
	"""Generates genotype data for n people for a particular gene
	   lociCount is the number of loci for the gene
	   Returns a n by K matrix (as a list)
	"""
	inputX = []
	for j in range(0, n):
		rawData = np.random.randint(3, size=lociCount)
		inputX.append(rawData.tolist())
	return inputX

def generateEta(inputX, K):
	"""Generates eta based on inputX. 
	   Each person has K eta.
	"""
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
