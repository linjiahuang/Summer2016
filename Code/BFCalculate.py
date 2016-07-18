import numpy as np
import scipy.special as special
#import GenerateData
from returnBeta import multivariateBeta # for returning the multivariate Beta function
import math

numOfTrials = 100000
Sigma_0 = 1
Sigma_1 = 1
test_position = 0

def returnBeta(inputEta):
	"""Returns the 1/B() term in the integral
	"""
	n = len(inputEta)

	result = 1
	for j in range(0, n):
		result = result * multivariateBeta(inputEta[j])

	return result

def returnY(inputY, inputEta):
	n = len(inputEta)
	K = len(inputEta[0])

	result = 1
	for j in range(0, n):
		for k in range(0, K):
			if inputY[j][k] == 0:
				return 0
			result = result * inputY[j][k]**(inputEta[j][k] - 1)

	return result

def returnEta_denom(inputY, sigma):
	n = len(inputY)
	K = len(inputY[0])

	Eta = []
	for j in range(0, n):
		Eta_j = []
		for k in range(0, K):
			eta = np.random.lognormal(0,sigma)
			Eta_j.append(eta)
		Eta.append(Eta_j)

	return Eta

def returnEta_numer(inputX, inputY, sigma_0, sigma_1):
	n = len(inputY)	
	K = len(inputY[0])

	Eta = []
	for j in range(0, n):
		Eta_j = []
		for k in range(0, K):
			eta = np.random.lognormal(0,math.sqrt((inputX[j][test_position]**2) * (sigma_0**2) + sigma_1**2))
			Eta_j.append(eta)
		Eta.append(Eta_j)

	return Eta

def BFDenominator(inputY):
	totalSum = 0

	for i in range(0, numOfTrials):
		inputEta = returnEta_denom(inputY, Sigma_0)
		totalSum = totalSum + returnBeta(inputEta) * returnY(inputY, inputEta)

	return totalSum/numOfTrials

def BFNumerator(inputX, inputY):
	totalSum = 0

	for i in range(0, numOfTrials):
		inputEta = returnEta_numer(inputX, inputY, Sigma_0, Sigma_1)
		totalSum = totalSum + returnBeta(inputEta) * returnY(inputY, inputEta)

	return totalSum/numOfTrials

"""
# first argument is lociCount, second is number of individuals
resultX = GenerateData.generateInputX(4, 10)
print(resultX)
print("\n")

# second argument is number of isoforms
resultEta = GenerateData.generateEta(resultX, 2)
resultY = GenerateData.generateY(resultEta)

print(resultY)
print("\n")
"""

fixResultX = [[1, 0]]
fixResultY = [[0.1, 0.9]]

numerator = BFNumerator(fixResultX, fixResultY)
denominator = BFDenominator(fixResultY)
print("Numerator:    ", numerator)
print("Denominator:  ", denominator)
print("Bayes Factor: ", numerator/denominator)





