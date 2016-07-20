import scipy.special as special

def multivariateBeta(inputEta):
	"""Calculates multivariate 1/Beta, not Beta.
	   inputEta is a list of floats
	"""
	currSum = inputEta[0]
	n = len(inputEta)

	currProduct = 1
	for i in range(1, n):
		value = special.beta(currSum,inputEta[i])
		"""
		if (value < 0.00000001):
			value = 0.00000001
		if (value > 1000000000):
			value = 1000000000
		"""
		currProduct = currProduct * (1/value)
		currSum = currSum + inputEta[i]

	return currProduct

def multivariateLogBeta(inputEta):
	"""Calculate the natural log of multivariate
	"""
	currSum = inputEta[0]
	n = len(inputEta)

	currValue = 0
	for i in range(1, n):
		value = -special.betaln(currSum,inputEta[i])
		"""
		if (value < 0.00000001):
			value = 0.00000001
		if (value > 1000000000):
			value = 1000000000
		"""
		currValue = currValue + value
		currSum = currSum + inputEta[i]

	return currValue

"""
inputEta = [10, 3, 2, 6, 8]

def compBeta(inputEta):
	#Function to test if my multivariateBeta is implemented correctly
	
	currSum = 0
	n = len(inputEta)
	for i in range(0, n):
		currSum = currSum + inputEta[i]

	currProduct = 1
	for i in range(0, n):
		currProduct = currProduct * special.gamma(inputEta[i])

	return special.gamma(currSum)/currProduct

print(multivariateBeta(inputEta))
print(compBeta(inputEta))
print(math.log(compBeta(inputEta)))
print(multivariateLogBeta(inputEta))
"""