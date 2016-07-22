import numpy as np
import math

#hyper-parameters for testing purposes
sigma_0 = 1.0              #std dev for beta
sigma_1 = 1.0              #std dev for eta
#beta_0k = [0.2, 0.7, 0.5]       #intercept

def generate_correlated_data(gene_data, sigma_0, sigma_1, beta_0k):
	"""Generates the K isoform proportions for n individuals, for use in the alternative hypothesis.
	   First creates a 1 by K list containing the effect.
	   Then creates an n by K list containing epsilon, the noise.
	   Then creates an n by K list containing eta, the Dirichlet pdf parameters.
	   Then use those parameters to generate the isoform proportions, according to the equation in the write-up.

	   Returns an n by K list of isoform proportions such that for each individual, the K isoforms sum to 1.0.
	"""
	n = len(gene_data) #number of individuals
	K = len(beta_0k)   #number of isoforms

	beta_normal = [] #1 by K list containing the effect
	for k in range(0, K):
		beta_normal.append(np.random.normal(0, sigma_0))

	Epsilon_normal = [] #n by K list containing the noise
	for j in range(0, n):
		epsilon_normal = []
		for k in range(0, K):
			epsilon_normal.append(np.random.normal(0, sigma_1))
		Epsilon_normal.append(epsilon_normal)

	Eta = [] #n by K list 
	for j in range(0, n):
		eta = [] #1 by K list
		for k in range(0, K):
			eta.append(math.exp(beta_0k[k] + gene_data[j]*beta_normal[k] + Epsilon_normal[j][k]))
		Eta.append(eta)

	Isoform = [] #n by K list
	for j in range(0, n):
		arrayEta = np.asarray(Eta[j])
		Isoform.append((np.random.dirichlet(arrayEta).tolist()))
	return Isoform

def generate_random_data(n, K, beta_0k):
	"""Generates the K isoform proportions for n individuals, for use in the null hypothesis.
	   First creates an n by K list containing eta, the Dirichlet pdf parameters
	   Then use those parameters to generate the isoform proportions.
	   Note that the n by K list for eta is independent of the genotype, since we are using this for the null hypothesis.

	   Returns a n by K list of isoform proportions such that for each individual, the K isoforms sum to 1.0.
	"""
	Eta = [] #n by K list 
	for j in range(0, n):
		eta = [] #1 by K list
		for k in range(0, K):
			eta.append(np.random.lognormal(beta_0k[k], sigma_1))
		Eta.append(eta)

	Isoform = [] #n by K list
	for j in range(0, n):
		arrayEta = np.asarray(Eta[j])
		Isoform.append((np.random.dirichlet(arrayEta).tolist()))
	return Isoform

def generate_gene_data(n, maf):
	"""Generates the gene data for n people. 

       E.g.
	   Suppose n=100 and maf=0.1.
	   Then total number of alleles=200 while total number of minor alleles is about 20.

       Returns a list of length n with entries 0, 1 or 2.
	"""
	x = np.random.binomial(2,maf,n)
	# ensures that x contains all 3 types of data (i.e. 0, 1 and 2)
	while (len(np.where(x==0)[0])*len(np.where(x==1)[0])*len(np.where(x==2)[0])==0):
		x =  np.random.binomial(2,maf,n)
	return x

def print_data(n, K, maf, sigma_0, sigma_1, beta_0k):
	gene_data = generate_gene_data(n, maf)
	random_data = generate_random_data(n, K, beta_0k)
	correlated_data = generate_correlated_data(gene_data, sigma_0, sigma_1, beta_0k)

	return [gene_data, random_data, correlated_data]
	"""
	for i in range(0, n):
		print(gene_data[i], end='\t')
		for k in range(0, K-1):
			print(correlated_data[i][k], end='\t')
		print(correlated_data[i][K-1])

	"""
#print_data(10, 2, 0.2, sigma_0, sigma_1, beta_0k)

































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