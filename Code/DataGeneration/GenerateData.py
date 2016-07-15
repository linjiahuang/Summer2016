import numpy as np
import math

#hyper-parameters for testing purposes
sigma_0 = 1                #std dev for beta
sigma_1 = 1.5              #std dev for eta
beta_0k = [2.0, 1.5, 1.0]  #intercept

def generate_correlated_data(gene_data, sigma_0, sigma_1, beta_0k):
	"""Incomplete
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

def generate_random_data(n, K):
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
			eta.append(np.random.lognormal(0, sigma_1))
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

       Given n, create a list of length n. For each entry in list, generate 2 random numbers independently and 
       uniformly in [0, 1]. The zygosity of this entry is equal to (random1 < maf) + (random2 < maf).

       Returns a list of length n with entries 0, 1 or 2.
	"""
	gene_data = []
	for i in range(0, n):
		random1 = np.random.uniform()
		random2 = np.random.uniform()
		gene_data.append((random1 < maf) + (random2 < maf)) #code detail: note that in python (True) + (True) = 2
	return gene_data

gene_data = generate_gene_data(5, 0.7)

print(gene_data)
print(generate_correlated_data(gene_data, sigma_0, sigma_1, beta_0k))




























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