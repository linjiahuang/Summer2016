import numpy as np
import math
from random import shuffle
import os as os

sigma_0 = 1.0                         #std dev for beta_k
sigma_1 = 1.5                         #std dev for eta_jk
lower_bound = -1.0                    #lower bound for beta_0k
upper_bound = 1.0                     #upper bound for beta_0k

def generate_gene_data(n, maf):
	"""Generates the gene data for n people with minor allele frequency maf. 

       E.g.
	   Suppose n=100 and maf=0.1.
	   Then total number of alleles=100*2=200 while total number of minor alleles is about 200*0.1=20.

       Returns a list of length n with entries 0, 1 or 2.
	"""
	x = np.random.binomial(2,maf,n)
	# ensures that x contains all 3 types of data (i.e. 0, 1 and 2)
	# note that np.where(x==0) is the same as (x==0).nonzero(), which returns a tuple of arrays (in our case
	# of the form array([1, 2, 3, 4]),) if indexes 1, 2, 3, 4 are nonzero.). Therefore, we need to take the 
	# first element [0].
	while (len(np.where(x==0)[0])*len(np.where(x==1)[0])*len(np.where(x==2)[0])==0):
		x =  np.random.binomial(2,maf,n)
	return x.tolist()

def generate_beta_0k(lower_bound, upper_bound, K):
	"""Helper method for generating intercept data. Uses uniform distribution.

	   Returns a list of length K.
	"""
	return np.random.uniform(lower_bound, upper_bound, K).tolist()

def generate_beta_k(sigma, K):
	"""Helper method for generating effect data. Uses normal distribution.

	   Returns a list of length K.
	"""
	return np.random.normal(0, sigma, K).tolist()

def generate_epsilon_j_k(sigma, n, K):
	"""Helper method for generating noise data. Uses normal distribution.

	   Returns a list of dimension n by K.
	"""
	return np.random.normal(0, sigma, (n, K)).tolist()

def generate_correlated_data(gene_data, beta_0k, beta_k, epsilon_j_k):
	"""Generates the K isoform proportions for n individuals, for use in the alternative hypothesis.

	   1) gene_data is a 1 by n list of minor allele count for the n individuals 
	   2) beta_k is a 1 by K list containing the effect, values drawn from normal with mean 0 and std dev sigma_0 
	   3) beta_0k is a 1 by K list containing the intercept, values drawn uniformly from -lower_bound to upper_bound
	   4) epsilon_j_k is an n by K list containing the noise, values drawn from drawn from normal with mean 0 and std dev sigma_1 
	   5) values from steps 1-4 are then used to create an n by K list containing eta, the Dirichlet pdf parameters
	   6) the Dirichlet pdf parameters then generate the isoform proportions, according to the equation in the write-up.

	   Returns an n by K list of isoform proportions such that for each individual, the K isoforms sum to 1.0.
	"""
	n = len(gene_data) #number of individuals
	K = len(beta_0k)   #number of isoforms

	Eta = [] #n by K list 
	for j in range(0, n):
		eta = [] #1 by K list
		for k in range(0, K):
			eta.append(math.exp(beta_0k[k] + gene_data[j]*beta_k[k] + epsilon_j_k[j][k]))
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

def create_folders():
	directory = '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/Simulation Result'
	n = [10, 100, 500]
	maf = [0.05, 0.1, 0.2]
	K = [3, 10, 100]

	for n_item in n:
		for maf_item in maf:
			for K_item in K:
				new_directory = directory + '/file_n=' + str(n_item) + '-maf=' + str(int(maf_item * 100)) + '-K=' + str(K_item) 
				os.mkdir(new_directory)


def print_data(trial_number, n, K, maf):
	#generate appropriate data
	gene_data = generate_gene_data(n, maf)
	beta_0k = generate_beta_0k(lower_bound, upper_bound, K)
	beta_k = generate_beta_k(sigma_0, K)
	epsilon_j_k = generate_epsilon_j_k(sigma_1, n, K)
	Isoform = generate_correlated_data(gene_data, beta_0k, beta_k, epsilon_j_k)

	filepath = '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/Simulation Result/file_n=' + str(n) + '-maf=' + str(int(maf * 100)) + '-K=' + str(K) + '/'

	#open file for simulation
	maf_formated = maf * 100
	filename_simulation = 'simulation' + str(trial_number) + '_n=' + str(n) + '-maf=' + str(int(maf_formated)) + '-K=' + str(K) + '.txt'
	simulation = open(filepath + filename_simulation, 'w')
	for i in range(0, n):
		simulation.write('{0:1d}'.format(gene_data[i]))
		simulation.write('\t')
		for k in range(0, K-1):
			simulation.write('{0:.10f}'.format(Isoform[i][k]))
			simulation.write('\t')
		simulation.write('{0:.10f}'.format(Isoform[i][K-1]))
		simulation.write('\n')
	simulation.close()

	#shuffles for permutation
	shuffle(Isoform)

	#opens file for permutation
	filename_permutation = 'permutation' + str(trial_number) + '_n=' + str(n) + '-maf=' + str(int(maf_formated)) + '-K=' + str(K) + '.txt'
	permutation = open(filepath + filename_permutation, 'w')
	for i in range(0, n):
		permutation.write('{0:1d}'.format(gene_data[i]))
		permutation.write('\t')
		for k in range(0, K-1):
			permutation.write('{0:.10f}'.format(Isoform[i][k]))
			permutation.write('\t')
		permutation.write('{0:.10f}'.format(Isoform[i][K-1]))
		permutation.write('\n')
	permutation.close()

	#opens files for specification
	#as of now only specifies intercept and effect, not noise
	filename_specification = 'specification' + str(trial_number) + '_n=' + str(n) + '-maf=' + str(int(maf_formated)) + '-K=' + str(K) + '.txt'
	specification = open(filepath + filename_specification, 'w')
	specification.write('Values used for beta_0k(intercept). Sampled from {0} to {1} uniformly:\n'.format(lower_bound, upper_bound))
	for k in range(0, K):
		specification.write('{0:4d}\t{1:.10f}'.format(k+1, beta_0k[k]))
		specification.write('\n')
	specification.write('\n')
	specification.write('Values used for beta_k(effect). Sampled from normal with mean 0 and std {0}:\n'.format(sigma_0))
	for k in range(0, K):
		specification.write('{0:4d}\t{1:.10f}'.format(k+1, beta_k[k]))
		specification.write('\n')
	specification.close()
	"""
	for i in range(0, n):
		print('{0:1d}'.format(gene_data[i]), end='\t')
		for k in range(0, K-1):
			print('{0:.10f}'.format(Isoform[i][k]), end='\t')
		print('{0:.10f}'.format(Isoform[i][K-1]))
	"""

def print_all_data():
	n = [10, 100, 500]
	maf = [0.05, 0.1, 0.2]
	K = [3, 10, 100]

	for n_item in n:
		for maf_item in maf:
			for K_item in K:
				for i in range(0, 500):
					print_data(i+1, n_item, K_item, maf_item)

#create_folders()
print_all_data()































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