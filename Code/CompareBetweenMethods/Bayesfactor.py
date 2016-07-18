"""MCMC calculate Bayesfactor, for checking with Laplace approximation
"""

import sys
sys.path.insert(0, '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/Code')

import numpy as np
import scipy.special as special
from MultivariateBeta import multivariateBeta # for returning the multivariate Beta function
import math

numOfTrials = 100000
Sigma_0 = 1
Sigma_1 = 1
test_position = 0 #0 to 3

y_1 = 0.001
y_2 = 0.234
y_3 = 0.676
y_4 = 1.0 - y_1 - y_2 - y_3

generated_data_inputX = [[1, 0, 0, 1]]
generated_data_inputY = [[y_1, y_2, y_3, y_4]]

def integrand_beta(inputEta):
	"""Returns the product of the 1/B() term in the integral
	"""
	n = len(inputEta)

	result = 1
	for j in range(0, n):
		result = result * multivariateBeta(inputEta[j])

	return result

def integrand_y(inputY, inputEta):
	n = len(inputEta)
	K = len(inputEta[0])

	result = 1
	for j in range(0, n):
		for k in range(0, K):
			if inputY[j][k] == 0:
				return 0
			result = result * inputY[j][k]**(inputEta[j][k] - 1)

	return result

def normal_rv_mu(inputY):
	"""Draws random variable mu from standard normal
	"""
	n = len(inputY)
	K = len(inputY[0])

	Mu = []
	for j in range(0, n):
		Mu_j = []
		for k in range(0, K):
			mu = np.random.normal(0,1)
			Mu_j.append(mu)
		Mu.append(Mu_j)

	return Mu

def mu_to_eta_null_hypothesis(inputMu, sigma_1):
	n = len(inputMu)
	K = len(inputMu[0])

	Eta = []
	for j in range(0, n):
		Eta_j = []
		for k in range(0, K):
			temp = inputMu[j][k] * sigma_1
			Eta_j.append(math.exp(temp))
		Eta.append(Eta_j)

	return Eta

def mu_to_eta_alt_hypothesis(inputX, inputMu, sigma_0, sigma_1):
	n = len(inputMu)
	K = len(inputMu[0])

	Eta = []
	for j in range(0, n):
		Eta_j = []
		for k in range(0, K):
			temp = inputMu[j][k] * math.sqrt((inputX[j][test_position]**2) * (sigma_0**2) + sigma_1**2)
			Eta_j.append(math.exp(temp))
		Eta.append(Eta_j)

	return Eta

def bayes_factor(inputX, inputY):
	totalSum_null = 0
	totalSum_alt = 0

	for i in range(0, numOfTrials):
		mu_standard_normal = normal_rv_mu(inputY)
		Eta_null_hypothesis = mu_to_eta_null_hypothesis(mu_standard_normal, Sigma_1)
		Eta_alt_hypothesis = mu_to_eta_alt_hypothesis(inputX, mu_standard_normal, Sigma_0, Sigma_1)

		totalSum_null = totalSum_null + integrand_beta(Eta_null_hypothesis) * integrand_y(inputY, Eta_null_hypothesis)
		totalSum_alt = totalSum_alt + integrand_beta(Eta_alt_hypothesis) * integrand_y(inputY, Eta_alt_hypothesis)

	numerator = totalSum_alt/numOfTrials
	denominator = totalSum_null/numOfTrials

	print("Numerator:    ", numerator)
	print("Denominator:  ", denominator)
	print("Bayes Factor: ", numerator/denominator)
bayes_factor(generated_data_inputX, generated_data_inputY)

