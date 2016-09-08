import sys
sys.path.insert(0, '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/Code')
sys.path.insert(0, '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/Code/DataGeneration')

import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize
import scipy.special as special
import MultivariateBeta as mb
import GenerateData as gd
import math
import numdifftools as nd

# Global constants
n = 10000000.0 # for use in Laplace approximation

"""
# variables for testing purposes
y_1 = 0.03
y_2 = 1.0 - y_1
x = 2.0
"""

sigma_0 = 1.0
sigma_1 = 1.0
beta_0k = [0.5, 0.4, -0.5, -0.1]

def objective_r_denom(input_value, *args):
	mu = [] # Dirichlete pdf Parameters (exp)
	eta = [] # Dirichlete pdf Parameters (eta in our paper)
	y = [] # isoform proportions

	for v in input_value:
		mu.append(v)
		eta.append(math.exp(v))
	for v in args:
		y.append(v)

	K = len(mu)
		
	log_beta_func = mb.multivariateLogBeta(eta)	
	log_y_exponents = 0
	log_exponentials = 0
	for i in range(0, K):
		log_y_exponents = log_y_exponents + math.log(y[i]) * (eta[i] - 1) 
		log_exponentials = log_exponentials - ((mu[i] - beta_0k[i])**2)/(2*(sigma_1**2))
	other_factors = 1.0/(math.sqrt(2 * math.pi) * sigma_1)**K

	return -1.0 * (1.0/n) * (log_beta_func + math.log(other_factors) + log_y_exponents + log_exponentials)

def objective_r_numer(input_value, *args):
	mu = []   # Dirichlete pdf Parameters (exp)
	eta = []  # Dirichlete pdf Parameters (eta in our paper)
	y = []    # isoform proportions

	for v in input_value:
		mu.append(v)
		eta.append(math.exp(v))
	for v in args[0]:
		y.append(v)
	x = args[1]

	K = len(input_value)

	variance = x**2 * sigma_0**2 + sigma_1**2

	log_beta_func = mb.multivariateLogBeta(eta)
	log_y_exponents = 0
	log_exponentials = 0
	for i in range(0, K):
		log_y_exponents = log_y_exponents + math.log(y[i]) * (eta[i] - 1) 
		log_exponentials = log_exponentials - ((mu[i] - beta_0k[i])**2)/(2*(variance))
	other_factors = 1.0/(math.sqrt(2 * math.pi) * variance)**K

	return -1.0 * (1.0/n) * (math.log(other_factors) + log_beta_func + log_y_exponents + log_exponentials)

def bf_denominator(y):
	n = len(y)
	input_value = [1.1] * n #starting values for optimization
	args = tuple(y)
	res = minimize(objective_r_denom, input_value, args=args, method='Nelder-Mead', tol=1e-12)
	H = nd.Hessian(objective_r_denom)(res.x, *args)
	return (2 * math.pi / n)**(n/2.0) * math.sqrt(1.0 / np.linalg.det(H)) * math.exp(-1.0 * n * objective_r_denom(res.x, *args)) 

def bf_numerator(y, x):
	n = len(y)
	input_value = [1.1] * n #starting values for optimization

	args = tuple([y, x])
	res = minimize(objective_r_numer, input_value, args=args, method='Nelder-Mead', tol=1e-12)
	H = nd.Hessian(objective_r_numer)(res.x, *args)
	return (2 * math.pi / n)**(n/2.0) * math.sqrt(1.0 / np.linalg.det(H)) * math.exp(-1.0 * n * objective_r_numer(res.x, *args)) 

def bayes_factor(y, x):
	n = len(y) #number of individuals

	BFdenominator = 1
	BFnumerator = 1
	BF = 1.0

	for j in range(0, n):
		currDenominator = bf_denominator(y[j])
		BFdenominator = BFdenominator * currDenominator
		currNumerator = bf_numerator(y[j], x[j])
		BFnumerator = BFnumerator * currNumerator
		BF = BF * (currNumerator/currDenominator)

	print("Bayes Factor is : ", BF)

for i in range(0, 10):
	generated_data = gd.print_data(100, 4, 0.1, sigma_0, sigma_1, beta_0k)
	print("Correlated:")
	bayes_factor(generated_data[2], generated_data[0])
	print("Random    :")
	bayes_factor(generated_data[1], generated_data[0])
	print("------------------------------------")













"""
x = [0, 1, 0, 1, 0, 0, 2, 1, 1, 0]
random_y = [[0.2184238486801664, 0.7815761513198336], [0.7187021998226819, 0.2812978001773181], [0.3914529281843244, 0.6085470718156757], [0.00016211308888103415, 0.9998378869111189], [0.24824453394912138, 0.7517554660508786], [0.519195682930413, 0.48080431706958693], [0.9135355892859972, 0.08646441071400271], [0.680937107703438, 0.3190628922965621], [0.1840592144568694, 0.8159407855431307], [0.04374333747271165, 0.9562566625272884]]

correlated_y = [[0.8987624623472313, 0.10123753765276869], [0.033675225728584564, 0.9663247742714154], [0.24531644877001227, 0.7546835512299876], [0.0023969272650335633, 0.9976030727349664], [0.0037476747393794183, 0.9962523252606206], [0.05130808639847684, 0.9486919136015232], [0.008217753143655216, 0.9917822468563448], [0.025393369426237336, 0.9746066305737627], [0.1385368435601203, 0.8614631564398797], [0.4815579929828904, 0.5184420070171096]]

print("Correlated Bayes Factor")
bayes_factor(correlated_y, x)
print("------------------------------------")
print("Random Bayes Factor")
bayes_factor(random_y, x)
"""


"""
res = minimize(objective_r_numer, input_value, args=args, method='Nelder-Mead', tol=1e-12)
print(res)

#print(objective_r_numer([4.51850662,  7.99117481], *args))

H = nd.Hessian(objective_r_numer)(res.x, *args)
print(H)

def result_of_laplace_method():
	return (2 * math.pi / n)**(len(input_value)/2.0) * math.sqrt(1.0 / np.linalg.det(H)) * math.exp(-1.0 * n * objective_r_numer(res.x, *args)) 

print("Result from Laplace method is: ")
print(result_of_laplace_method())
print("---------------------------------------")


def f(a, b):
	mu = [a, b]   # Dirichlete pdf Parameters (exp)
	eta = [math.exp(a), math.exp(b)]  # Dirichlete pdf Parameters (eta in our paper)
	y = [y_1, y_2]    # isoform proportions

	K = 2
	variance = x**2 * sigma_0**2 + sigma_1**2

	beta_func = mb.multivariateBeta(eta)
	y_exponents = 1.0
	exponentials = 1.0
	for i in range(0, K):
		y_exponents = y_exponents * (y[i]**(eta[i] - 1))
		exponentials = exponentials * math.exp(-(mu[i] - beta_0k[i])**2/(2*(variance)))
	other_factors = 1.0/(math.sqrt(2 * math.pi))**K * 1.0/((variance)**K)

	return beta_func * y_exponents * exponentials * other_factors

result = integrate.dblquad(f, -20, 6, lambda x: -20, lambda x: 6)

print("from integration: ", result)
"""

