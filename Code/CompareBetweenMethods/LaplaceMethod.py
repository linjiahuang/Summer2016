import sys
sys.path.insert(0, '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/Code')

import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize
import scipy.special as special
import MultivariateBeta as mb
import math
import numdifftools as nd

# Global constants
n = 10000000.0 # for use in Laplace approximation


# variables for testing purposes
y_1 = 0.03
y_2 = 1.0 - y_1
x = 2.0
sigma_0 = 1.0
sigma_1 = 1.0
beta_0k = [-6.0, 8]

input_value = [0.06076677,  0.47356138] #first K mu, which is log(eta), next K beta_0k
args = ([y_1, y_2], x)                       #isoform proportions

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
	"""
	print("beta_func is ", beta_func)
	print("y_exponents is ", log_y_exponents)
	print("exponentials is", log_exponentials)
	print("other_factors is", other_factors)
	"""
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
