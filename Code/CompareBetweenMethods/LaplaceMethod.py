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
n = 100000.0 # for use in Laplace approximation


# variables for testing purposes
y_1 = 0.1
y_2 = 1- y_1
x = 2.0
sigma_0 = 1.0
sigma_1 = 1.0

input_value = [0.06076677,  0.47356138, 1.0, 1.0] #first K mu, which is log(eta), next K beta_0k
args = (y_1, y_2, x, sigma_0, sigma_1)                       #isoform proportions

def objective_r_denom(input_value, *args):
	mu = [] # Dirichlete pdf Parameters (exp)
	eta = [] # Dirichlete pdf Parameters (eta in our paper)
	y = [] # isoform proportions

	for v in input_value:
		mu.append(v)
		eta.append(math.exp(v))
	for v in args:
		y.append(v)
		
	beta_func = mb.multivariateBeta(eta)
	y_exponents = 1
	exponentials = 1
	for i in range(0, len(args)):
		y_exponents = y_exponents * (y[i]**(eta[i] - 1))
		exponentials = exponentials * math.exp(-(mu[i]**2)/(2*(sigma**2)))
	other_factors = 1.0/(math.sqrt(2 * math.pi) * sigma)**(len(args))

	return -1.0 * (1.0/n) * math.log(beta_func * y_exponents * exponentials * other_factors)

def objective_r_numer(input_value, *args):
	mu = []   # Dirichlete pdf Parameters (exp)
	eta = []  # Dirichlete pdf Parameters (eta in our paper)
	beta_0k = [] # Beta_0k
	y = []    # isoform proportions

	K = int(len(input_value)/2.0)
	for k in range(0, K):
		mu.append(input_value[k])
		eta.append(math.exp(input_value[k]))
	for k in range(K, 2*K):
		beta_0k.append(input_value[k])
	for k in range(0, K):
		y.append(args[k])
	x = args[K]
	sigma_0 = args[K+1]
	sigma_1 = args[K+2]
	variance = x**2 * sigma_0**2 + sigma_1**2

	beta_func = mb.multivariateBeta(eta)
	y_exponents = 1
	exponentials = 1
	for i in range(0, K):
		y_exponents = y_exponents * (y[i]**(eta[i] - 1))
		exponentials = exponentials * math.exp(-(mu[i] - beta_0k[i])**2/(2*variance))
	other_factors = 1.0/(math.sqrt(2 * math.pi) * variance)**K * 1.0/(variance**K)

	return -1.0 * (1.0/n) * math.log(beta_func * y_exponents * exponentials * other_factors)
"""
res = minimize(objective_r_numer, input_value, args=args, method='Nelder-Mead', tol=1e-12)
print(res)

H = nd.Hessian(objective_r_numer)(res.x, *args)
print(H)

def result_of_laplace_method():
	return (2 * math.pi / n)**(len(args)/2.0) * math.sqrt(1.0 / np.linalg.det(H)) * math.exp(-1.0 * n * objective_r_numer(res.x, *args)) 

print("Result from Laplace method is: ")
print(result_of_laplace_method())
print("---------------------------------------")


"""

def f(a, b, c, d):
	mu = [a, b]   # Dirichlete pdf Parameters (exp)
	eta = [math.exp(a), math.exp(b)]  # Dirichlete pdf Parameters (eta in our paper)
	beta_0k = [c, d] # Beta_0k
	y = [y_1, y_2]    # isoform proportions

	K = 2

	variance = x**2 * sigma_0**2 + sigma_1**2

	beta_func = mb.multivariateBeta(eta)
	y_exponents = 1
	exponentials = 1
	for i in range(0, K):
		y_exponents = y_exponents * (y[i]**(eta[i] - 1))
		exponentials = exponentials * math.exp(-(mu[i] - beta_0k[i])**2/(2*variance))
	other_factors = 1.0/(math.sqrt(2 * math.pi) * variance)**K * 1.0/(variance**K)

	return -1.0 * (1.0/n) * math.log(beta_func * y_exponents * exponentials * other_factors)

result = integrate.nquad(f, [[-5,5],[-5,5],[-5,5],[-5,5]])

print(result)