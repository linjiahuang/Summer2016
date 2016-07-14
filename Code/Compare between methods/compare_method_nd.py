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
sigma = 1.0
y_1 = 0.001
y_2 = 0.234
y_3 = 0.676
y_4 = 1.0 - y_1 - y_2 - y_3

input_value = [0.06076677,  0.47356138, 0.341, 0.111]
args = (y_1, y_2, y_3, y_4)

def objective_r(input_value, *args):
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

res = minimize(objective_r, input_value, args=args, method='Nelder-Mead', tol=1e-12)
print(res)

H = nd.Hessian(objective_r)(res.x, *args)
print(H)

def result_of_laplace_method():
	return (2 * math.pi / n)**(len(args)/2.0) * math.sqrt(1.0 / np.linalg.det(H)) * math.exp(-1.0 * n * objective_r(res.x, *args)) 

print("Result from Laplace method is: ")
print(result_of_laplace_method())
print("---------------------------------------")