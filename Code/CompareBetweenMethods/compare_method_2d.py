import sys
sys.path.insert(0, '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/Code')

import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize
import scipy.special as special
import MultivariateBeta as mb
import math
import numdifftools as nd

y_1 = 0.1
y_2 = 1.0 - y_1
sigma = 1.0

input_value = [0.06076677,  0.47356138] # starting value for optimization
args = (100000.0, y_1, y_2, sigma)

def integrand(y, x):
	beta_input = [math.exp(x), math.exp(y)]
	beta_func = mb.multivariateBeta(beta_input)
	first_y = y_1 ** (math.exp(x) - 1.0)
	sec_y = y_2 ** (math.exp(y) - 1.0)
	exp_1 = math.exp(-(x**2)/(2*(sigma**2)))
	exp_2 = math.exp(-(y**2)/(2*(sigma**2)))

	return beta_func * first_y * sec_y * exp_1 * exp_2 / (2 * math.pi * (sigma**2))

def func_g(y, x):
	beta_input = [math.exp(x), math.exp(y)]
	beta_func = mb.multivariateBeta(beta_input)
	first_y = y_1 ** (math.exp(x) - 1.0)
	sec_y = y_2 ** (math.exp(y) - 1.0)

	return beta_func * first_y * sec_y

result = integrate.dblquad(integrand, -6, 6, lambda x: -6, lambda x: 6)

print("From integration     : " + str(result))

def mc(n):
	summation = 0
	for i in range(0, n):
		x = np.random.normal(0, sigma)
		y = np.random.normal(0, sigma)
		summation = summation + func_g(y, x)
	return summation/n

print("Mean from MC         : " + str(mc(100000)))

def objective_r(input_value, *args):
	mu_1, mu_2 = input_value
	n, y_1, y_2, sigma = args

	x = math.exp(mu_1)
	y = math.exp(mu_2)

	beta_input = [x, y]
	beta_func = mb.multivariateBeta(beta_input)
	first_y = y_1 ** (x - 1.0)
	sec_y = y_2 ** (y - 1.0)
	exp_1 = math.exp(-(mu_1**2)/(2*(sigma**2)))
	exp_2 = math.exp(-(mu_2**2)/(2*(sigma**2)))

	return -1.0 * (1.0/n) * math.log(beta_func * first_y * sec_y * exp_1 * exp_2 / (2 * math.pi * sigma**2))

res = minimize(objective_r, input_value, args=args, method='Nelder-Mead', tol=1e-12)
H = nd.Hessian(objective_r)(res.x, *args)

def result_of_laplace_method(input_value, *args):
	mu_1, mu_2 = input_value
	n, y_1, y_2, sigma = args

	return (2 * math.pi / n) * math.sqrt(1.0 / np.linalg.det(H)) * math.exp(-1.0 * n * objective_r(res.x, *args)) 
print("Result from Laplace method is: ")
print(result_of_laplace_method(res.x, *args))
print("---------------------------------------")



