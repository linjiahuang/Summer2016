"""Implemented the Laplace approximation in integrating a 2D function. 
   Contains both the Hessian from numdifftools and from closed form solution
"""

import sys
sys.path.insert(0, '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/Code')

import numpy as np
from scipy.optimize import minimize
import scipy.special as special
import MultivariateBeta as mb
import math
import numdifftools as nd

input_value = [0.06076677,  0.47356138]
args = (100000.0, 0.0001, 0.9999, 1.0)

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

def grad1_objective_r(input_value, *args):
	mu_1, mu_2 = input_value
	n, y_1, y_2, sigma = args

	x = math.exp(mu_1)
	y = math.exp(mu_2)

	return -1.0 * (1.0/n) * (special.digamma(x + y)*x + math.log(y_1)*x - mu_1/sigma**2 - special.digamma(x)*x)

def grad2_objective_r(input_value, *args):
	mu_1, mu_2 = input_value
	n, y_1, y_2, sigma = args

	x = math.exp(mu_1)
	y = math.exp(mu_2)

	return -1.0 * (1.0/n) * (special.digamma(x + y)*y + math.log(y_2)*y - mu_2/sigma**2 - special.digamma(y)*y)
"""
grad = nd.Gradient(objective_r)(input_value, *args)
print("-------------------------------------------")
print("numdifftools gradient is: ")
print(grad)
print("--------------------------")
print("analytical gradient is: ")
print(grad1_objective_r(input_value, *args))
print("--------------------------")
"""
res = minimize(objective_r, input_value, args=args, method='Nelder-Mead', tol=1e-12)
"""
print(res)
print("--------------------------")
func = objective_r(res.x, *args)
print("analytical function is: " + str(func))
print("--------------------------")
"""
def matrix_one_one(input_value, *args):
	mu_1, mu_2 = input_value
	n, y_1, y_2, sigma = args

	x = math.exp(mu_1)
	y = math.exp(mu_2)

	first_gamma = special.polygamma(1, x + y) * x * x + special.digamma(x + y) * x  
	mid_terms = x * math.log(y_1) - 1.0/sigma**2
	second_gamma = special.polygamma(1, x) * x * x + special.digamma(x) * x 

	return -1.0 * (1.0/n) * (first_gamma + mid_terms - second_gamma)

def matrix_two_two(input_value, *args):
	mu_1, mu_2 = input_value
	n, y_1, y_2, sigma = args

	x = math.exp(mu_1)
	y = math.exp(mu_2)

	first_gamma = special.polygamma(1, x + y) * y * y + special.digamma(x + y) * y
	mid_terms = y * math.log(y_2) - 1.0/sigma**2
	second_gamma = special.polygamma(1, y) * y * y + special.digamma(y) * y

	return -1.0 * (1.0/n) * (first_gamma + mid_terms - second_gamma)

def matrix_other_entry(input_value, *args):
	mu_1, mu_2 = input_value
	n, y_1, y_2, sigma = args

	x = math.exp(mu_1)
	y = math.exp(mu_2)

	return -1.0 * (1.0/n) *(special.polygamma(1, x + y) * x * y)
"""
print("--------------------------")
print("first entry of matrix is : " + str(matrix_one_one(res.x, *args)))

print("second entry of matrix is: " + str(matrix_two_two(res.x, *args)))

print("other entry of matrix is : " + str(matrix_other_entry(res.x, *args)))

H = nd.Hessian(objective_r)(res.x, *args)
print(H)
print("--------------------------")
"""




"""
matrix_r = np.array([[matrix_one_one(res.x, *args), matrix_other_entry(res.x, *args)], [matrix_other_entry(res.x, *args), matrix_two_two(res.x, *args)]])

def result_of_laplace_method(input_value, *args):
	mu_1, mu_2 = input_value
	n, y_1, y_2, sigma = args

	return (2 * math.pi / n) * math.sqrt(1.0 / np.linalg.det(matrix_r)) * math.exp(-1.0 * n * objective_r(res.x, *args)) 
print("Result from Laplace method is: ")
print(result_of_laplace_method(res.x, *args))
print("---------------------------------------")
"""




"""
def f(input_value, *args):
	x, y = input_value
	a, b = args

	return x**2 + y**2

input_value = [1, 2]
args = (3, 4)

H = nd.Hessian(f)(input_value, *args)
print(H)

grad = nd.Gradient(f)(input_value, *args)
print(grad)
"""