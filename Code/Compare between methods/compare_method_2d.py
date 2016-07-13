import sys
sys.path.insert(0, '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/Code')

import numpy as np
import scipy.integrate as integrate
import scipy.special as special
import MultivariateBeta as mb
import math

y_1 = 0.58
y_2 = 0.42
sigma = 1

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

#result = integrate.dblquad(integrand, -6, 6, lambda x: -6, lambda x: 6)

#print("From integration     : " + str(result))

def mc(n):
	summation = 0
	for i in range(0, n):
		x = np.random.normal(0, sigma)
		y = np.random.normal(0, sigma)
		summation = summation + func_g(y, x)
	return summation/n

#print("Mean from MC         : " + str(mc(100000)))




