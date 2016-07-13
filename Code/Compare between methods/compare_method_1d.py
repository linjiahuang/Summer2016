"""Tested a simple integral (defined in integrand) using direct integration, MC, MCMC and the Laplace method.
   Note that this integral is in 1 dimension.
"""

import numpy as np
import scipy.integrate as integrate
import math

sigma = 1

def func_f(x):
	return math.exp(-x**2/(2 * sigma**2)) / (math.sqrt(2 * math.pi) * sigma)

def integrand(x):
	return x**2 * math.exp(-x**2/(2 * sigma**2)) / (math.sqrt(2 * math.pi) * sigma)

result = integrate.quad(integrand, -np.inf, np.inf)

def mc(n):
	summation = 0
	for i in range (0, n):
		x = np.random.normal(0, sigma)
		summation = summation + x**2
	return summation/n

# MCMC global variables
b = 1

def mcmc(n):
	summation = 0
	x_old = np.random.normal(0, b)
	for i in range(0, n):
		summation = summation + x_old**2
		x_new = np.random.normal(x_old, b)
		r = min(func_f(x_new)/func_f(x_old), 1)
		u = np.random.uniform()
		if (u < r):
			x_old = x_new
	return summation/n

print("From integration: " + str(result))

def get_average(inputList):
	return sum(inputList) / len(inputList)

def laplace_method(n):
	return math.sqrt(2 * math.pi / n) * math.sqrt((n * sigma**2) / 2) * math.exp(-1 * (1/n) * math.log(2 * sigma / (math.exp(1) * math.sqrt(2 * math.pi))))

def compare_methods(m, n):
	mc_list = []
	mcmc_list = []
	for i in range(0, m):
		mc_list.append(mc(n))
		mcmc_list.append(mcmc(n))
	mc_mean = get_average(mc_list)
	mcmc_mean = get_average(mcmc_list)
	mc_var = 0
	mcmc_var = 0
	for i in range(0, m):
		mc_var = mc_var + (mc_list[i] - mc_mean)**2
		mcmc_var = mcmc_var + (mcmc_list[i] - mcmc_mean)**2
	print("Mean from MC         : " + str(get_average(mc_list)))
	print("Mean from MCMC       : " + str(get_average(mcmc_list)))
	print("Var from MC          : " + str(mc_var/m))
	print("Var from MCMC        : " + str(mcmc_var/m))
	print("Ans from Laplace     : " + str(laplace_method(1000)))

compare_methods(100, 10000)

