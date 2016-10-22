import sys
sys.path.insert(0, '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/Code')

import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize
import scipy.special as special
import MultivariateBeta as mb
import math
import numdifftools as nd

y_11 = 0.1
y_12 = 1.0 - y_11
y_21 = 0.2
y_22 = 1.0 - y_21
sigma = 1.0

def func_g(mu_11, mu_12, mu_21, mu_22):
	beta_input_1 = [math.exp(mu_11), math.exp(mu_12)]
	beta_input_2 = [math.exp(mu_21), math.exp(mu_22)]

	beta_func_1 = mb.multivariateBeta(beta_input_1)
	beta_func_2 = mb.multivariateBeta(beta_input_2)

	first_y = y_11 ** (math.exp(mu_11) - 1.0)
	sec_y = y_12 ** (math.exp(mu_12) - 1.0)
	third_y = y_21 ** (math.exp(mu_21) - 1.0)
	fourth_y = y_22 ** (math.exp(mu_22) - 1.0)

	return beta_func_1 * beta_func_2 * first_y * sec_y * third_y * fourth_y

def mc_given_beta(beta_01, beta_02, total_iterations):
	summation = 0
	for i in range(0, total_iterations):
		mu_11 = np.random.normal(beta_01, sigma)
		mu_12 = np.random.normal(beta_02, sigma)
		mu_21 = np.random.normal(beta_01, sigma)
		mu_22 = np.random.normal(beta_02, sigma)
		summation = summation + float(func_g(mu_11, mu_12, mu_21, mu_22))
	return summation/total_iterations

def mc_total(lower_bound_for_beta, upper_bound_for_beta, total_iterations):
	summation = 0
	for i in range(0, total_iterations):
		beta_01 = np.random.uniform(lower_bound_for_beta, upper_bound_for_beta)
		beta_02 = np.random.uniform(lower_bound_for_beta, upper_bound_for_beta)
		summation = summation + mc_given_beta(beta_01, beta_02, total_iterations)
	return summation/total_iterations * (upper_bound_for_beta - lower_bound_for_beta)**2

for i in range(0, 3):
	print("[-1, 1]: " + str(mc_total(-1, 1, 1000)))

for i in range(0, 3):
	print("[-2, 2]: " + str(mc_total(-2, 2, 1000)))

for i in range(0, 3):
	print("[-3, 3]: " + str(mc_total(-3, 3, 1000)))

for i in range(0, 3):
	print("[-4, 4]: " + str(mc_total(-4, 4, 1000)))
