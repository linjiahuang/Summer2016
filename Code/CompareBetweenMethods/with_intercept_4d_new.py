import sys
sys.path.insert(0, '/home/jiahuang/Dropbox/Princeton/Sophomore Summer/Summer Project/Summer2016/Code')

import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize
import scipy.special as special
import MultivariateBeta as mb
import math
import numdifftools as nd

y_11 = 0.5
y_12 = 1.0 - y_11
y_21 = 0.5
y_22 = 1.0 - y_21
sigma_0 = 1.0
sigma_1 = 1.0

# for use with Python's nquad method (without manually integrating beta_0k)
def f(beta_2, beta_1, beta_02, beta_01, mu_22, mu_21, mu_12, mu_11):

	# sets the values for x
	x_1 = 1.0
	x_2 = 1.0

	beta_input1 = [math.exp(mu_11), math.exp(mu_12)]
	beta_input2 = [math.exp(mu_21), math.exp(mu_22)]

	beta_func1 = mb.multivariateBeta(beta_input1)
	beta_func2 = mb.multivariateBeta(beta_input2)

	exp_y_11 = y_11 ** (math.exp(mu_11) - 1.0)
	exp_y_12 = y_12 ** (math.exp(mu_12) - 1.0)
	exp_y_21 = y_21 ** (math.exp(mu_21) - 1.0)
	exp_y_22 = y_22 ** (math.exp(mu_22) - 1.0)

	exp_11 = math.exp(-(mu_11 - beta_01 - x_1*beta_1)**2 / (2.0 * sigma_1**2))
	exp_12 = math.exp(-(mu_12 - beta_02 - x_1*beta_2)**2 / (2.0 * sigma_1**2))
	exp_21 = math.exp(-(mu_21 - beta_01 - x_2*beta_1)**2 / (2.0 * sigma_1**2))
	exp_22 = math.exp(-(mu_22 - beta_02 - x_2*beta_2)**2 / (2.0 * sigma_1**2))

	exp_beta_1 = math.exp(-(beta_1)**2 / (2.0 * sigma_0**2))
	exp_beta_2 = math.exp(-(beta_2)**2 / (2.0 * sigma_0**2))

	constant_factor = (1.0/(4 * math.pi**2 * sigma_1**4)) * (1.0/(2 * math.pi * sigma_0**2))

	return constant_factor * beta_func1 * beta_func2 * exp_y_11 * exp_y_12 * exp_y_21 * exp_y_22 * exp_11 * exp_12 * exp_21 * exp_22 * exp_beta_1 * exp_beta_2

#result = integrate.nquad(f, [[1, 2],[1, 2],[1, 2],[1, 2],[1, 2],[1, 2],[1, 2],[1, 2]])
#print(result)
#Ran for about 12 minutes with no result
#Then ran with Mathematica and got some warning. 

# for use with Python's nquad method (manually integrated beta_0k)
def f(beta_2, beta_1, mu_22, mu_21, mu_12, mu_11):

	# sets the values for x
	x_1 = 1.0
	x_2 = 1.0

	beta_input1 = [math.exp(mu_11), math.exp(mu_12)]
	beta_input2 = [math.exp(mu_21), math.exp(mu_22)]

	beta_func1 = mb.multivariateBeta(beta_input1)
	beta_func2 = mb.multivariateBeta(beta_input2)

	exp_y_11 = y_11 ** (math.exp(mu_11) - 1.0)
	exp_y_12 = y_12 ** (math.exp(mu_12) - 1.0)
	exp_y_21 = y_21 ** (math.exp(mu_21) - 1.0)
	exp_y_22 = y_22 ** (math.exp(mu_22) - 1.0)

	exp_11 = math.exp(-(mu_11 - x_1*beta_1)**2 / (2.0 * sigma_1**2))
	exp_12 = math.exp(-(mu_12 - x_1*beta_2)**2 / (2.0 * sigma_1**2))
	exp_21 = math.exp(-(mu_21 - x_2*beta_1)**2 / (2.0 * sigma_1**2))
	exp_22 = math.exp(-(mu_22 - x_2*beta_2)**2 / (2.0 * sigma_1**2))
	exp_1 = math.exp((mu_11 - x_1*beta_1 + mu_21 - x_2*beta_1)**2 / (4.0 * sigma_1**2))
	exp_2 = math.exp((mu_12 - x_1*beta_2 + mu_22 - x_2*beta_2)**2 / (4.0 * sigma_1**2))

	exp_beta_1 = math.exp(-(beta_1)**2 / (2.0 * sigma_0**2))
	exp_beta_2 = math.exp(-(beta_2)**2 / (2.0 * sigma_0**2))

	constant_factor = (1.0/(4 * math.pi * sigma_1**2)) * (1.0/(2 * math.pi * sigma_0**2))

	return constant_factor * beta_func1 * beta_func2 * exp_y_11 * exp_y_12 * exp_y_21 * exp_y_22 * exp_11 * exp_12 * exp_21 * exp_22 * exp_1 * exp_2 * exp_beta_1 * exp_beta_2

#result = integrate.nquad(f, [[1, 2],[1, 2],[1, 2],[1, 2],[1, 2],[1, 2]])
#print(result)
#Ran for about 20 minutes and no result




