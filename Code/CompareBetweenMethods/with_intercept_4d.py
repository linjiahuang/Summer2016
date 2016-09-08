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
sigma = 1.0
n = 100000.0

# for use with Python's nquad method (without manually integrating beta_0k)
def f(beta_02, beta_01, mu_22, mu_21, mu_12, mu_11):
	beta_input1 = [math.exp(mu_11), math.exp(mu_12)]
	beta_input2 = [math.exp(mu_21), math.exp(mu_22)]

	beta_func1 = mb.multivariateBeta(beta_input1)
	beta_func2 = mb.multivariateBeta(beta_input2)

	exp_y_11 = y_11 ** (math.exp(mu_11) - 1.0)
	exp_y_12 = y_12 ** (math.exp(mu_12) - 1.0)
	exp_y_21 = y_21 ** (math.exp(mu_21) - 1.0)
	exp_y_22 = y_22 ** (math.exp(mu_22) - 1.0)

	exp_11 = math.exp(-(mu_11 - beta_01)**2 / (2.0 * sigma**2))
	exp_12 = math.exp(-(mu_12 - beta_02)**2 / (2.0 * sigma**2))
	exp_21 = math.exp(-(mu_21 - beta_01)**2 / (2.0 * sigma**2))
	exp_22 = math.exp(-(mu_22 - beta_02)**2 / (2.0 * sigma**2))

	return beta_func1 * beta_func2 * exp_y_11 * exp_y_12 * exp_y_21 * exp_y_22 * exp_11 * exp_12 * exp_21 * exp_22 / (4 * math.pi * sigma**4)
#result = integrate.nquad(f, [[1, 2],[1, 2],[1, 2],[1, 2],[1, 2],[1, 2]])
#print(result)


# for use with Python's nquad method (manually integrated beta_0k)
def integrand(mu_22, mu_21, mu_12, mu_11):
	beta_input1 = [math.exp(mu_11), math.exp(mu_12)]
	beta_input2 = [math.exp(mu_21), math.exp(mu_22)]

	beta_func1 = mb.multivariateBeta(beta_input1)
	beta_func2 = mb.multivariateBeta(beta_input2)

	exp_y_11 = y_11 ** (math.exp(mu_11) - 1.0)
	exp_y_12 = y_12 ** (math.exp(mu_12) - 1.0)
	exp_y_21 = y_21 ** (math.exp(mu_21) - 1.0)
	exp_y_22 = y_22 ** (math.exp(mu_22) - 1.0)

	exp_1 = math.exp(-((mu_11 - mu_21)/2.0)**2)
	exp_2 = math.exp(-((mu_12 - mu_22)/2.0)**2)

	return beta_func1 * beta_func2 * exp_y_11 * exp_y_12 * exp_y_21 * exp_y_22 * exp_1 * exp_2 / (4 * math.pi * sigma**4)
#result = integrate.nquad(integrand, [[-5, 3],[-5, 3],[-5, 3],[-5, 3]])
#print(result)

"""
num_trial = 100000

for i in range(0, num_trial):
	result = 0
	mu = np.random.uniform(-5, 5, 4)
	result = result + integrand(mu[0], mu[1], mu[2], mu[3])

print(result/num_trial)
"""

# for use with Laplace method (without manually integrating beta_0k)
def objective_rr(input_value, *args):

	beta_02, beta_01, mu_22, mu_21, mu_12, mu_11 = input_value

	beta_input1 = [math.exp(mu_11), math.exp(mu_12)]
	beta_input2 = [math.exp(mu_21), math.exp(mu_22)]

	beta_func1 = mb.multivariateBeta(beta_input1)
	beta_func2 = mb.multivariateBeta(beta_input2)

	exp_y_11 = y_11 ** (math.exp(mu_11) - 1.0)
	exp_y_12 = y_12 ** (math.exp(mu_12) - 1.0)
	exp_y_21 = y_21 ** (math.exp(mu_21) - 1.0)
	exp_y_22 = y_22 ** (math.exp(mu_22) - 1.0)

	exp_11 = math.exp(-(mu_11 - beta_01)**2 / (2.0 * sigma**2))
	exp_12 = math.exp(-(mu_12 - beta_02)**2 / (2.0 * sigma**2))
	exp_21 = math.exp(-(mu_21 - beta_01)**2 / (2.0 * sigma**2))
	exp_22 = math.exp(-(mu_22 - beta_02)**2 / (2.0 * sigma**2))

	return -1.0 * (1.0/n) * math.log(beta_func1 * beta_func2 * exp_y_11 * exp_y_12 * exp_y_21 * exp_y_22 * exp_11 * exp_12 * exp_21 * exp_22/ (4 * math.pi * sigma**4))

#res = minimize(objective_rr, input_value, args=args, method='Nelder-Mead', tol=1e-12)
#print(res)


# for use with Laplace method (manually integrated beta_0k)
def objective_r(input_value, *args):

	mu_22, mu_21, mu_12, mu_11 = input_value
	"""
	beta_input1 = [math.exp(mu_11), math.exp(mu_12)]
	beta_input2 = [math.exp(mu_21), math.exp(mu_22)]
	
	beta_func1 = 1.0/(special.beta(math.exp(mu_11), math.exp(mu_12)))
	beta_func2 = 1.0/(special.beta(math.exp(mu_21), math.exp(mu_22)))

	exp_y_11 = y_11 ** (math.exp(mu_11) - 1.0)
	exp_y_12 = y_12 ** (math.exp(mu_12) - 1.0)
	exp_y_21 = y_21 ** (math.exp(mu_21) - 1.0)
	exp_y_22 = y_22 ** (math.exp(mu_22) - 1.0)
	
	exp_11 = math.exp(-mu_11**2 / 2.0)
	exp_12 = math.exp(-mu_21**2 / 2.0)
	exp_21 = math.exp(-mu_12**2 / 2.0)
	exp_22 = math.exp(-mu_22**2 / 2.0)

	exp_1 = math.exp(((mu_11 + mu_21)/2.0)**2)
	exp_2 = math.exp(((mu_12 + mu_22)/2.0)**2)
	"""

	log_beta_func1 = -special.betaln(math.exp(mu_11), math.exp(mu_12))
	log_beta_func2 = -special.betaln(math.exp(mu_21), math.exp(mu_22))

	log_exp_y_11 = math.log(y_11) * (math.exp(mu_11) - 1.0)
	log_exp_y_12 = math.log(y_12) * (math.exp(mu_12) - 1.0)
	log_exp_y_21 = math.log(y_21) * (math.exp(mu_21) - 1.0)
	log_exp_y_22 = math.log(y_22) * (math.exp(mu_22) - 1.0)

	log_exp_1 = -((mu_11 - mu_21)/2.0)**2
	log_exp_2 = -((mu_12 - mu_22)/2.0)**2

	return -1.0 * (1.0/n) * (log_beta_func1 + log_beta_func2 + log_exp_y_11 + log_exp_y_12 + log_exp_y_21 + log_exp_y_22 + log_exp_1 + log_exp_2 - math.log(4 * math.pi * sigma**4))


input_value = [5, 5, 5, 5]
args = (n, 1.1)

#print(objective_r(input_value, *args))

#res = minimize(objective_r, input_value, args=args, method='powell', tol=1e-12)
#print(res)
#H = nd.Hessian(objective_r)(res.x, *args)
#print(np.linalg.det(H))

def result_of_laplace_method(input_value, *args):
	mu_22, mu_21, mu_12, mu_11 = input_value

	return (2 * math.pi / n)**2 * math.sqrt(1.0 / np.linalg.det(H)) * math.exp(-1.0 * n * objective_r(res.x, *args)) 

#print(result_of_laplace_method(res.x, *args))

def f(input_value, *args):
	x, y = input_value
	a, b = args
	first_term = 1.0/special.beta(x, y)
	second_term = a**(x-1) * b**(y-1)
	return first_term * second_term

input_value = [4, 5]
args = (0.4, 0.6)

res = minimize(f, input_value, args=args, method='Nelder-Mead', tol=1e-12)
print(res)

