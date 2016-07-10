import scipy.special as special
import numpy as np
from scipy import optimize
import math

args = (0.4, 0.6)

def f(x, *args):
	alpha, beta = x
	y_1, y_2 = args
	first_term = special.gamma(alpha + beta) * y_1**(alpha-1) * y_2**(beta-1) / (special.gamma(alpha) * special.gamma(beta))
	second_term = (1/alpha)*math.exp(-1*(math.log(alpha))**2)*(1/beta)*math.exp(-1*(math.log(beta))**2)
	return -1 * math.log(first_term * second_term)

def samef(alpha, beta):
	y_1 = 0.4
	y_2 = 0.6
	first_term = special.gamma(alpha + beta) * y_1**(alpha-1) * y_2**(beta-1) / (special.gamma(alpha) * special.gamma(beta))
	second_term = (1/alpha)*math.exp(-1*(math.log(alpha))**2)*(1/beta)*math.exp(-1*(math.log(beta))**2)
	return -1 * math.log(first_term * second_term)

def gradfalpha(alpha, beta, y_1):
	return special.digamma(alpha + beta)/special.gamma(alpha + beta) + math.log(y_1) - special.digamma(alpha)/special.gamma(alpha) - 1/alpha - 2*math.log(alpha)/alpha

def gradfbeta(alpha, beta, y_2):
	return special.digamma(alpha + beta)/special.gamma(alpha + beta) + math.log(y_2) - special.digamma(beta)/special.gamma(beta) - 1/beta - 2*math.log(beta)/beta

x0 = np.asarray((0.6, 0.7))

res1 = optimize.fmin_cg(f, x0, args=args)
#t1 , t2 = res1

print(res1)
print(args)
print(samef(res1[0], res1[1]))
print("derivative: " + str(gradfalpha(res1[0], res1[1], 0.4)))
print("derivative: " + str(gradfbeta(res1[0], res1[1], 0.6)))




