import numpy as np
import scipy.special as special

N = 10000
sumOfTerms = 0
y_1 = 0.1
y_2 = 1 - y_1
eta_2 = 7
sigma = 1

def term(eta_1, eta_2):
	return special.gamma(eta_1+eta_2)/(special.gamma(eta_1)*special.gamma(eta_2)) * (y_1**(eta_1-1)) * (y_2**(eta_2-1))

for i in range(0, N):
	eta_1 = np.random.lognormal(0,sigma)
	eta_2 = np.random.lognormal(0,sigma)
	sumOfTerms = sumOfTerms + term(eta_1, eta_2)

print(sumOfTerms/N)

