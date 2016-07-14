import numpy as np
import scipy.integrate as integrate
import scipy.special as special
import math

y_1 = 0.02
y_2 = 0.98
sigma = 1

def integrand(y, x):
	gamma_top = special.gamma(y + x)
	gamma_bot = special.gamma(x) * special.gamma(y)
	first_y = y_1 ** (x - 1)
	sec_y = y_2 ** (y - 1)
	exp_1 = math.exp(-(math.log(x))**2/(2*(sigma**2)))
	exp_2 = math.exp(-(math.log(y))**2/(2*(sigma**2)))

	return gamma_top * first_y * sec_y * exp_1 * exp_2 / (2 * math.pi * gamma_bot * (sigma**2) * x * y)

result = integrate.dblquad(integrand, 0, 80, lambda x: 0, lambda x: 80)

print(result)