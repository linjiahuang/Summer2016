import numpy as np
import scipy.integrate as integrate
import scipy.special as special
import math

result = integrate.quad(lambda x: (special.gamma(x+7)/special.gamma(x))*(0.1**(x-1))*(1/x)*((1/math.sqrt(2*math.pi))*math.exp(-((math.log(x))**2)/2)), 0, 10000)
print(result)

# 