import numpy as np
import scipy.integrate as integrate
from scipy.optimize import minimize
import scipy.special as special
import math
import numdifftools as nd

y_1 = 0.5
y_2 = 1.0 - y_1

def f(b, a):

	logbeta = -special.betaln(math.exp(a), math.exp(b))
	other_logbeta = (math.exp(a) + math.exp(b))* math.log(math.exp(a) + math.exp(b)) - math.exp(a)*a - math.exp(b)*b
	#print("log is: ", logbeta)
	#print("other is: ", other_logbeta)
	log_y_1 = (math.exp(a) - 1.0) * math.log(y_1)
	#print(log_y_1)
	log_y_2 = (math.exp(b) - 1.0) * math.log(y_2)
	#print(log_y_2)
	#print(logbeta + log_y_1 + log_y_2 )

	return math.exp(a)*(math.log(math.exp(a) + math.exp(b)) - a + math.log(y_1)) + math.exp(b)*(math.log(math.exp(a) + math.exp(b)) - b + math.log(y_2))

#result = integrate.dblquad(f, -np.inf,6.235, lambda x: -np.inf, lambda x: 6.235)
#print(result)

for i in range(2, 200, 20):
	print(f(i, i))
