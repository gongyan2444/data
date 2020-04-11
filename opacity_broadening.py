import numpy as np
import scipy.optimize

"""
the script is used to show how to derive opacity from line width ratios 
This is also a good example to show how to solve non-linear equations with Python
"""
def G(x,r):
	return 1./np.sqrt(np.log(2.))*(np.log(x/np.log(2./(np.exp(-x)+1))))**(0.5)-r 


def F(x,r):
	return np.log(x/np.log(2./(np.exp(-x)+1)))-r**2*np.log(2) 

def F2(x):
	return np.exp(x)+x**2.-10


def F3(x):
	return np.log(x/np.log(2./(np.exp(-x)+1)))-2**2*np.log(2) 
# test the function
print(G(11.09011069,2))  ### should be close to 0 

# two examples to solve a non-linear equation
x = scipy.optimize.broyden1(F2, [2], f_tol=1e-14)
print(x)

y = scipy.optimize.broyden1(F3, [2], f_tol=1e-14)
print(y)



# This provides an example to solve a non-linear equation in a loop

for i in np.arange(11):
	def F4(x):
		return np.log(x/np.log(2./(np.exp(-x)+1)))-(i*0.1+1)**2*np.log(2) 
	z = scipy.optimize.broyden1(F4, [2], f_tol=1e-14)
	print("results for a ratio of ", i*0.1+1)
	print(z)




