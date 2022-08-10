import numpy as np
import matplotlib.pyplot as plt 

# reference: https://ui.adsabs.harvard.edu/abs/2001ApJ...561..254T/abstract
# see fig. 1.
# kinetic temperature 
#tt = np.arange(300)+5

def cal_opr_h2(tt):
	rotb =87.6
	s1 = 0
	s2 = 0 
	## use 100 as a large number
	for i in np.arange(100):
		j1 = 2*i+1
		j2 = 2*i
		tmp1 = rotb/tt*(j1*(j1+1))
		tmp2 = rotb/tt*(j2*(j2+1))
		s1 = s1 + (2*j1+1)*np.exp(-tmp1)
		s2 = s2 + (2*j2+1)*np.exp(-tmp2)
	opr = 3*s1/s2
	return opr


mytem = np.arange(295)+5

myopr = cal_opr_h2(mytem)

plt.plot(mytem, myopr,'r-')
plt.yscale("log")
plt.show()
