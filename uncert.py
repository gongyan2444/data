import numpy as np


tex = 50 
h = 6.6260755e-27
k = 1.380658e-16
mu = 220.39868420 *1e9
bb = 55101.01 * 1e6
tbg = 2.73          # cosmic microwave background 
mj1 = h*mu/k*(np.exp(h*mu/k/tex)-1)**(-1)
mj2 = h*mu/k*(np.exp(h*mu/k/tbg)-1)**(-1)
eu = 15.9
Qrot = k*tex/(h*bb)+1./3



factor = Qrot*(np.exp(eu/tex))/(np.exp(h*mu/k/tex)-1)*1./(mj1-mj2)

print(factor)



