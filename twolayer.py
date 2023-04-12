## This model is taken from Di Francesco+2001
## The model is based on Myers+1996 but including the free-free central source. 
##

import numpy as np
import matplotlib.pyplot as plt
#import emcee
import multiprocessing
#import corner

ncpu = multiprocessing.cpu_count()

def jfunc(t, nu):
    """
    t- kelvin
    nu - Hz
    """
    H=6.6260755e-27
    K=1.380658e-16
    to = H*nu/K
    return to/(np.exp(to/t)-1.0)



def twolayer(v, v_lsr, freq, tf, tr, tau_f0, tau_r0, sigma1, sigma2):
	tb = 2.73
	tau_f = tau_f0*np.exp(-(v-v_lsr)**2/(2*sigma1**2))
	tau_r = tau_r0*np.exp(-(v-v_lsr)**2/(2*sigma2**2))
	brightness = jfunc(tf, freq)*(1-np.exp(-tau_f)) + jfunc(tr, freq)*(1-np.exp(-tau_r))*np.exp(-tau_f) - jfunc(tb,freq)*(1.-np.exp(-tau_r-tau_f))
	return brightness



xx = np.arange(1001)*0.1-50
v_lsr = 0 
freq = 89.18852470*1e9
tf = 2.73
tr = 4.5 
sigma1 = 0.2 
sigma2 = 0.1
tau_f0 = 4
tau_r0 = 5

yy = twolayer(xx, v_lsr, freq, tf, tr, tau_f0, tau_r0, sigma1, sigma2)

plt.step(xx, yy, where='mid')
plt.show()

