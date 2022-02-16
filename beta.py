from phy import * 
import numpy as np

c_cgs       = 29979245800.0
mu450  = c_cgs/0.0450
mu850  = c_cgs/0.0850

ratio = 1.7 

up  = np.log(ratio) - np.log(Planckfunc_cgs(mu450, 10)/Planckfunc_cgs(mu850, 10))
low = np.log(mu450/mu850)

beta = up/low

print(beta)