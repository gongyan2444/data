import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def lfun(x, a, b):
	return a*x+b


tq = np.array([9.375, 18.75, 37.5, 75.0, 150.0, 225.0, 300.0, 500.0, 1000., 2000.])
q  = np.array([3.7435, 7.1223, 13.8965, 27.4545, 54.5814, 81.7184, 108.8651, 181.3025, 362.6910, 726.7430])

popt, pcov = curve_fit(lfun, np.log(tq), np.log(q))

tex = np.arange(2000)

#qex=popt[0]*tex+popt[1]

qex = np.exp(np.log(tex)*popt[0]+popt[1])

plt.plot(np.log(tq),np.log(q),'ro', np.log(tex), np.log(qex),'b-')