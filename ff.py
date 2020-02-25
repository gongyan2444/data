import numpy as np

def fftime(rho):
	"""
	calculate free-fall timescale of molecular clouds
	input:
	number denseity in units of cm^-3
	return :
	free-fall timescale in units of Myr
	"""
	G = 6.67259e-8
	t = np.sqrt(3.0*np.pi/(32.*G*rho*2.8*1.6733e-24))
	t = t/3600./24./365./1e6
	return t

def ffv(m, r):
	"""
	calculate free-fall velocity of molecular clouds
	input:
	mass in units of M$_\odot$
	r    in units of pc
	return :
	free-fall velocity in units of km/s
	"""
	G = 6.67259e-8
	v = np.sqrt(2*G*m*1.99e33/(r*3.086e18))
	v = v/1e5
	return v

