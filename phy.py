import numpy as np


"""
the calculation is mainly carried out in the cgs units
"""

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

def Fctime(rho, A):
	"""
	calculate collapse timescale of filaments (Clarke & Whitworth 2015)
	input:
	rho: number denseity in units of cm^-3
        A: aspect ratio defined as half length divided by width 
	return :
	collapse timescale in units of Myr
	"""
	G = 6.67259e-8
	t = (0.49+0.26*A)/np.sqrt(G*rho*2.8*1.6733e-24)
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


def tau(tex, tmb, mu):
	"""
	calculate opacity of spectra with the radiative transfer equation
	input:
	tex in units of K   (excitation temperature)
	tmb in units of K   (main beam temperature)
	mu  in units of GHz (rest frequency)
	"""
	h = 6.6260755e-27
	k = 1.380658e-16
	mu = mu *1e9
	tbg = 2.73          # cosmic microwave background 
	mj1 = h*mu/k*(np.exp(h*mu/k/tex)-1)**(-1)
	mj2 = h*mu/k*(np.exp(h*mu/k/tbg)-1)**(-1)
	opacity = - np.log(1-tmb/(mj1-mj2))
	return opacity

def alfven(bt, rho):
	"""
	calculate Alfvenic velocity.
	Inputs:
		bt:  magnetic field in units of micro Gauss
		rho: number density in units of cm-2
	return:
		Alfvenic velocity in units of km/s.
	"""
	avelo = bt*1e-6/np.sqrt(4*np.pi*rho*2.8*1.6733e-24)/1e5

	return avelo



def Planckfunc_cgs(freq, temperature):
    """
    Calculate Planck function. Taken from Astrobao
  
    Inputs:
       freq: frequency, in Hz
       temperature: temperature in Kelvin
     
    Return:
       Intensity: in cgs unit ( erg s^-1 sr^-1 cm^-2 Hz-1 )
    """

    # defining physical constants
    c_cgs       = 29979245800.0   # light speed
    h_cgs       = 6.62606885e-27  # planck constant
    kB_cgs      = 1.38064852e-16  # Boltzmann constant

    A = ( 2.0 * h_cgs * (freq**3.0) ) /  ( c_cgs ** 2.0 )
    B = np.exp( (h_cgs * freq) / (kB_cgs * temperature) )
    return A * ( 1.0 / (B - 1.0) )


## this set up gives you the option to run (or not run) a chunk of code when imported from another python file ##

if __name__ == '__main__':
	print("run from python")
else:
	print("run from import")
