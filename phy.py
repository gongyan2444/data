import numpy as np
import scipy.optimize

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

def cal_opr_h2(tt):
	"""
	calculate H2 ortho to para ratio. Taken from Fujun Du's code.
	The formulate is based on Eq. 1 in https://ui.adsabs.harvard.edu/abs/2001ApJ...561..254T/abstract
	The paper's Fig. 1 is well reproduced. 
	"""
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


def vlsr_to_vhelio(ra_deg,dec_deg,vlsr):
    """
        LSR defined by removing peculiar Solar Motion
        of 20.0 km/s toward 18.0 hours, +30.0 degree(RA, DEC)
        (Defined in 1900.0 system, so precess to J2000)
    """
    Vsun = 20.0 #kms
    RA_Vsun = np.deg2rad(18.0640000*15.0) #18 hours in J1900
    DE_Vsun = np.deg2rad(30.0024000)      #30 deg in J1900

    cos_RA =np.cos(RA_Vsun)
    sin_RA =np.sin(RA_Vsun)

    cos_DE =np.cos(DE_Vsun)
    sin_DE =np.sin(DE_Vsun)

    ra_rad = np.deg2rad(ra_deg)
    dec_rad = np.deg2rad(dec_deg)

    cos_alpha = np.cos(ra_rad)
    sin_alpha = np.sin(ra_rad)

    cos_delta = np.cos(dec_rad)
    sin_delta = np.sin(dec_rad)

    #solar pecular motion in Equatorial Cartesian frame
    Xo = Vsun*cos_RA*cos_DE        # toward RA=0h
    Yo = Vsun*sin_RA*cos_DE        # toward RA=6h in equat plane
    Zo = Vsun*sin_DE               # toward Dec=90 (North)

    v_proj = -Xo*cos_alpha*cos_delta - Yo*sin_alpha*cos_delta - Zo*sin_delta
    vhelio = vlsr+v_proj

    return vhelio




def mcratio(a, b, da, db):
    '''
    estimate the uncertainties of the ratio (a/b) using Monte Carlo
    a        Numerator  
    b        denominator
    da       1 sigma error of Numerator
    db       1 sigma error of denominator
    '''
    ratio = a/b
    s1 = np.random.normal(a,da,10000)
    s2 = np.random.normal(b,db,10000)
    sample = s1/s2
    err = np.std(sample)
    err2 = 0.5*(np.percentile(sample, 84)- np.percentile(sample, 16))
    return ratio, err, err2

def wmean(a,b, da, db):
    '''
    estimate the uncertainties of the weighted mean using Monte Carlo
    a        Numerator  
    b        denominator
    da       1 sigma error of Numerator
    db       1 sigma error of denominator
    '''
    val = (a/da**2+b/db**2)/(1/da**2+1/db**2)
    s1 = np.random.normal(a,da,10000)
    s2 = np.random.normal(b,db,10000)
    sample = (s1/da**2+s2/db**2)/(1/da**2+1/db**2)
    err = np.std(sample)
    err2 = 0.5*(np.percentile(sample, 84)- np.percentile(sample, 16))
    return val, err, err2


def caltau(iso, ratio):
	def F(x):
		return (1-np.exp(-x))/(1-np.exp(-x/iso))-ratio

	tau = scipy.optimize.broyden1(F, [0.1], f_tol=1e-14)
	return tau
## this set up gives you the option to run (or not run) a chunk of code when imported from another python file ##

if __name__ == '__main__':
	print("run from python")
else:
	print("run from import")
