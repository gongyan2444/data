# try hfs 

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import emcee
import corner


def gauss(x,a,b,c):
        return a*np.exp(-(x-b)**2/2/c**2)

def spec(x,a,b,c):
	tau = gauss(x,a,b,c)
	result = 1- np.exp(-tau)
	return result


# 1-0  
# 2/9 -501    
# 4/9 -293    0
# 3/9 0.724
### J= 1-0
v1 = -(-0.501+0.293)/(112359.285-0.293)*2.99792458e5
v2 = 0 
v3 = -(0.724+0.293)/(112359.285-0.293)*2.99792458e5

tau = 1 
dis = 0.4

tau1 = 2/9*tau
tau2 = 4/9*tau 
tau3 = 3/9*tau

xx= np.arange(1001)*0.01-5

yy1 = spec(xx, tau1, v1, dis)
yy2 = spec(xx, tau2, v2, dis)
yy3 = spec(xx, tau3, v3, dis)
yy = yy1 + yy2 + yy3

plt.plot(xx, yy,  label="sum")
plt.plot(xx, yy1, label="F=2/3-5/2")
plt.plot(xx, yy2, label="F=7/2-5/2")
plt.plot(xx, yy3, label="F=5/2-5/2")
plt.legend()
print("ratio")
print(np.max(yy)/np.max(yy2))
#print(np.max(yy2))
#print("ratio: "+np.max(yy)/np.max(yy2))
plt.show()

### J=2-1

#−867 1/25
#−323 64/525
#−213 6/35
#-169 1/3
#-154 1/15
# 358 7/75
# 694 1/63
# 804 2/21
# 902 14/225
vv1 = -(-0.867+0.323) /(224714.368-0.323)*2.99792458e5
vv2 = 0 
vv3 = -(-0.213+0.323) /(224714.368-0.323)*2.99792458e5
vv4 = -(-0.169+0.323) /(224714.368-0.323)*2.99792458e5
vv5 = -(-0.154+0.323) /(224714.368-0.323)*2.99792458e5
vv6 = -(0.358+0.323) /(224714.368-0.323)*2.99792458e5
vv7 = -(0.694+0.323) /(224714.368-0.323)*2.99792458e5
vv8 = -(0.804+0.323) /(224714.368-0.323)*2.99792458e5
vv9 = -(0.902+0.323) /(224714.368-0.323)*2.99792458e5


ttau1 = 1/25.*tau 
ttau2 = 64/525.*tau
ttau3 = 6/35.*tau
ttau4 = 1/3.*tau
ttau5 = 1/15.*tau
ttau6 = 7/75.*tau
ttau7 = 1/63.*tau
ttau8 = 2/21.*tau
ttau9 = 14/225.*tau


yb1 = spec(xx, ttau1, vv1, dis)
yb2 = spec(xx, ttau2, vv2, dis)
yb3 = spec(xx, ttau3, vv3, dis)
yb4 = spec(xx, ttau4, vv4, dis)
yb5 = spec(xx, ttau5, vv5, dis)
yb6 = spec(xx, ttau6, vv6, dis)
yb7 = spec(xx, ttau7, vv7, dis)
yb8 = spec(xx, ttau8, vv8, dis)
yb9 = spec(xx, ttau9, vv9, dis)

yb = yb1+yb2+yb3+yb4+yb5+yb6+yb7+yb8+yb9

plt.plot(xx, yb,  label="sum")
#plt.plot(xx, yb1, label="F=2/3-5/2")
plt.plot(xx, yb2+yb3+yb4+yb5, label="F=5/2-5/2")
#plt.plot(xx, yb3, label="F=7/2-5/2")
#plt.plot(xx, yb4, label="F=9/2-7/2")
#plt.plot(xx, yb5, label="F=1/2-3/2")
#plt.plot(xx, yb6, label="F=3/2-3/2")
#plt.plot(xx, yb7, label="F=5/2-7/2")
#plt.plot(xx, yb8, label="F=7/2-7/2")
#plt.plot(xx, yb9, label="F=5/2-3/2")
plt.legend()
print("ratio")
print(np.max(yb)/np.max(yb2+yb3+yb4+yb5))
#print(np.max(yy2))
#print("ratio: "+np.max(yy)/np.max(yy2))
plt.show()

