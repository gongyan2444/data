"""
This script follows the method of Stephens+2017
This used Monte Carlo simulations to show 
the expected observed distribution of angles 
of two vectors projected into three dimensions
"""

import numpy as np
import matplotlib.pyplot as plt 

u1 = np.random.random(1000000)*2-1
theta1 = np.random.random(1000000)*2*np.pi
x1 = np.sqrt(1-u1**2)*np.cos(theta1)
y1 = np.sqrt(1-u1**2)*np.sin(theta1)
z1 = u1

u2 = np.random.random(1000000)*2-1
theta2 = np.random.random(1000000)*2*np.pi
x2 = np.sqrt(1-u2**2)*np.cos(theta2)
y2 = np.sqrt(1-u2**2)*np.sin(theta2)
z2 = u2

value3d = (x1*x2+y1*y2+z1*z2)/np.sqrt(x1**2+y1**2+z1**2)/np.sqrt(x2**2+y2**2+z2**2)
value2d = (y1*y2+z1*z2)/np.sqrt(y1**2+z1**2)/np.sqrt(y2**2+z2**2)

alpha3d = np.rad2deg(np.arccos(value3d))
alpha3d[np.where(alpha3d >90)] = 180- alpha3d[np.where(alpha3d >90)] 
alpha2d = np.rad2deg(np.arccos(value2d))
alpha2d[np.where(alpha2d > 90)] = 180 - alpha2d[np.where(alpha2d > 90)]


## random distribution 
values, base = np.histogram(alpha2d, bins=40)
cumulative = np.cumsum(values)      
cumulative = np.hstack((0, cumulative))
plt.plot(base, cumulative/np.max(cumulative),c='r')

## 0-20 degree in 3d 
new1 = alpha2d[np.where((alpha3d >= 0) & (alpha3d <=20))] 
values1, base1 = np.histogram(new1, bins=40)
cumu1 = np.cumsum(values1)      
cumu1 = np.hstack((0, cumu1))
plt.plot(base1, cumu1/np.max(cumu1),c='g')


## 70-90 degree in 3d
new2 = alpha2d[np.where((alpha3d >= 70) & (alpha3d <=90))] 
values2, base2 = np.histogram(new2, bins=40)
cumu2 = np.cumsum(values2)      
cumu2 = np.hstack((0, cumu2))
plt.plot(base2, cumu2/np.max(cumu2),c='b')


for i in np.arange(9):
	idx1 = int((i+1)/10.*60000)
	idx2 = int((1-(i+1)/10) *60000)
	bim1 = np.concatenate((new1[0:idx1],new2[0:idx2])) 	
	bval1, bbas1 = np.histogram(bim1, bins=40)
	bcumu1 = np.cumsum(bval1)      
	bcumu1 = np.hstack((0, bcumu1))
	plt.plot(bbas1, bcumu1/np.max(bcumu1),c='tab:orange',ls='--')


plt.show()




