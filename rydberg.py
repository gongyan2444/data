import numpy as np


## verify the Rydberg contstant calculation

## nu = Z**2*rm*(1/n**2-1/(n+1)**2)
## rm = r_in/(1+me/mn)  # me is the electron mass, mn the the nucleus mass
## or rm = r_in*(1-me/(me+mn))

## for hydrogen

me = 9.109384e-31

mh = 1.66053906660e-27*1.007825 # from wiki

RH = 3.28984202e15*(mh-me)/mh   # Hz
# 3.2880512911578896 consistent with Wilson+2009

## for Deuterium

md = 1.66053906660e-27*2.01410177811 # from wiki

RD = 3.28984202e15*(md-me)/md   # Hz



try:
    os.system('rm -rf D-alpha.dat')
    infile = open('D-alpha.dat','w')
except:
    infile = open('D-alpha.dat','w')


infile.write("quantum"+' '+"Frequency(MHz)"+'\n')    
    
for i in np.arange(200)+1:
    f = RD*(1./i**2 -1/(i+1)**2)
    infile.write("{:.0f}".format(i)+' '+"{:.3f}".format(f/1e6)+'\n')


infile.close()
print("Done")


