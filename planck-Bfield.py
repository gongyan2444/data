#
# This file is part of Magnetar
#
# Copyright (C) 2013-2018 Juan Diego Soler
#
# This file is used to show how to make a LIC map with Magnetar.
# 
# This can be used to show large-scale gradients as well.

import sys
import numpy as np
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from tqdm import tqdm
from astropy.wcs import WCS
from scipy import interpolate
import os


def planckct():

   colombi1_cmap = matplotlib.colors.ListedColormap(np.loadtxt('/homes/ygong/soft/magnetar/Planck_Parchment_RGB.txt')/255.)
   colombi1_cmap.set_bad('white') # color of missing pixels
   colombi1_cmap.set_under("blue")

   return colombi1_cmap


# ================================================================================================================================
def lic(vx, vy, length=8, niter=1, normalize=True, amplitude=False, level=0.1, scalar=1, interpolation='nearest', inputmap=None):
   # Calculates the line integral convolution representation of the 2D vector field represented by Vx and Vy.
   # INPUTS
   # Vx     - X
   # Vy     - Y
   # length - L

   vxbad=np.isnan(vx).nonzero()
   vybad=np.isnan(vy).nonzero()

   vx[vxbad]=0. 
   vy[vybad]=0.

   sz=np.shape(vx)

   ni=sz[0]
   nj=sz[1]

   uu=np.sqrt(vx**2+vy**2)
   ii=(uu == 0.).nonzero()

   if (np.size(ii) > 0):
      uu[ii]=1.0
   
   if (normalize):
      ux=vx/uu
      uy=vy/uu
   else: 
      ux=vx/np.max(uu)
      uy=vy/np.max(uu)

   if (inputmap is None):
      vl=np.random.rand(ni,nj)
   else:
      vl=inputmap

   xi=np.arange(ni)
   xj=np.arange(nj)

   outvl=np.zeros([niter,ni,nj])
   
   for i in range(0,niter):

      print('iter {:.0f} / {:.0f}'.format(i+1, niter))

      texture=vl
      vv=np.zeros([ni,nj])

      pi0, pj0 = np.meshgrid(xi, xj, indexing ='ij')
      pi, pj   = np.meshgrid(xi, xj, indexing ='ij')
      mi=pi 
      mj=pj
        
      ppi=1.*pi
      ppj=1.*pj
      mmi=1.*mi
      mmj=1.*mj

      pbar = tqdm(total=length)

      for l in range(0,length):

         ppi0=ppi
         ppj0=ppj
         points   =np.transpose(np.array([pi0.ravel(),pj0.ravel()]))
         outpoints=np.transpose(np.array([ppi.ravel(),ppj.ravel()]))
         dpi=interpolate.griddata(points, uy.ravel(), outpoints, method=interpolation)
         dpj=interpolate.griddata(points, ux.ravel(), outpoints, method=interpolation)

         ppi=ppi0+0.25*np.reshape(dpi,[ni,nj])
         ppj=ppj0+0.25*np.reshape(dpj,[ni,nj])

         mmi0=mmi
         mmj0=mmj
         points   =np.transpose(np.array([pi0.ravel(),pj0.ravel()]))
         outpoints=np.transpose(np.array([mmi.ravel(),mmj.ravel()]))
         dmi=interpolate.griddata(points, uy.ravel(), outpoints, method=interpolation)
         dmj=interpolate.griddata(points, ux.ravel(), outpoints, method=interpolation)

         mmi=mmi0-0.25*np.reshape(dmi,[ni,nj])
         mmj=mmj0-0.25*np.reshape(dmj,[ni,nj])

         pi=(np.fix(ppi) + ni) % ni
         pj=(np.fix(ppj) + nj) % nj
         mi=(np.fix(mmi) + ni) % ni
         mj=(np.fix(mmj) + nj) % nj

         ppi=pi+(ppi.copy()-np.fix(ppi.copy()))
         ppj=pj+(ppj.copy()-np.fix(ppj.copy()))
         mmi=mi+(mmi.copy()-np.fix(mmi.copy()))
         mmj=mj+(mmj.copy()-np.fix(mmj.copy()))

         points   =np.transpose(np.array([pi0.ravel(),pj0.ravel()]))
         outpoints=np.transpose(np.array([ppi.ravel(),ppj.ravel()]))
         tempA=interpolate.griddata(points, texture.ravel(), outpoints, method=interpolation)
   
         points   =np.transpose(np.array([pi0.ravel(),pj0.ravel()]))
         outpoints=np.transpose(np.array([mmi.ravel(),mmj.ravel()]))
         tempB=interpolate.griddata(points, texture.ravel(), outpoints, method=interpolation)

         vv=vv.copy() + np.reshape(tempA,[ni,nj]) + np.reshape(tempB,[ni,nj])

         pbar.update()

      pbar.close()
     
      vl=0.25*vv/length

      outvl[i,:,:]=vl

   vl[vxbad]=np.nan
   vl[vybad]=np.nan

   return outvl


# ================================================================================================================================
def vectors(image, vx, vy, pitch=10, normalize=True, cmap='binary', savefile=False):
   # Calculates the line integral convolution representation of the 2D vector field represented by Vx and Vy.
   # INPUTS
   # Vx     - X
   # Vy     - Y
   # pitch  - 

   sz=np.shape(image)

   nx=sz[0]
   ny=sz[1]

   uu=np.sqrt(vx**2+vy**2)
   ii=(uu == 0.).nonzero()

   if (np.size(ii) > 0):
      uu[ii]=1.0

   if (normalize):
      ux=vx/uu
      uy=vy/uu
   else:
      ux=vx/np.max(uu)
      uy=vy/np.max(uu)
 
   X, Y = np.meshgrid(np.arange(0, sz[1]-1, pitch), np.arange(0, sz[0]-1, pitch))
   ux0=ux[Y,X]
   uy0=uy[Y,X]
   
   return X, Y, ux0, uy0

# ================================================================================================================================

prefix=sys.argv[1]
hdu=fits.open(prefix+"-I.fits")
Imap=hdu[0].data
refhdr1=hdu[0].header 
hdu=fits.open(prefix+'-Q.fits')
Qmap=hdu[0].data
hdu=fits.open(prefix+'-U.fits')
Umap=hdu[0].data

psi=np.arctan2(-Umap,Qmap)

ex=np.cos(psi)
ey=np.sin(psi)
bx=ey
by=-ex

sz=np.shape(Qmap)
length=int(0.1*sz[0]) 

licmap=lic(bx, by, length=length, niter=1)
x, y, ux, uy = vectors(Imap, bx, by, pitch=15)

ax1=plt.subplot(1,1,1, projection=WCS(refhdr1))
im=ax1.imshow(Imap, origin='lower', cmap=planckct(), vmin=0.01, vmax=0.04)
ax1.imshow(licmap.squeeze(), origin='lower', alpha=0.4, cmap='binary', clim=[np.mean(licmap)-np.std(licmap),np.mean(licmap)+np.std(licmap)])
arrows=plt.quiver(x, y, ux, uy, units='width', color='black', pivot='middle', headlength=0, headwidth=0)  
plt.colorbar(im)
plt.show()
