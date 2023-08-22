from astropy.io import fits
import numpy as np
from regions import PixCoord
from regions.shapes import PolygonPixelRegion
from astropy.nddata import Cutout2D
from regions import Regions
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from matplotlib.path import Path

# read fits files
dat, hd = fits.getdata("N97L_m0.fits", header=True)
# read ds9 polygons
reg1 = Regions.read('test2.reg', format='ds9')
# creating masks
wcs = WCS(hd)
vertices_pixel =  reg1[0].vertices.to_pixel(wcs)
mask = np.zeros_like(dat, dtype=bool)

vertices_pixel = np.array(vertices_pixel).T

path = Path(vertices_pixel)

x, y = np.meshgrid(np.arange(mask.shape[1]), np.arange(mask.shape[0]))
x_flat, y_flat = x.flatten(), y.flatten()

points = np.column_stack((x_flat, y_flat))
mask_indices = path.contains_points(points)
mask_indices = mask_indices.reshape(mask.shape)

mask[mask_indices] = 1

### masked data ###
maskdata = dat*mask
fits.writeto("test.fits", maskdata, hd, overwrite=True)


