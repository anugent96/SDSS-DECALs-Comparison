"""
Takes a SDSS csv file and DECaLS fits file and plots histograms of RA, DEC, half-light radius and 
magnitudes.
"""

import sys
image1 = str(sys.argv[1]) # SDSS csv file
image2 = str(sys.argv[2]) # DECaLS tractor.fits file

"""
Importing the SDSS csv file into a readable python array
"""
import csv

import pandas
df = pandas.read_csv(image1, sep = ',', header=None)

RA = df.values[:,0]
dec = df.values[:,1]
g = df.values[:,2]
r = df.values[:,3]
z = df.values[:,4]

ra2 = [float(x) for x in RA if x != 'ra']
dec2 = [float(x) for x in dec if x != 'dec']
g_mag = [float(x) for x in g if x != 'g']
r_mag = [float(x) for x in r if x != 'r']
z_mag = [float(x) for x in z if x != 'z']

"""
Importing DECaLS fits file into readable python array
"""
from astropy.io import fits
tractor = fits.open(image2)

tbl = tractor[1].data # data stored within the table
ra_d = tbl.field('ra') # RA values
dec_d = tbl.field('dec') # DEC values
rE_d = tbl.field('shapeExp_r') # half-light radius, exponential model
rD_d = tbl.field('shapeDev_r') # half-light radius, devaucouleurs model
fD_d = tbl.field('decam_flux') # model flux, decam
eE1_d = tbl.field('shapeExp_e1') # ellipticity component 1, exponential model
eD1_d = tbl.field('shapeDev_e1') # ellipticity component 1, devaucouleurs model
eE2_d = tbl.field('shapeExp_e2') # ellipticity component 2, exponential model
eD2_d = tbl.field('shapeDev_e2') # ellipticity component 2, devaucouleurs model
type1_d = tbl.field('type') # type: EXP, DEV, PSF, SIMP
Mask_d = tbl.field('decam_anymask') #check this- any mask (for errors in telescope)
PSF_d = tbl.field('decam_psfsize') # PSF (Point Spread Function) Size
d_nobs = tbl.field('decam_nobs') # Decam number of exposures


import math
import itertools
from numpy import mean
def magnitude(f): # nanomaggies to magnitudes
   if f <= 0:
      m = 35
   else:
      m = 22.5 - 2.5*math.log10(f)
   return m

aMask = [sum(x) for x in Mask_d]

rE_1 = [word for (word, mask1, mask2) in zip(rE_d, type1_d, aMask) if mask1 == 'EXP' and  mask2 == 0]
rD_1 = [word for (word, mask1, mask2) in zip(rD_d, type1_d, aMask) if mask1 == 'DEV' and  mask2 == 0]

g_d = fD_d[:,1]
gMAG = [magnitude(f) for f in g_d]
gMAG_1 = [x for (x, mask) in zip(gMAG, aMask) if mask == 0]

r_d = fD_d[:,2]
rMAG = [magnitude(f) for f in r_d]
rMAG_1 = [x for (x, mask) in zip(rMAG, aMask) if mask == 0]

z_d = fD_d[:,4]
zMAG = [magnitude(f) for f in z_d]
zMAG_1 = [x for (x, mask) in zip(zMAG, aMask) if mask == 0]

"""
Histograms
"""

import matplotlib.pyplot as plt

# RA histogram
plt.hist(ra_d, bins=60, color='c', label='DECaLS')
plt.hist(ra2, bins=60, color='g', label='SDSS', alpha=0.7, histtype = 'stepfilled')
plt.title('SDSS and DECaLS RA Histograms')
plt.legend()
plt.show()

# dec histogram
plt.hist(dec_d, bins=60, color='c', label='DECaLS')
plt.hist(dec2, bins=60, color='g', label='SDSS', alpha=0.7, histtype = 'stepfilled')
plt.title('SDSS and DECaLS DEC Histograms')
plt.legend()
plt.show()

# g_mag histogram
plt.hist(gMAG_1, bins=60, color='c', label='DECaLS')
plt.hist(g_mag, bins=60, color='g', label='SDSS', alpha=0.7, histtype = 'stepfilled')
plt.title('SDSS and DECaLS g-Magnitude Histograms')
plt.legend()
plt.show()

# r_mag histogram
plt.hist(rMAG_1, bins=60, color='c', label='DECaLS')
plt.hist(r_mag, bins=60, color='g', label='SDSS', alpha=0.7, histtype = 'stepfilled')
plt.title('SDSS and DECaLS r-Magnitude Histograms')
plt.legend()
plt.show()

# z_mag histogram
plt.hist(zMAG_1, bins=60, color='c', label='DECaLS')
plt.hist(z_mag, bins=60, color='g', label='SDSS', alpha=0.7, histtype = 'stepfilled')
plt.title('SDSS and DECaLS z-Magnitude Histograms')
plt.legend()
plt.show()