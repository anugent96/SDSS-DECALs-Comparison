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
import itertools
import pandas
df = pandas.read_csv(image1, sep = ',', header=None)
RA = df.values[:,0] #ra
dec = df.values[:,1] #dec
g = df.values[:,2] #g-Magnitude
r = df.values[:,3] #r-Magnitude
z = df.values[:,4] #z-Magnitude
#de Vaucouleurs radius based on filter
rD_g = df.values[:,5]
rD_r = df.values[:,6]
rD_z = df.values[:,7]
#Exponential radius based on filter
rE_g = df.values[:,8]
rE_r = df.values[:,9]
rE_z = df.values[:,10]
#Ellipticity based on filter
deVAB_g = df.values[:,11]
deVAB_r = df.values[:,12]
deVAB_z = df.values[:,13]

expAB_g = df.values[:,14]
expAB_r = df.values[:,15]
expAB_z = df.values[:,16] 


#Strings to floats
ra2 = [float(x) for x in RA if x != 'ra']
dec2 = [float(x) for x in dec if x != 'dec']

g_mag = [float(x) for x in g if x != 'g']
r_mag = [float(x) for x in r if x != 'r']
z_mag = [float(x) for x in z if x != 'z']

rD_g1 = [float(x) for x in rD_g if x != 'deVRad_g']
rD_r1 = [float(x) for x in rD_r if x != 'deVRad_r']
rD_z1 = [float(x) for x in rD_z if x != 'deVRad_z']

rE_g1 = [float(x) for x in rE_g if x != 'expRad_g']
rE_r1 = [float(x) for x in rE_r if x != 'expRad_r']
rE_z1 = [float(x) for x in rE_z if x != 'expRad_z']

dAB_g = [float(x) for x in deVAB_g if x != 'deVAB_g']
dAB_r = [float(x) for x in deVAB_r if x != 'deVAB_r']
dAB_z = [float(x) for x in deVAB_z if x != 'deVAB_z']

eAB_g = [float(x) for x in expAB_g if x != 'expAB_g']
eAB_r = [float(x) for x in expAB_r if x != 'expAB_r']
eAB_z = [float(x) for x in expAB_z if x != 'expAB_z']

# Getting ellipticity value from b/a
from operator import truediv
def ellipticity(x):
    e = truediv((1-x),(1+x))
    return e

dev_e_g = [ellipticity(x) for x in dAB_g if x !=0]
dev_e_r = [ellipticity(x) for x in dAB_r if x !=0]
dev_e_z = [ellipticity(x) for x in dAB_z if x !=0]

exp_e_g = [ellipticity(x) for x in eAB_g if x !=0]
exp_e_r = [ellipticity(x) for x in eAB_r if x !=0]
exp_e_z = [ellipticity(x) for x in eAB_z if x !=0]



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
import numpy as np
from numpy import mean
def magnitude(f): # nanomaggies to magnitudes
   if f <= 0:
      m = 35
   else:
      m = 22.5 - 2.5*math.log10(f)
   return m

aMask = [sum(x) for x in Mask_d]

rE_1 = [word for (word, mask1, mask2) in zip(rE_d, type1_d, aMask) if (mask1 == 'EXP' or mask1 == 'COMP') and  mask2 == 0]
rD_1 = [word for (word, mask1, mask2) in zip(rD_d, type1_d, aMask) if (mask1 == 'DEV'or mask1 == 'COMP') and  mask2 == 0]

g_d = fD_d[:,1]
gMAG = [magnitude(f) for f in g_d]
gMAG_1 = [x for (x, mask) in zip(gMAG, aMask) if mask == 0]

r_d = fD_d[:,2]
rMAG = [magnitude(f) for f in r_d]
rMAG_1 = [x for (x, mask) in zip(rMAG, aMask) if mask == 0]

z_d = fD_d[:,4]
zMAG = [magnitude(f) for f in z_d]
zMAG_1 = [x for (x, mask) in zip(zMAG, aMask) if mask == 0]

eE1_1 = [word for (word, mask1, mask2) in zip(eE1_d, type1_d, aMask) if (mask1 == 'EXP'or mask1 == 'COMP')and  mask2 == 0]
eE2_1 = [word for (word, mask1, mask2) in zip(eE2_d, type1_d, aMask) if (mask1 == 'EXP'or mask1 == 'COMP')and  mask2 == 0]
eD1_1 = [word for (word, mask1, mask2) in zip(eD1_d, type1_d, aMask) if (mask1 == 'DEV'or mask1 == 'COMP')and  mask2 == 0]
eD2_1 = [word for (word, mask1, mask2) in zip(eD2_d, type1_d, aMask) if (mask1 == 'DEV'or mask1 == 'COMP')and  mask2 == 0]

eE = [math.sqrt((x**2)+(y**2)) for (x,y) in zip(eE1_1, eE2_1)]
eD = [math.sqrt((x**2)+(y**2)) for (x,y) in zip(eD1_1, eD2_1)]

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
plt.yscale('log')
plt.show()

# r_mag histogram
plt.hist(rMAG_1, bins=60, color='c', label='DECaLS')
plt.hist(r_mag, bins=60, color='g', label='SDSS', alpha=0.7, histtype = 'stepfilled')
plt.title('SDSS and DECaLS r-Magnitude Histograms')
plt.legend()
plt.yscale('log')
plt.show()

# z_mag histogram
plt.hist(zMAG_1, bins=60, color='c', label='DECaLS')
plt.hist(z_mag, bins=60, color='g', label='SDSS', alpha=0.7, histtype = 'stepfilled')
plt.title('SDSS and DECaLS z-Magnitude Histograms')
plt.legend()
plt.yscale('log')
plt.show()

# Exponential Radius Histograms
plt.hist(rE_r1, bins=np.logspace(0, 2.5, 50), color='g', label='SDSS r', alpha=0.4, histtype = 'stepfilled')
plt.hist(rE_z1, bins=np.logspace(0, 2.5, 50), color='b', label='SDSS z', alpha=0.5, histtype = 'stepfilled')
plt.hist(rE_g1, bins=np.logspace(0, 2.5, 50), color='r', label='SDSS g', alpha=0.4, histtype = 'stepfilled')
plt.hist(rE_1, bins=np.logspace(0, 2.5, 50), color='c', label='DECaLS', alpha=0.5)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.title('SDSS and DECaLS Half Light Radius (EXP) Histograms')
plt.show()

# de Vaucouleurs Radius Histograms
plt.hist(rD_r1, bins=np.logspace(0, 2.5, 50), color='g', label='SDSS r', alpha=0.4, histtype = 'stepfilled')
plt.hist(rD_z1, bins=np.logspace(0, 2.5, 50), color='b', label='SDSS z', alpha=0.5, histtype = 'stepfilled')
plt.hist(rD_g1, bins=np.logspace(0, 2.5, 50), color='r', label='SDSS g', alpha=0.4, histtype = 'stepfilled')
plt.hist(rD_1, bins=np.logspace(0, 2.5, 50), color='c', label='DECaLS', alpha=0.5)
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.title('SDSS and DECaLS Half Light Radius (DEV) Histograms')
plt.show()

#Exponential Ellipticity Hist
plt.hist(exp_e_g, bins=60, color='g', label='SDSS g', alpha=0.4, histtype = 'stepfilled')
plt.hist(exp_e_r, bins=60, color='b', label='SDSS r', alpha=0.5, histtype = 'stepfilled')
plt.hist(exp_e_z, bins=60, color='r', label='SDSS z', alpha=0.4, histtype = 'stepfilled')
plt.hist(eE, bins=60, color='c', label='DECaLS', alpha=0.5)
plt.legend()
plt.yscale('log')
plt.title('SDSS and DECaLS Elliptcity (EXP) Histograms')
plt.show()

#de Vaucouleurs Ellipticity Hist
plt.hist(dev_e_g, bins=60, color='g', label='SDSS g', alpha=0.4, histtype = 'stepfilled')
plt.hist(dev_e_r, bins=60, color='b', label='SDSS r', alpha=0.5, histtype = 'stepfilled')
plt.hist(dev_e_z, bins=60, color='r', label='SDSS z', alpha=0.4, histtype = 'stepfilled')
plt.hist(eD, bins=60, color='c', label='DECaLS', alpha=0.5)
plt.legend()
plt.yscale('log')
plt.title('SDSS and DECaLS Elliptcity (DEV) Histograms')
plt.show()

