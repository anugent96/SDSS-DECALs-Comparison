"""
Finds differences in ellipticity between SDSS and DECaLS and prints the average difference and the average SDSS error.
"""

import sys
image1 = str(sys.argv[1]) # SDSS csv file
image2 = str(sys.argv[2]) # DECaLS fits file

"""
Importing the SDSS csv file into a readable python array
"""
import csv
import itertools
import pandas

df = pandas.read_csv(image1, sep = ',', header=None)
ra_s1 = df.values[:,0] # sdss ra
dec_s1 = df.values[:,1] # sdss dec
deVAB_r = df.values[:,2] # b/a (semi-minor/semi-major axis) ellipticity 
lnDeV_r = df.values[:,3] # ln(likelihood DEV)
deVABErr_r = df.values[:,4] # error in b/a

#Strings to floats
import numpy as np

ra_s = np.array([float(x) for x in ra_s1 if x != 'ra'])
dec_s = np.array([float(x) for x in dec_s1 if x != 'dec'])
dAB_r = np.array([float(x) for x in deVAB_r if x != 'deVAB_r'])
ln_DeV_r = np.array([float(x) for x in lnDeV_r if x != 'lnLDeV_r'])
dABErr_r = np.array([float(x) for x in deVABErr_r if x != 'deVABErr_r'])

from operator import truediv
def ellipticity(x): # b/a to ellipticity constant
    e = truediv((1-x),(1+x))
    return e

from math import exp

DeVLike_r = [exp(x) for x in ln_DeV_r] # likelihod (probability range: 0 > 1)

# If an object has a 90% chance or higher of being DEV and a reasonable uncertainty (< 22)
ra_s2 = [x for (x, y, z) in zip(ra_s, DeVLike_r, dABErr_r) if y > 0.9 and (z < 22 and z > 0)] 
dec_s2 = [x for (x, y, z) in zip(dec_s, DeVLike_r, dABErr_r) if y > 0.9 and (z < 22 and z > 0)]
deV_r = [x for (x, y, z) in zip(dAB_r, DeVLike_r, dABErr_r) if y > 0.9 and (z < 22 and z > 0)]
dabErr_r = [x for (x, y, z) in zip(dABErr_r, DeVLike_r, dABErr_r) if y > 0.9 and (z < 22 and z > 0)]


"""
Importing DECaLS fits file into readable python array
"""
from astropy.io import fits
tractor = fits.open(image2)

tbl = tractor[1].data # data stored within the table
ra_d1 = tbl.field('ra') # DECaLS RA values
dec_d1 = tbl.field('dec') # DECaLS DEC values
eD1_d = tbl.field('shapedev_e1') # DECaLS DEV ellipticity component 1
eD2_d = tbl.field('shapedev_e2') # DECaLS DEV ellipticity component 2
type1 = tbl.field('type')
Mask_d = tbl.field('decam_anymask') # DECaLS masking values

import math
import numpy as np
from numpy import mean

aMask = [sum(x) for x in Mask_d] # sum of masking values in each filter

# Taking out objects with masking values > 0 
# only objects with type DEV or COMP will be considered
ra_d = [x for (x,y,z) in zip(ra_d1, type1, aMask) if (y == 'DEV'or y == 'COMP')and z == 0]
dec_d = [x for (x,y,z) in zip(dec_d1, type1, aMask) if (y == 'DEV'or y == 'COMP')and z == 0]

eD1 = [x for (x, y, z) in zip(eD1_d, type1, aMask) if (y == 'DEV'or y == 'COMP')and z == 0]
eD2 = [x for (x, y, z) in zip(eD2_d, type1, aMask) if (y == 'DEV'or y == 'COMP')and z == 0]

# Finding ellipticity constant from e1 and e2
eD = [math.sqrt((x**2)+(y**2)) for (x,y) in zip(eD1, eD2)]

"""
Matching SDSS and DECaLS and viewing their recorded ellipticity values
"""
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky

# Finding nearest points in DECaLS to SDSS (1st pair)
decals = SkyCoord(ra = ra_d*u.degree, dec = dec_d*u.degree)
sdss = SkyCoord(ra = ra_s2*u.degree, dec = dec_s2*u.degree)
idx, d2d, d3d = match_coordinates_sky(sdss, decals, 1) # Idx = indices in decals to get coordiantes for SDSS

i = 0 
e_diff = [] # difference in ellipticities btwn DECaLS and SDSS
e_s_new = []# new ellipticity list for SDSS (only for matched objects)
e_d_new = []# new ellipticity list for DECaLS (only for matched objects)
sdss_error = [] # new SDSS ellipticity error list

while i < len(idx)-1:
    if d2d[i] <= 0.00027*u.degree: # if objects are within an arcsecond
        e_diff.append(abs(deV_r[i] - eD[idx[i]]))
        e_s_new.append(deV_r[i])
        e_d_new.append(eD[idx[i]])
        sdss_error.append(dabErr_r[i])
    i += 1

# Printing difference and SDSS error
from numpy import mean    
print('avg difference', np.mean(e_diff))
print('max difference', max(e_diff))
print('min difference', min(e_diff))
print('avg SDSS error (b/a)', mean(sdss_error))
