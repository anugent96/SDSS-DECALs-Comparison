import sys
image1 = str(sys.argv[1]) # SDSS csv file
image2 = str(sys.argv[2]) # DECaLS tractor.fits file
image3 = str(sys.argv[3]) # SDSS csv file
image4 = str(sys.argv[4]) # DECaLS tractor.fits file
image5 = str(sys.argv[5]) # SDSS csv file
image6 = str(sys.argv[6]) # DECaLS tractor.fits file
"""
Importing the SDSS csv file into a readable python array
"""
import csv
import itertools
import pandas

# 1st SDSS csv file
df = pandas.read_csv(image1, sep = ',', header=None)
ra_s1 = df.values[:,0] # sdss ra
dec_s1 = df.values[:,1] # sdss dec
g_s = df.values[:,2] #sdss g-Magnitude
r_s = df.values[:,3] #sdss r-Magnitude
z_s = df.values[:,4] #sdss z-Magnitude

#Strings to floats
import numpy as np

ra_s = np.array([float(x) for x in ra_s1 if x != 'ra'])
dec_s = np.array([float(x) for x in dec_s1 if x != 'dec'])

g_mag_s = np.array([float(x) for x in g_s if x != 'g'])
r_mag_s = np.array([float(x) for x in r_s if x != 'r'])
z_mag_s = np.array([float(x) for x in z_s if x != 'z'])

# 2nd SDSS csv file
csv = pandas.read_csv(image3, sep = ',', header=None)

ra_s2 = csv.values[:,0] # sdss ra
dec_s2 = csv.values[:,1] # sdss dec
g_s2 = csv.values[:,2] #sdss g-Magnitude
r_s2 = csv.values[:,3] #sdss r-Magnitude
z_s2 = csv.values[:,4] #sdss z-Magnitude

#Strings to floats

ra_s3 = np.array([float(x) for x in ra_s2 if x != 'ra'])
dec_s3 = np.array([float(x) for x in dec_s2 if x != 'dec'])

g_mag_s3 = np.array([float(x) for x in g_s if x != 'g'])
r_mag_s3 = np.array([float(x) for x in r_s if x != 'r'])
z_mag_s3 = np.array([float(x) for x in z_s if x != 'z'])


# 3rd SDSS csv file
sql = pandas.read_csv(image5, sep = ',', header=None)

ra_s4 = sql.values[:,0] # sdss ra
dec_s4 = sql.values[:,1] # sdss dec
g_s4 = sql.values[:,2] #sdss g-Magnitude
r_s4 = sql.values[:,3] #sdss r-Magnitude
z_s4 = sql.values[:,4] #sdss z-Magnitude

#Strings to floats

ra_s5 = np.array([float(x) for x in ra_s4 if x != 'ra'])
dec_s5 = np.array([float(x) for x in dec_s4 if x != 'dec'])

g_mag_s5 = np.array([float(x) for x in g_s4 if x != 'g'])
r_mag_s5 = np.array([float(x) for x in r_s4 if x != 'r'])
z_mag_s5 = np.array([float(x) for x in z_s4 if x != 'z'])

"""
Importing DECaLS fits file into readable python array
"""
from astropy.io import fits
tractor = fits.open(image2)

# 1st DECaLS fits file
tbl = tractor[1].data # data stored within the table
ra_d1 = tbl.field('ra') # DECaLS RA values
dec_d1 = tbl.field('dec') # DECaLS DEC values
fD_d = tbl.field('decam_flux') # DECaLS model flux, decam
Mask_d = tbl.field('decam_anymask') # DECaLS masking values
objid = tbl.field('objid')

import math
import numpy as np
from numpy import mean
def magnitude(f): # nanomaggies to magnitudes
   if f <= 0:
      m = 35
   else:
      m = 22.5 - 2.5*math.log10(f)
   return m

aMask = [sum(x) for x in Mask_d] # sum of masking values in each filter

# Taking out objects with masking values > 0
ra_d = np.array([x for (x,y) in zip(ra_d1, aMask) if y == 0])
dec_d = np.array([x for (x,y) in zip(dec_d1, aMask) if y == 0])

g_d = fD_d[:,1]
g_mag_d = np.array([magnitude(x) for (x, y) in zip(g_d, aMask) if y == 0]) #DECaLS g-magnitude

r_d = fD_d[:,2]
r_mag_d = np.array([magnitude(x) for (x, y) in zip(r_d, aMask) if y == 0]) #DECaLS r-magnitude

z_d = fD_d[:,4]
z_mag_d = np.array([magnitude(x) for (x, y) in zip(z_d, aMask) if y == 0]) #DECaLS z-magnitude

# 2nd DECaLS fits file
hdulist = fits.open(image4)
data = hdulist[1].data # data stored within the table
ra_d2 = data.field('ra') # DECaLS RA values
dec_d2 = data.field('dec') # DECaLS DEC values
fD_d2 = data.field('decam_flux') # DECaLS model flux, decam
Mask_d2 = data.field('decam_anymask') # DECaLS masking values
objid2 = data.field('objid')

aMask2 = [sum(x) for x in Mask_d2] # sum of masking values in each filter

# Taking out objects with masking values > 0
ra_d3 = np.array([x for (x,y) in zip(ra_d2, aMask2) if y == 0])
dec_d3 = np.array([x for (x,y) in zip(dec_d2, aMask2) if y == 0])

g_d2 = fD_d2[:,1]
g_mag_d2 = np.array([magnitude(x) for (x, y) in zip(g_d2, aMask2) if y == 0]) #DECaLS g-magnitude

r_d2 = fD_d2[:,2]
r_mag_d2 = np.array([magnitude(x) for (x, y) in zip(r_d2, aMask2) if y == 0]) #DECaLS r-magnitude

z_d2 = fD_d2[:,4]
z_mag_d2 = np.array([magnitude(x) for (x, y) in zip(z_d2, aMask2) if y == 0]) #DECaLS z-magnitude

# 3rd DECaLS fits file
decam = fits.open(image6)
lists = decam[1].data # data stored within the table
ra_d4 = lists.field('ra') # DECaLS RA values
dec_d4 = lists.field('dec') # DECaLS DEC values
fD_d4 = lists.field('decam_flux') # DECaLS model flux, decam
Mask_d4 = lists.field('decam_anymask') # DECaLS masking values
objid4 = lists.field('objid')

aMask4 = [sum(x) for x in Mask_d4] # sum of masking values in each filter

# Taking out objects with masking values > 0
ra_d5 = np.array([x for (x,y) in zip(ra_d4, aMask4) if y == 0])
dec_d5 = np.array([x for (x,y) in zip(dec_d4, aMask4) if y == 0])

g_d4 = fD_d4[:,1]
g_mag_d4 = np.array([magnitude(x) for (x, y) in zip(g_d4, aMask4) if y == 0]) #DECaLS g-magnitude

r_d4 = fD_d4[:,2]
r_mag_d4 = np.array([magnitude(x) for (x, y) in zip(r_d4, aMask4) if y == 0]) #DECaLS r-magnitude

z_d4 = fD_d4[:,4]
z_mag_d4 = np.array([magnitude(x) for (x, y) in zip(z_d4, aMask4) if y == 0]) #DECaLS z-magnitude
""" 
Comparing magnitudes by like RA and DEC
"""
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.coordinates import match_coordinates_sky

# Finding nearest points in DECaLS to SDSS (1st pair)
decals = SkyCoord(ra = ra_d*u.degree, dec = dec_d*u.degree)
sdss = SkyCoord(ra = ra_s*u.degree, dec = dec_s*u.degree)
idx, d2d, d3d = match_coordinates_sky(sdss, decals, 1) # Idx = indices in decals to get coordiantes for SDSS

i = 0 
g_diff = [] # difference in g-magnitude values between SDSS and DECaLS
r_diff = [] # difference in r-magnitude values between SDSS and DECaLS
z_diff = [] # difference in z-magnitude values between SDSS and DECaLS
# New magnitude lists for SDSS- to ensure that we are only using the points that were in the match
g_s_new = []
r_s_new = []
z_s_new = []
# New magnitude lists for DECaLS
g_d_new = []
r_d_new = []
z_d_new = []

while i < len(idx)-1:
    if d2d[i] <= 0.00027*u.degree:
        g_diff.append(g_mag_s[i] - g_mag_d[idx[i]])
        r_diff.append(r_mag_s[i] - r_mag_d[idx[i]])
        z_diff.append(z_mag_s[i] - z_mag_d[idx[i]])
        g_s_new.append(g_mag_s[i])
        r_s_new.append(r_mag_s[i])
        z_s_new.append(z_mag_s[i])
        g_d_new.append(g_mag_d[idx[i]])
        r_d_new.append(r_mag_d[idx[i]])
        z_d_new.append(z_mag_d[idx[i]])
    i += 1
    
# Finding nearest points in DECaLS to SDSS (2nd pair)
# Adds to lists in first pair
decals2 = SkyCoord(ra = ra_d3*u.degree, dec = dec_d3*u.degree)
sdss2 = SkyCoord(ra = ra_s3*u.degree, dec = dec_s3*u.degree)
idr, r2r, r3r = match_coordinates_sky(sdss2, decals2, 1) # Idr = indices in decals to get coordiantes for SDSS
j = 0
while j < len(idr)-1:
    if r2r[j] <= 0.00027*u.degree:
        g_diff.append(g_mag_s3[j] - g_mag_d2[idr[j]])
        r_diff.append(r_mag_s3[j] - r_mag_d2[idr[j]])
        z_diff.append(z_mag_s3[j] - z_mag_d2[idr[j]])
        g_s_new.append(g_mag_s3[j])
        r_s_new.append(r_mag_s3[j])
        z_s_new.append(z_mag_s3[j])
        g_d_new.append(g_mag_d2[idr[j]])
        r_d_new.append(r_mag_d2[idr[j]])
        z_d_new.append(z_mag_d2[idr[j]])
    j += 1

# Finding nearest points in DECaLS to SDSS (2nd pair)
# Adds to lists in first and second pair
decals3 = SkyCoord(ra = ra_d5*u.degree, dec = dec_d5*u.degree)
sdss3 = SkyCoord(ra = ra_s5*u.degree, dec = dec_s5*u.degree)
idy, y2y, y3y = match_coordinates_sky(sdss3, decals3, 1) # Idy = indices in decals to get coordiantes for SDSS

k = 0
while k < len(idy)-1:
    if y2y[k] <= 0.00027*u.degree:
        g_diff.append(g_mag_s5[k] - g_mag_d4[idy[k]])
        r_diff.append(r_mag_s5[k] - r_mag_d4[idy[k]])
        z_diff.append(z_mag_s5[k] - z_mag_d4[idy[k]])
        g_s_new.append(g_mag_s5[k])
        r_s_new.append(r_mag_s5[k])
        z_s_new.append(z_mag_s5[k])
        g_d_new.append(g_mag_d4[idy[k]])
        r_d_new.append(r_mag_d4[idy[k]])
        z_d_new.append(z_mag_d4[idy[k]])
    k += 1
"""
2D-Histograms of SDSS magnitudes vs. the difference between SDSS and DECaLS magnitudes
"""
import matplotlib.pyplot as plt
plt.hist2d(g_diff, g_s_new, bins=60, range=np.array([(-2, 2), (15, 27)]))
plt.xlabel('Difference between SDSS and DECaLS Mag')
plt.ylabel('SDSS Mag')
plt.title('g-filter')
plt.colorbar()
plt.show()

import matplotlib.pyplot as plt
plt.hist2d(r_diff, r_s_new, bins=60, range=np.array([(-2, 2), (15, 27)]))
plt.xlabel('Difference between SDSS and DECaLS Mag')
plt.ylabel('SDSS Mag')
plt.title('r-filter')
plt.colorbar()
plt.show()

import matplotlib.pyplot as plt
plt.hist2d(z_diff, z_s_new, bins=60, range=np.array([(-2, 2), (15, 27)]))
plt.xlabel('Difference between SDSS and DECaLS Mag')
plt.ylabel('SDSS Mag')
plt.title('z-filter')
plt.colorbar()
plt.show()
