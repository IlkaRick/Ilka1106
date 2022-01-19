# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 09:24:14 2022

@author: Yi Gan
"""

import numpy as np
import scipy as sp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import astropy
from astropy.table import Table
from astropy.io import fits
import astropy.cosmology as cosmo
import matplotlib as mpl


hdulist = fits.open(r"C:\Users\Yi Gan\Documents\Y3 Lab\A1 - astronomical image processing\A1_mosaic\A1_mosaic.fits")
#%%
image = hdulist[0].data
plt.hist(image)

pix = []
for i in range(len(image)):
    for j in range(len(image[0])):
        pix.append(image[i,j])
        
plt.hist(pix,bins=1000)
#%%
zoom = []
for p in range(len(pix)):
    if pix[p] >=3300 and pix[p] <= 3600:
        zoom.append(pix[p])


mean = np.mean(zoom)
sigma = []
for k in range(len(zoom)):
    sigma1 = (zoom[k]-mean)**2/len(zoom)
    sigma.append(sigma1)
guess_sigma = np.mean(sigma)
values, bins, patches = plt.hist(zoom, bins = 300)
test_mean = mean
test_sigma = 50
test_norm = 1e8
def gauss_function(x, sigma, mu, norm):
    gaus = norm*(1/(sigma*(sp.sqrt(2*sp.pi))))*(sp.exp(-(1/2)*((x-mu)/sigma)**2))
    return(gaus)
p0 = [test_sigma,test_mean,test_norm]
fit_gauss, fit_cov = curve_fit(gauss_function, xdata = bins[0:300], ydata = values, p0 = p0, maxfev=5000)
fitted_gauss = gauss_function(bins[0:300],*fit_gauss)

print(fit_gauss)
print(" sigma = " + str(fit_gauss[0])+ '+/-' +str (np.sqrt(fit_cov[0,0])))
print(" mean = " + str(fit_gauss[1])+ '+/-' +str (np.sqrt(fit_cov[1,1])))
print(" normalisation = " + str(fit_gauss[2])+ '+/-' +str (np.sqrt(fit_cov[2,2])))



plt.plot(bins[0:300], fitted_gauss)
plt.show()

zoom = []
for p in range(len(pix)):
    if pix[p] >=4000: #and pix[p] <= 4000:
        zoom.append(pix[p])

plt.hist(zoom, bins = 100)
#%%
mask_ones = np.ones((4611,2570))
for i in range (len(image)):
    for j in range (len(image[0])):
        if image[i,j] <= 3421:
            mask_ones[i,j] = 0