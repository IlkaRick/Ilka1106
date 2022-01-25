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

#masked right regions
x, y, width, height, angle = np.loadtxt(r"C:\Users\Yi Gan\Documents\Y3 Lab\A1 - astronomical image processing\mask_1.txt", delimiter=' ', unpack = True)
mask= np.ones((4611,2570))
for k in range(len(x)):
    xcentrefloat = float(x[k]) # Convert data from file to floats 
    ycentrefloat = float(y[k])
    widthfloat = float(width[k])
    heightfloat = float(height[k])
    xcentre = int(xcentrefloat) - 1 # Convert floats to intgers
    ycentre = int(ycentrefloat) - 1
    width1 = (int(widthfloat) + 1)/2 # Divide wdith and height by 2 to add and subtract from centres
    height1 = (int(heightfloat) + 1)/2    
    widthedge1 = xcentre - width1 # Find edges of regions
    widthedge2 = xcentre + width1
    heightedge1 = ycentre - height1
    heightedge2 = ycentre + height1
    for i in range (4611):
        for j in range (2570):
            if i >= heightedge1 and i <=heightedge2:
                if j >= widthedge1 and j <= widthedge2:
                    mask[i][j] = 0
print(mask)
print('right regions succeed')
#masked left regions
x, y, width, height, angle = sp.loadtxt(r"C:\Users\Yi Gan\Documents\Y3 Lab\A1 - astronomical image processing\left_regions.txt",dtype = str, unpack = True)
for k in range(len(x)):
    xcentrefloat = float(x[k]) # Convert data from file to floats 
    ycentrefloat = float(y[k])
    widthfloat = float(width[k])
    heightfloat = float(height[k])
    xcentre = int(xcentrefloat) - 1 # Convert floats to intgers
    ycentre = int(ycentrefloat)- 1
    width1 = (int(widthfloat) + 1)/2 # Divide wdith and height by 2 to add and subtract from centres
    height1 = (int(heightfloat) + 1)/2    
    widthedge1 = xcentre - width1 # Find edges of regions
    widthedge2 = xcentre + width1
    heightedge1 = ycentre - height1
    heightedge2 = ycentre + height1
    for i in range (4611):
        for j in range (2570):
            if i >= heightedge1 and i <=heightedge2:
                if j >= widthedge1 and j <= widthedge2:
                    mask[i][j] = 0
print(mask)
print('run left region succeed')
#masked bigger area of the edges to get rid of the strange peak at 3421
x, y, width, height, angle = sp.loadtxt(r"C:\Users\Yi Gan\Documents\Y3 Lab\A1 - astronomical image processing\mask_3.txt",dtype = str, unpack = True)
for k in range(len(x)):
    xcentrefloat = float(x[k]) # Convert data from file to floats 
    ycentrefloat = float(y[k])
    widthfloat = float(width[k])
    heightfloat = float(height[k])
    xcentre = int(xcentrefloat) - 1# Convert floats to intgers
    ycentre = int(ycentrefloat) - 1
    width1 = (int(widthfloat) + 1)/2 # Divide wdith and height by 2 to add and subtract from centres
    height1 = (int(heightfloat) + 1)/2    
    widthedge1 = xcentre - width1 # Find edges of regions
    widthedge2 = xcentre + width1
    heightedge1 = ycentre - height1
    heightedge2 = ycentre + height1
    for i in range (4611):
        for j in range (2570):
            if i >= heightedge1 and i <=heightedge2:
                if j >= widthedge1 and j <= widthedge2:
                    mask[i][j] = 0
print(mask)
print('run background/edges succeed')

#setting a cutting edge for the pixels
for i in range(4611):
    for j in range(2570): # For loop through whole image
        impix = image[i][j] # Assign variable to pixel
        if impix <= 3500: # Exclude values below chosen background
            mask[i][j] = 0 # Change 1 to 0
#%%
pot_sources = [] # potential sources
coor = []
x_axis = []
y_axis = []
for i in range(4611):
    for j in range(2570):
        x = j
        y = i
        if mask[i][j] == 1:
            x_axis.append(x)
            y_axis.append(y)
            coor.append([y,x])
            pot_sources.append(image[i][j])
            maximum = max(pot_sources)
print(maximum)   #check if the first 3 loops work  
#%%     
x_range = []
y_range = []
total_pixel = []
lol = 50
lel = 50
for num in range(len(pot_sources)):
    if pot_sources[num] == maximum:
        x_centre = x_axis[num] - 1
        y_centre = y_axis[num] - 1
        y_up = y_centre + 20
        y_low = y_centre - lol
        x_right = x_centre + 25
        x_left = x_centre - lel
        print(x_centre)
        print(y_centre)
        print(x_left)
        print(x_right)
        print(y_up)
        print(y_low)
        for i in range(4611):
            for j in range(2570):
                if x_left <= j <= x_right:
                    if y_low <= i <= y_up:
                        if mask[i][j] == 1:
                            x_range.append(j)
                            y_range.append(i)
                            total_pixel.append(image[i][j])
                            mask[i][j] = 0

#%%
plt.plot(x_range,y_range,'+')
#%%
plt.hist(total_pixel)
