# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 09:57:21 2022

@author: ykg19
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

#%%
hdulist = fits.open(r"\\icnas4.cc.ic.ac.uk\ykg19\Y3 lab\Astronomical Imaging\A1_mosaic.fits")
#%%
image = hdulist[0].data


#masked right regions
x, y, width, height, angle = np.loadtxt(r"\\icnas4.cc.ic.ac.uk\ykg19\Y3 lab\Astronomical Imaging\mask_1.txt", delimiter=' ', unpack = True)
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
x, y, width, height, angle = sp.loadtxt(r"\\icnas4.cc.ic.ac.uk\ykg19\Y3 lab\Astronomical Imaging\left_regions.txt",dtype = str, unpack = True)
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
x, y, width, height, angle = sp.loadtxt(r"\\icnas4.cc.ic.ac.uk\ykg19\Y3 lab\Astronomical Imaging\mask_3.txt",dtype = str, unpack = True)
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
############# DEFINE CIRCLE AND ELLIPSE ############
def Circle_radius(x, y, x0, y0):
    r = ((x - x0)**2 + (y - y0)**2)**(1/2)
    return r

def Ellipse(x, y, x0, y0, a, b):
    one = ((x - x0)**2)/(a**2) + ((y - y0)**2)/(b**2)
    return one
#%%
pot_sources = [] # potential sources
coor = []
x_pos = []
y_pos = []
x_axis = []
y_axis = []
x_range = []
y_range = []
#total_pixel = []
for i in range(len(image)):
    for j in range(len(image[0])):
        x_pos.append(j)
        y_pos.append(i)
        coor.append([j,i])
        if mask[i][j] == 1:
            x_axis.append(j)
            y_axis.append(i)
            #coor.append([j,i])
            pot_sources.append(image[i][j])
            maximum = max(pot_sources)
            index_max = pot_sources.index(maximum)
            x_ctr = x_axis[index_max] - 1 #centre coordinate
            y_ctr = y_axis[index_max] - 1
#%%
        for num in range(len(coor)):
            x_coor = x_pos[num] - 1
            y_coor = y_pos[num] - 1
            x_rad = abs(x_coor - x_centre)
            y_rad = abs(y_coor - y_centre)
            if coor[num] == [y_ctr,x_ctr]:                 
                
                 
#%%        
        total_pixel = []
        for num in range(len(pot_sources)):
            if num == index_max:
                x_centre = x_axis[num] - 1
                y_centre = y_axis[num] - 1
            x_rad = abs(x_axis[num] - 1 - x_centre)
            y_rad = abs(y_axis[num] - 1 - y_centre)
            x_coor = x_axis[num] - 1
            y_coor = y_axis[num] - 1
            radius = Circle_radius(x_coor, y_coor, x_centre, y_centre) 
            #print('x = ', x_centre)
            #print(x_rad)
            
            if x_rad <= 10:
                if y_rad <= 10:
                    if radius  < 10:
                        print(pot_sources[num])
                        total_pixel.append(pot_sources[num])
print('total_pixel')
#%%
        if x == x_centre:
            if y == y_centre:
                x_rad1 = abs(x[ - 1 - x_centre)
                y_rad1 = abs(y_axis[num] - 1 - y_centre)                
                        
                        print('x =', x_coor)
                        print(x_centre)
                        print('y =', y_coor)
                        print(y_centre)
                        print(radius)
#print(maximum)   #check if the first 3 loops work  
#%%     

lol = 20
lel = 20
for num in range(len(pot_sources)):
    if pot_sources[num] == maximum:
        x_centre = x_axis[num] - 1
        y_centre = y_axis[num] - 1
        y_up = y_centre + lol
        y_low = y_centre - lol
        x_right = x_centre + lel
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
                            #mask[i][j] = 0

#%%
plt.plot(x_range,y_range,'+')
#%%
plt.hist(total_pixel)