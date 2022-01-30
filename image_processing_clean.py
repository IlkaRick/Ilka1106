# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:11:27 2022

@author: anugi
"""


from astropy.table import Table
from astropy.io import fits
import astropy.cosmology as cosmo
from astropy.units import quantity
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

hdulist = fits.open(r'C:\Users\anugi\OneDrive\Documents\Physics\Year 3\Labs\Astronomical Image Processing\A1_mosaic\A1_mosaic.fits')

pixel = hdulist[0].data
print(pixel[0])

columns = len(pixel[0])
rows = len(pixel)
all_pixels = []
for i in range(rows):
    pix_arr = pixel[i]
    for j in range(columns):
        all_pixels.append(pix_arr[j])
print(all_pixels)

plt.hist(all_pixels, bins=900, color='deeppink')
plt.show()

############# LOAD REGIONS TO BE MASKED FOR LEFT SIDE OF IMAGE ################
x, y, width, height, z = sp.loadtxt('C:/Users/anugi/OneDrive/Documents/Physics/Year 3/Labs/Astronomical Image Processing/left_regions.txt', dtype=str, unpack=True)
print(x)
print(width)
print(height)

########### DO LEFT SIDE ##############
mask = np.ones((4611,2570)) # Create mask of 1s to indicate valid objects and background

for k in range(len(x)):
    xcentrefloat = float(x[k])
    ycentrefloat = float(y[k])
    widthfloat = float(width[k])
    heightfloat = float(height[k])
    #widthfloat = width[k]
    #heightfloat = height[k]
    xcentre = int(xcentrefloat) - 1
    ycentre = int(ycentrefloat) - 1
    width1 = (int(widthfloat) + 1)/2
    height1 = (int(heightfloat) + 1)/2    
    widthedge1 = xcentre - width1
    widthedge2 = xcentre + width1
    heightedge1 = ycentre - height1
    heightedge2 = ycentre + height1
    for i in range (4611):
        for j in range (2570):
            if i >= heightedge1 and i <=heightedge2:
                if j >= widthedge1 and j <= widthedge2:
                    mask[i][j] = 0
print(mask)
#yote = mask[0]

####################### LOAD REGIONS FOR RIGHT SIDE ####################
x, y, width, height, z = sp.loadtxt('C:/Users/anugi/OneDrive/Documents/Physics/Year 3/Labs/Astronomical Image Processing/mask_1.txt', dtype=str, unpack=True)
print(x)
print(width)
print(height)

############## DO RIGHT SIDE ###############
for k in range(len(x)):
    xcentrefloat = float(x[k])
    ycentrefloat = float(y[k])
    widthfloat = float(width[k])
    heightfloat = float(height[k])
    #widthfloat = width[k]
    #heightfloat = height[k]
    xcentre = int(xcentrefloat) - 1
    ycentre = int(ycentrefloat) - 1
    width1 = (int(widthfloat) + 1)/2
    height1 = (int(heightfloat) + 1)/2    
    widthedge1 = xcentre - width1
    widthedge2 = xcentre + width1
    heightedge1 = ycentre - height1
    heightedge2 = ycentre + height1
    for i in range (4611):
        for j in range (2570):
            if i >= heightedge1 and i <=heightedge2:
                if j >= widthedge1 and j <= widthedge2:
                    mask[i][j] = 0
print(mask)

############# LOAD EDGES #######################
x, y, width, height, z = sp.loadtxt('C:/Users/anugi/OneDrive/Documents/Physics/Year 3/Labs/Astronomical Image Processing/mask_3.txt', dtype=str, unpack=True)
print(x)
print(width)
print(height)

############## EDGES ###############
for k in range(len(x)):
    xcentrefloat = float(x[k])
    ycentrefloat = float(y[k])
    widthfloat = float(width[k])
    heightfloat = float(height[k])
    #widthfloat = width[k]
    #heightfloat = height[k]
    xcentre = int(xcentrefloat) - 1
    ycentre = int(ycentrefloat) - 1
    width1 = (int(widthfloat) + 1)/2
    height1 = (int(heightfloat) + 1)/2    
    widthedge1 = xcentre - width1
    widthedge2 = xcentre + width1
    heightedge1 = ycentre - height1
    heightedge2 = ycentre + height1
    for i in range (4611):
        for j in range (2570):
            if i >= heightedge1 and i <=heightedge2:
                if j >= widthedge1 and j <= widthedge2:
                    mask[i][j] = 0
print(mask)

###### Choose to cut off background from 3500 #######
################## CUT OFF BACKGROUND #######################

for i in range(4611):
    for j in range(2570): # For loop through whole image
        # girl what the fuck
        impix = pixel[i][j] # Assign variable to pixel
        if impix <= 3500: # Exclude values below chosen background
            mask[i][j] = 0 # Change 1 to 0
print(mask)
######################### FINDING THE SOURCES #####################
pix_value = []
pix_pos = []
copy_pixvalue = []
copy_pixpos = []
for i in range(4611):
    for j in range(2570):
        impix = pixel[i][j]
        maskpix = mask[i][j]
        if maskpix == 1:
            pix_value.append(impix)
            copy_pixvalue.append(impix)
            pix_pos.append(([i,j]))
            copy_pixpos.append(([i,j]))

#%%

'''
NOTE: Indices are reversed:
    Takes form of [y, x]
'''
#pepa = [2,4,1,5,7,3,4,6]
#felix = [[2,0], [4,0], [1,0], [5,0], [7,0], [3,0], [4,0], [6,0]]
pix_value, pix_pos = (list(t) for t in zip(*sorted(zip(pix_value, pix_pos), reverse=True)))
#pix_value = copy_pixvalue
#pix_pos = copy_pixpos
max_value = max(pix_value)
max_index = pix_value.index(max_value) # Find index 
sepindex = pix_pos[max_index]
print(max_value)
print(max_index)
print(pix_pos[max_index])
print(sepindex[0])

############# DEFINE CIRCLE AND ELLIPSE ############
def Circle(x, y, x0, y0):
    r = ((x - x0)**2 + (y - y0)**2)**(1/2)
    return r

def Ellipse(x, y, x0, y0, a, b):
    one = ((x - x0)**2)/(a**2) + ((y - y0)**2)/(b**2)
    return one

   

################## FIND SMALLEST SOURCES ###############
sources = []
for k in range(len(pix_value)):
    print('yo')
    k_index = pix_pos[k]
    kx = k_index[1]
    ky = k_index[0]
    if mask[ky][kx] == 0: # If mask is 0 at pix_value position, skip pix_value
        print('yos')
        continue
    print('Ya maximum is:', max(pix_value))
    #if pix_value[k] == max(pix_value):
    max_value = pix_value[k]
    max_index = pix_pos[k] # make sure index of value in loop matches brightest pixel
    x0 = max_index[1] # Separate x and y coordinates
    y0 = max_index[0]
    r1 = 10.5 # Choose aperture radius from region file containing radius of brightest source
    r2 = 16 # Choose bg aperture size
    bg_aperture = []
    aperture = []
    bg_indices = []
    indices = []
    print('pixvalue:', len(pix_value))
    for i in range(4611):
        for j in range(2570):
            impix = pixel[i][j]
            maskpix = mask[i][j]
            index = [i, j]
            radius = Circle(j, i, x0, y0)
            if radius > r1 and radius <= r2: # Check if pixel is within background aperture
                bg_aperture.append(impix)
                bg_indices.append(index)
            elif radius <= r1:
                aperture.append(impix)
                indices.append(index)
                mask[i][j] == 0
                    #if index in pix_pos:
                        #pixvalue_index = pix_pos.index(index)
                        #pix_value.pop(pixvalue_index)
                        #pix_pos.pop(pixvalue_index)
                        #print('new pixvalue:', len(pix_value), 'and pi:', pi, 'and pv:', pv)
    print('yus')
    fileindex = str(k)
    f = open('Sources/source'+fileindex, 'x') # Create file to put data in
    np.savetxt(f, aperture)
    f.close()
    file = open('Index/index'+fileindex, 'x')
    np.savetxt(file, indices)
    file.close()
    bg = open('bg_aperture/bg_aperture'+fileindex, 'x') # Create file to put data in
    np.savetxt(bg, bg_aperture)
    bg.close()
    gb= open('bg_index/bg_index'+fileindex, 'x')
    np.savetxt(gb, bg_indices)
    gb.close()
    print('yeet')
print('And the maximum pixel count is...', max(pix_value), '!')
