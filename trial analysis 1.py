# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 01:03:04 2022

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


hdulist = fits.open(r"C:\Users\Yi Gan\Documents\Y3 Lab\A1 - astronomical image processing\A1_mosaic.fits")
image = hdulist[0].data
#%%
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
print('slicing succeeded')
#%%
pix_value = []
pix_pos = []
copy_pixvalue = []
copy_pixpos = []
for i in range(4611):
    for j in range(2570):
        impix = image[i][j]
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
    if mask[ky][kx] == 0:
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
            impix = image[i][j]
            maskpix = mask[i][j]
            index = [i, j]
            radius = Circle(j, i, x0, y0)
            if radius > r1 and radius <= r2: # Check if pixel is within background aperture
                bg_aperture.append(impix)
                bg_indices.append(index)
            elif radius <= r1:
                aperture.append(impix)
                indices.append(index)
                mask[i][j] = 0
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
#%%
import os.path
emptylist = [] 
finallist = []
for j in range(112956):
    num = str(j)
    lel = 'Sources\source'+ num
    emptylist.append(lel)
    
for i in range(len(emptylist)):
    emp = emptylist[i]
    file_exist = os.path.exists(emp)
    print(file_exist)
    if file_exist == True:
        finallist.append(i)
        


#%%
sum_sources = [] #sum of pixel values of sources
num_sources = [] #number of pixels present in source file
max_sources = [] #value of maximum pixel in the source file
x_coor = [] #x coordinate of the maximum pixel in the source file
y_coor = [] #y coordinate of the maximum pixel in the source file

max_bg = [] #maximum pixel value in the background pixel file
median_bg = [] #median pixel value in the background pixel file
num_bg = [] #number of pixels present in background file
bgx_coor = [] #x coordinate of the maximum pixel in the bg file
bgy_coor = [] #y coordinate of the maximum pixel in the bg file


for file in range(len(finallist)):
    num = str(finallist[file])
    source_file = np.loadtxt(r'C:\Users\Yi Gan\Documents\Y3 Lab\A1 - astronomical image processing\Sources\source'+num)
    index_y, index_x = np.loadtxt(r'C:\Users\Yi Gan\Documents\Y3 Lab\A1 - astronomical image processing\Index\index'+num, delimiter = ' ', unpack = True)
    source_file = source_file.tolist() #change array to a list to get the index of maximum value
    
    print('sum of the sources pixels =', sum(source_file)) # check for value appended is correct
    sum_sources.append(sum(source_file)) #append to a list
    
    print('number of pixels present = ', len(source_file)) #just to double check the area of aperture
    num_sources.append(len(source_file))
    
    print('maximum pixel value = ', max(source_file)) 
    maximum = max(source_file)
    max_sources.append(maximum)
    index_max = source_file.index(maximum)
    #print(index_max)
    print('x_coor=', index_x[index_max]) #check for correct coordinate appended
    print('y_coor=', index_y[index_max])
    x_coor.append(index_x[index_max])
    y_coor.append(index_y[index_max])
    
    plt.hist(source_file, bins=400)
    plt.title('Plot count histogram source'+num)
    plt.savefig('histogram\plot'+num)
    plt.show()
    
#for file in range(len(first)):
    #num = str(first[file])
    bg_file = np.loadtxt(r'C:\Users\Yi Gan\Documents\Y3 Lab\A1 - astronomical image processing\bg_aperture\bg_aperture'+num)
    bgindex_y, bgindex_x = np.loadtxt(r'C:\Users\Yi Gan\Documents\Y3 Lab\A1 - astronomical image processing\bg_index\bg_index'+num, delimiter = ' ', unpack = True)
    bg_file = bg_file.tolist()
    
    print('number of pixels present = ', len(bg_file)) #just to double check the area of aperture
    num_bg.append(len(bg_file))
    
    print('maximum background=' , max(bg_file)) #check for the value appended
    bg_maximum = max(bg_file)
    max_bg.append(bg_maximum)
    bgindex_max = bg_file.index(bg_maximum)
    
    print('x_coor=', bgindex_x[bgindex_max]) #check for correct coordinate appended
    print('y_coor=', bgindex_y[bgindex_max])
    bgx_coor.append(bgindex_x[bgindex_max])
    bgy_coor.append(bgindex_y[bgindex_max])
    
    median = np.median(bg_file) #local background median
    print('background median',[num],'=', median)
    median_bg.append(median) 
   
    plt.hist(bg_file, bins=400)
    plt.title('Plot count histogram background'+num)
    plt.savefig('background histogram\plot'+num)
    plt.show()

final_sum = []
for i in range (len(sum_sources)):
    background = median_bg[i] # local background pixel
    num_pix = num_sources[i] #number of pixel present in the source
    new_sum = sum_sources[i] - num_pix*background
    final_sum.append(new_sum)

finalsum_arr = np.asarray(final_sum)

cal_mag = [] # calibrated magnitudes
mag_ZP = 2.530E+01
error_ZP = 2.000E-02

for r in range(len(finalsum_arr)):
    mag = mag_ZP - 2.5*np.log10(finalsum_arr[r])
    print('calibrated magnitude = ', mag)
    cal_mag.append(mag)
    
for r in range(len(cal_mag)):
    string = str(cal_mag[r])
    if string == 'nan':
        print('index=', r)
        
data = np.array([finallist, x_coor, y_coor, cal_mag, final_sum, sum_sources, num_sources, max_sources, median_bg, num_bg, max_bg, bgx_coor, bgy_coor ])

data = data.T

np.savetxt('catalog_2.txt', data, delimiter = ',')

#%%
finalsum_arr = np.asarray(final_sum)
cal_mag = []
mag_ZP = 2.530E+01
error_ZP = 2.000E-02

for r in range(len(finalsum_arr)):
    mag = mag_ZP - 2.5*np.log10(finalsum_arr[r])
    print('calibrated magnitude = ', mag)
    cal_mag.append(mag)
  #%%  
for r in range(len(cal_mag)):
    string = str(cal_mag[r])
    if string == 'nan':
        print('index=', r)
#%%
finallist, x_coor, y_coor, cal_mag, final_sum, sum_sources, num_sources, max_sources, median_bg, num_bg, max_bg, bgx_coor, bgy_coor = np.loadtxt(r'C:\Users\Yi Gan\Documents\Y3 Lab\A1 - astronomical image processing\catalog_3.txt', delimiter = ',', unpack=True)        

#%%
plt.hist(cal_mag)
#%%
N_6 = []
N_8 = []
N_10 = []
N_12 = []
N_14 = []
N_16 = []
N_18 = []
for i in range(len(finallist)):
    if cal_mag[i] >= 6:
        N_6.append(cal_mag[i])
    if cal_mag[i] >= 8:
        N_8.append(cal_mag[i])
    if cal_mag[i] >= 10:
        N_10.append(cal_mag[i])
    if cal_mag[i] >= 12:
        N_12.append(cal_mag[i])
    if cal_mag[i] >= 14:
        N_14.append(cal_mag[i])
    if cal_mag[i] >= 16:
        N_16.append(cal_mag[i])
    if cal_mag[i] >= 18:
        N_18.append(cal_mag[i])
len_6= len(N_6)
len_8 = len(N_8)
len_10 = len(N_10)
len_12 = len(N_12)
len_14 = len(N_14)
len_16 = len(N_16)
len_18 = len(N_18)        
#%%
x = [6,8,10,12,14,16,18]
y = [len_6,len_8,len_10,len_12,len_14,len_16,len_18]
log_y = []
for i in range(len(y)):
    ylog = np.log10(y[i])
    print(ylog)
    log_y.append(ylog)
plt.plot(x,log_y)
plt.plot(x,log_y,'*')

#%%
print(max(cal_mag), 'is the maximum!')
print(min(cal_mag), 'is the minimum!')
maglist = np.linspace(9, 20, num=11, endpoint=True) # Create magnitude list
print(maglist)

number = [] # number of sources
magmedian = [] # median of bin
numbin = []
error_up = []
error_down = []
for l in range(len(maglist)-1):
    mag_bin = []
    binend1 = maglist[l] # assign variables to either end of the bin
    binend2 = maglist[l+1]
    for q in range(len(cal_mag)):
        mag = cal_mag[q]
        if binend1 <= mag < binend2: # check if magnitude of source falls within bin
            mag_bin.append(mag)
    
    binnum = len(mag_bin)
    sum_mag = np.log10(binnum) # log of number of sources in bin    
    error = np.sqrt(binnum)
    up = abs(np.log10(binnum+error) - sum_mag)
    error_up.append(up)
    down = abs(sum_mag - np.log10(binnum - error))
    error_down.append(down)
    number.append(sum_mag)
    mmedian = (binend1 + binend2)/2 # find middle of bin
    magmedian.append(mmedian)
    numbin.append(binnum)
r = 7
fit, cov = sp.polyfit(magmedian[0:r], number[0:r], 1, w= error_down[0:r], cov=True)
print('fit coeffs:', fit)
print('cov matrix:', cov)

sig_0 = sp.sqrt(cov[0,0]) #The uncertainty in the slope
sig_1 = sp.sqrt(cov[1,1]) #The uncertainty in the intercept

print('Slope = %.3e +/- %.3e' %(fit[0],sig_0))# Note the %.3e forces the values to be printed in scientific notation with 3 decimal places.
print('Intercept = %.3e +/- %.3e' %(fit[1],sig_1))

polycule=sp.poly1d(fit)
print('polynomial throuple:', polycule)

print(sum(numbin))
plt.plot(magmedian, number, marker='x', linestyle='')
plt.errorbar(magmedian, number, yerr = (error_down, error_up), capsize =5, ls='')
plt.plot(magmedian[0:r],polycule(magmedian[0:r]), color='gold', lw=2)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlabel('Absolute magnitude', fontsize = 30)
plt.ylabel('Number of sources $\log10{N(m)}$',fontsize = 30)
plt.show()
