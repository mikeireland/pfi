"""A Script to compute the effective overlap between a laser and a star Airy disk,
when the starlight is dispersed, in a heterodyne laser frequency comb detection
scheme."""

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

#Start with a y(x) jinc function (i.e. an Airy disk) - plot the cross section.
x = 0.01*((np.arange(1000)+0.5) - 500)
y = 2*(sp.jn(1,x*np.pi)/(x*np.pi))
plt.clf()
plt.plot(x,y**2)
plt.plot(x-2,y**2)
plt.plot(x+2,y**2)
plt.plot(x+0.667,y**2,'--')
plt.plot(x-0.667,y**2,'--')
plt.plot([-1,-1,1,1],[-1,1,1,-1], ':', linewidth=2)
plt.axis((-2,2,-0.1,1))
plt.xlabel('Dispersion Axis')
plt.ylabel('Intensity')

#Now for some calculations
hw = 60                                     #Half-width
x = np.arange(2*hw)-hw
xy = np.meshgrid(x,x)                       #A 2D co-ordinate grid.
p = (xy[0]**2 + xy[1]**2) < (hw/12.0)**2    #A circle
im = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(p)))   #A 2D Fourier transform of a circle (i.e. an image)
win = np.zeros((2*hw,2*hw))                             
win[hw-12:hw+12,hw-24:hw+24]=1                          #A square window (i.e. a "pixel")

#Compute the maximum throughput into a pixel
print("Max Star Throughput: {0:6.2f}".format(np.sum(win*np.abs(im)**2)/np.sum(np.abs(im)**2)))

#Compute the minimum throughput into a pixel, when the dispersion places the image 2/3 of the way
#to the edge of the pixel, which we'll take as the limit of our electrical bandwidth.
print("Min Star Throughput: {0:6.2f}".format(np.sum(np.abs(np.roll(im,8,axis=0))**2*win)/np.sum(np.abs(im)**2)))

#Now compute the overlap between laser and starlight at this point - this is the real throughput
#at the edge of the electrical BW.
overlap = np.sum(np.real(np.roll(im,8,axis=0)*np.conj(im)*win))/np.sum(win*np.abs(im)**2)
print("Min Star Overlap: {0:6.2f}".format(overlap))
print("Effective min Throughput: {0:6.2f}".format(overlap**2)) 