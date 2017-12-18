"""A script to plot a circular array geometry with 20 telescopes, and output the array elements"""

from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.signal as sig

import scipy.ndimage as nd
import pdb

#Parameters here
outfile = "circ_array.txt"
ntel = 12
short_bl = 550
fov_bl = 30 #To define the field of view.
jitter = -1 ## 400 ##-1 for redundant arrays...
extra_jitter = 5
r_circ = ntel*short_bl/2.0/np.pi
dphi = 2.0*np.pi/ntel
wave=11e-6  #Wavelength in m

xc = np.zeros(ntel)
yc = np.zeros(ntel)

nrot = 59
rotdeg = 1
zooms = [1.02,1.04,1.06,1.08,1.10]
sz = 1024
ring_rad = 0.15/140.
ring_flux = 2.5 #Ratio of ring to star

#----- Automatic from here ------
fov = np.degrees(wave/fov_bl)*3600.

#Create the telescope positions xc and yc.
for i in range(ntel):
    #f(r) = 2r
    #int[f(r)] = r^2
    if (jitter > 0):
        rand_dr   = np.sqrt(np.random.random())*jitter
        rand_dphi = np.random.random()*2.0*np.pi
        xc[i] = np.cos(i*dphi)*r_circ + np.cos(rand_dphi)*rand_dr
        yc[i] = np.sin(i*dphi)*r_circ + np.sin(rand_dphi)*rand_dr
    elif (jitter < 0):
        if (i % 4 == 1): 
            dphi_offset = -1 + 4*4/13
        elif (i % 4 == 2):
            dphi_offset = -2 + 5*4/13
        elif (i % 4 == 3):
            dphi_offset = -3 + 7*4/13
        else:
            dphi_offset = 0.0
        rand_dphi = dphi*dphi_offset
        xc[i] = np.cos(i*dphi + rand_dphi)*r_circ 
        yc[i] = np.sin(i*dphi + rand_dphi)*r_circ
    else:
        xc[i] = np.cos(i*dphi)*r_circ 
        yc[i] = np.sin(i*dphi)*r_circ
    xc[i] += np.random.normal()*extra_jitter
    yc[i] += np.random.normal()*extra_jitter

#Add these telescopes to a UV plane.
uv_tel = np.zeros((sz,sz))
uv_tel[(xc/fov_bl).astype(int) + sz//2,(yc/fov_bl).astype(int) + sz//2]=1

#Create a uv plane from this
uv = np.fft.fftshift(np.abs(np.fft.ifft2(np.abs(np.fft.fft2(uv_tel))**2)))
a = sig.convolve(uv,[[0.5,1,0.5],[1,1,1],[0.5,1,0.5]], mode='same')

#Now create a uv plane with sky rotation.
a0 = a.copy()
for i in range(1,nrot):
    a += nd.interpolation.rotate(a0, rotdeg*i, reshape=False, order=1)

a1 = a.copy()
for i in range(1,len(zooms)):
    az = nd.interpolation.zoom(a1, zooms[i], order=1)
    sz0 = az.shape[0]
    a += az[sz0//2-512:sz0//2+512,sz0//2-512:sz0//2+512]

f1 = plt.figure(1)
f1.clf()
ax1 = f1.add_subplot(111, aspect='equal')

plt.imshow(np.log(np.maximum(a,.5)), cmap=cm.gray, interpolation='nearest')

ax1.set_xticks([])
ax1.set_yticks([])


f3 = plt.figure(3)
f3.clf()
ax = f3.add_subplot(111, aspect='equal')
im_psf = np.fft.fft2(a)
im_psf = np.abs(im_psf)**2
im_psf /= np.max(im_psf)

nth=100
thetas = np.arange(nth)/nth*2*np.pi
imring = np.zeros(im_psf.shape)
for theta in thetas:
    xd = int(sz//2 + 0.8*ring_rad/(fov/sz)*np.cos(theta))
    yd = int(sz//2 + ring_rad/(fov/sz)*np.sin(theta))
    imring[xd,yd] += ring_flux/nth

#Convolve the ring image with the PSF.
imring = np.fft.irfft2(np.fft.rfft2(imring)*np.fft.rfft2(np.fft.fftshift(im_psf)))

#For a planet of circumstellar disk magnitude 11.5 and star of magnitude 8.5, we 
#just add a copy of the PSF and add a copy of the ring.
im = im_psf.copy()
im += np.roll(np.roll(im_psf,int(0.035/fov*sz/np.sqrt(2)), axis=1),int(0.035/fov*sz/np.sqrt(2)), axis=0)*10**(-0.4*3.0)
im += imring

iml = np.log10(np.maximum(im,1e-4))
plt.imshow(np.fft.fftshift(iml), interpolation='nearest', cmap=cm.gist_heat, extent=[-fov/2,fov/2,-fov/2,fov/2])
#plt.axis([512-200,512+200,512-200,512+200])
plt.colorbar()
ax.set_xlabel('Delta RA (arcsec)')
ax.set_ylabel('Delta Dec (arcsec)')


f2 = plt.figure(2)
f2.clf()
ax = f2.add_subplot(111, aspect='equal')
ax.plot(yc,xc,'o')
ax.set_xlabel('X position (m)')
ax.set_ylabel('Y position (m)')

dists1 = np.sqrt( (xc[1:]-xc[:-1])**2 + (yc[1:]-yc[:-1])**2)
dists2 = np.sqrt( (xc[2:]-xc[:-2])**2 + (yc[2:]-yc[:-2])**2)
dists10 = np.sqrt( (xc[10:]-xc[:-10])**2 + (yc[10:]-yc[:-10])**2)
print(np.min(dists1))
print(np.max(dists2))

f = open(outfile, "w")
for i in range(ntel):
    f.write("{0:6.1f} {1:6.1f}\n".format(xc[i],yc[i]))
f.close()

plt.show()

#Save the circular array.
np.savetxt('circ_array.txt',np.array( [xc,yc] ).T,fmt="%.1f")