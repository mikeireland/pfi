"""This plots an example 30 telescope array, in a near-Y configuration

It is a script - so change the numbers in the Settings section to
change the array parameters."""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.signal as sig
import scipy.ndimage as nd


#---- Settings ----
#Number of telescopes divided by 3
nt_3 = 7

#Length difference between successive elements in m
dl = 130.0

#Angular 
dphi = 0.0 #np.radians(10.0)

#Logarithmic distance between successive telescopes
l = 2**(np.arange(nt_3)*0.5) * dl

#---- More or less automatic from here ----
xc = np.zeros(nt_3)
yc = np.zeros(nt_3)
phi = 0
xc[0] = dl/np.sqrt(3)
for i in range(1,nt_3):
  xc[i] = xc[i-1] + l[i-1]*np.cos(phi)
  yc[i] = yc[i-1] + l[i-1]*np.sin(phi)
  phi += dphi
xy = np.array([xc,yc])  
R = np.array([[-0.5,np.sqrt(3)/2],[-np.sqrt(3)/2,-0.5]])
xy1 = np.dot(R,xy)
xc = np.append(xc,xy1[0])
yc = np.append(yc,xy1[1])
xy1 = np.dot(np.linalg.inv(R),xy)
xc = np.append(xc,xy1[0])
yc = np.append(yc,xy1[1])

print "Max baseline: ", np.sqrt(((xc[nt_3-1]-xc[2*nt_3-1])**2 + (yc[nt_3-1]-yc[2*nt_3-1])**2))

uv = np.zeros((1024,1024))
uv[(xc/35).astype(int) + 512,(yc/35).astype(int) + 512]=1
test = np.fft.fftshift(np.abs(np.fft.ifft2(np.abs(np.fft.fft2(uv))**2)))
a = sig.convolve(test,[[0.5,1,0.5],[1,1,1],[0.5,1,0.5]])

plt.figure(1)
plt.clf()
plt.imshow(np.minimum(a,1), cmap=cm.gray, interpolation='nearest')

plt.figure(2)
plt.clf()
im = np.fft.fft2(test)
im = np.abs(im)**2
im = np.fft.fftshift(im)
im0 = im.copy()
#Add some sky rotation
for i in range(30):
    im += nd.interpolation.rotate(im0,i*2,reshape=False)
im /= np.max(im)
plt.imshow(im, interpolation='nearest')
plt.colorbar()

f = plt.figure(3)
f.clf()
ax = f.add_subplot(111, aspect='equal')
ax.plot(yc,xc,'o')
ax.set_xlabel('X position (m)')
ax.set_ylabel('Y position (m)')
plt.show()

#Save the Y-shaped array.
np.savetxt('Y_array.txt',np.array( [xc,yc] ).T,fmt="%.1f")