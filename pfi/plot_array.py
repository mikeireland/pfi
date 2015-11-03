"""This plots an example 30 telescope array, in a near-Y configuration"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.signal as sig

dl = 82.2
dphi = np.radians(10)
l = 2**(np.arange(10)*0.5) * dl
xc = np.zeros(10)
yc = np.zeros(10)
phi = 0
xc[0] = dl/np.sqrt(3)
for i in range(1,10):
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

print "Max baseline: ", np.sqrt(((xc[9]-xc[19])**2 + (yc[9]-yc[19])**2))

uv = np.zeros((1024,1024))
uv[(xc/35).astype(int) + 512,(yc/35).astype(int) + 512]=1
test = np.fft.fftshift(np.abs(np.fft.ifft2(np.abs(np.fft.fft2(uv))**2)))
a = sig.convolve(test,[[0.5,1,0.5],[1,1,1],[0.5,1,0.5]])
plt.figure(1)
plt.clf()

plt.imshow(np.minimum(a,1), cmap=cm.gray, interpolation='nearest')
im = np.fft.fft2(test)
im = np.abs(im)**2
plt.imshow(np.fft.fftshift(np.minimum(im,410000)), interpolation='nearest')

f = plt.figure(2)
f.clf()
ax = f.add_subplot(111, aspect='equal')
ax.plot(yc,xc,'o')
ax.set_xlabel('X position (m)')
ax.set_ylabel('Y position (m)')
plt.show()