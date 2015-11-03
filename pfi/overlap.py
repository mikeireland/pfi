import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

x = 0.01*((np.arange(1000)+0.5) - 500)
y = 2*(sp.jn(1,x*pi)/(x*pi))
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
hw = 60
x = np.arange(2*hw)-hw
xy = np.meshgrid(x,x)
p = (xy[0]**2 + xy[1]**2) < (hw/12.0)**2
im = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(p)))
win = np.zeros((2*hw,2*hw))
win[hw-12:hw+12,hw-24:hw+24]=1
print "Max Star Throughput: ", np.sum(win*np.abs(im)**2)/np.sum(np.abs(im)**2)
print "Min Star Throughput: ", np.sum(np.abs(np.roll(im,8,axis=0))**2*win)/np.sum(np.abs(im)**2)
overlap = np.sum(np.real(np.roll(im,8,axis=0)*np.conj(im)*win))/np.sum(win*np.abs(im)**2)
print "Min Star Overlap: ", overlap
print "Effective min Throughput: ", overlap**2 