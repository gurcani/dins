import numpy as np
import matplotlib.pylab as plt
import h5py

fl=h5py.File('test.h5');
k=fl['/fields/k'].value[0]
u=fl['/fields/u'].value
u=u['real']+1j*u['imag']
Nt=6800
vk=u[Nt,:,:]
[xx,yy]=np.meshgrid(np.arange(-2*np.pi,2*np.pi,np.pi/100),np.arange(-2*np.pi,2*np.pi,np.pi/100))
vx=np.zeros(xx.shape,dtype=complex)
vy=np.zeros(xx.shape,dtype=complex)
for l in range(k.shape[1]):
    vx=vx+vk[0,l]*np.exp(1j*(k[0,l]*xx+k[1,l]*yy))
    vy=vy+vk[1,l]*np.exp(1j*(k[0,l]*xx+k[1,l]*yy))
E=np.abs(vx)**2+np.abs(vy)**2
plt.figure(dpi=200)
plt.pcolormesh(xx,yy,E,shading='gouraud',rasterized=True,vmin=0,vmax=5)
plt.xlim(xmin=-6,xmax=6)
plt.ylim(ymin=-6,ymax=6)
plt.colorbar()
plt.show()
#plt.savefig('out0.png', bbox_inches='tight')
