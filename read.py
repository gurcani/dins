# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 14:00:36 2016

@author: ogurcan
"""

import numpy as np
import h5py as h5
import matplotlib.pylab as plt
f=h5.File("test.h5")
u=(f['fields/u'].value)['real']+1j*(f['fields/u'].value)['imag']
k=f['fields/k'].value
f.close()
N=np.shape(k)[2]
np1=12
np2=20
Ns=int(2*N/(np2+np1))

Nt=np.arange(5000,6800)
Ntt=6800
Nr=Nt[-1]-Nt[0]+1

kn=np.zeros(Ns)
Enav=np.zeros((Nr,Ns))
En=np.zeros(Ns)

for n in np.arange(0,Ns):
    n0=int((n%2)*np1+np.floor(n/2)*(np2+np1))
    n1=int(((n+1)%2)*np1+np.floor((n+1)/2)*(np2+np1))
    kn[n]=np.sqrt(np.sum(k[0,:,n0]**2,0))
    En[n]=np.mean(np.sum(np.abs(u[Ntt,:,n0:n1])**2,0),0)
    Enav[:,n]=np.mean(np.sum(np.abs(u[Nt,:,n0:n1])**2,1),1)
plt.loglog(kn,np.mean(Enav,0)/kn,kn,En/kn,'x-',kn,5e-2*kn**(-5/3))
plt.ylim([1e-12,1e-1])
plt.xlim([0.8,2e4])
plt.xlabel('k')
plt.ylabel('E(k)')
#plt.show()
plt.savefig("spectrum.svg")
