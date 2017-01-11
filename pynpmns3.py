# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 15:54:02 2016

@author: ogurcan
"""
import numpy as np
import h5py as h5
import matplotlib.pyplot as plt
from scipy.integrate import ode
flname="nsout.h5"
fl=h5.File(flname,"w")
grp=fl.create_group("fields")

#dset=g.create_dataset("init",data)
#Shell[0] is assumed to be a dodecahedron.
sh0=1
#i.e. sh0=1 means it is an icosahedron

links_ico=np.array([
    [[[10,10],[11,11],[7,12],[8,13],[9,14]]
    ,[[5,10],[6,11],[7,12],[8,13],[9,14]]
    ,[[15,10],[16,11],[17,7],[18,8],[19,9]]],
    [[[3,10],[4,14],[11,15],[6,7],[8,19]]
    ,[[1,10],[3,14],[18,15],[12,7],[16,19]]
    ,[[11,3],[13,4],[8,11],[2,6],[6,8]]],
    [[[5,10],[4,11],[9,15],[7,16],[6,8]]
    ,[[4,10],[2,11],[17,15],[19,16],[13,8]]
    ,[[14,5],[12,4],[7,9],[9,7],[3,6]]],
    [[[1,11],[5,12],[10,16],[8,17],[6,9]]
    ,[[0,11],[3,12],[18,16],[15,17],[14,9]]
    ,[[10,1],[13,5],[8,10],[5,8],[4,6]]],
    [[[6,5],[1,13],[2,12],[9,18],[11,17]]
    ,[[10,5],[4,13],[1,12],[16,18],[19,17]]
    ,[[0,6],[14,1],[11,2],[6,9],[9,11]]],
    [[[6,6],[7,18],[2,14],[3,13],[10,19]]
    ,[[11,6],[15,18],[0,14],[2,13],[17,19]]
    ,[[1,6],[5,7],[10,2],[12,3],[7,10]]],
    [[[4,0],[5,1],[1,2],[2,3],[3,4]]
    ,[[15,0],[16,1],[17,2],[18,3],[19,4]]
    ,[[5,4],[6,5],[7,1],[8,2],[9,3]]],
    [[[9,0],[10,4],[5,5],[0,17],[2,9]]
    ,[[11,0],[13,4],[8,5],[2,17],[6,9]]
    ,[[1,9],[3,10],[18,5],[12,0],[16,2]]],
    [[[11,0],[10,1],[3,5],[1,6],[0,18]]
    ,[[14,0],[12,1],[7,5],[9,6],[3,18]]
    ,[[4,11],[2,10],[17,3],[19,1],[13,0]]],
    [[[7,1],[11,2],[4,6],[2,7],[0,19]]
    ,[[10,1],[13,2],[8,6],[5,7],[4,19]]
    ,[[0,7],[3,11],[18,4],[15,2],[14,0]]],
    [[[0,15],[7,3],[8,2],[3,8],[5,7]]
    ,[[0,15],[14,3],[11,2],[6,8],[9,7]]
    ,[[10,0],[4,7],[1,8],[16,3],[19,5]]],
    [[[0,16],[1,8],[8,4],[9,3],[4,9]]
    ,[[1,16],[5,8],[10,4],[12,3],[7,9]]
    ,[[11,0],[15,1],[0,8],[2,9],[17,4]]]
    ])

links_dod=np.array([
    [[[15,6],[11,7],[14,8]]
    ,[[4,6],[9,7],[11,8]]
    ,[[10,15],[3,11],[5,14]]],
    [[[16,6],[12,8],[10,9]]
    ,[[5,6],[10,8],[7,9]]
    ,[[11,16],[4,12],[1,10]]],
    [[[17,6],[13,9],[11,10]]
    ,[[1,6],[11,9],[8,10]]
    ,[[7,17],[5,13],[2,11]]],
    [[[18,6],[14,10],[12,11]]
    ,[[2,6],[7,10],[9,11]]
    ,[[8,18],[1,14],[3,12]]],
    [[[19,6],[13,7],[10,11]]
    ,[[3,6],[10,7],[8,11]]
    ,[[9,19],[4,13],[2,10]]],
    [[[8,7],[7,8],[10,4]]
    ,[[5,7],[3,8],[6,4]]
    ,[[11,8],[9,7],[0,10]]],
    [[[9,8],[8,9],[11,5]]
    ,[[1,8],[4,9],[6,5]]
    ,[[7,9],[10,8],[0,11]]],
    [[[12,1],[5,9],[9,10]]
    ,[[6,1],[2,9],[5,10]]
    ,[[0,12],[8,5],[11,9]]],
    [[[13,2],[6,10],[5,11]]
    ,[[6,2],[3,10],[1,11]]
    ,[[0,13],[9,6],[7,5]]],
    [[[6,7],[14,3],[7,11]]
    ,[[2,7],[6,3],[4,11]]
    ,[[8,6],[0,14],[10,7]]],
    [[[5,0],[1,1],[4,2]]
    ,[[10,0],[3,1],[5,2]]
    ,[[4,5],[9,1],[11,4]]],
    [[[6,0],[2,2],[0,3]]
    ,[[11,0],[4,2],[1,3]]
    ,[[5,6],[10,2],[7,0]]],
    [[[7,0],[3,3],[1,4]]
    ,[[7,0],[5,3],[2,4]]
    ,[[1,7],[11,3],[8,1]]],
    [[[8,0],[4,4],[2,5]]
    ,[[8,0],[1,4],[3,5]]
    ,[[2,8],[7,4],[9,2]]],
    [[[9,0],[3,1],[0,5]]
    ,[[9,0],[4,1],[2,5]]
    ,[[3,9],[10,3],[8,0]]],
    [[[18,1],[17,2],[0,10]]
    ,[[11,1],[9,2],[0,10]]
    ,[[5,18],[3,17],[6,0]]],
    [[[19,2],[18,3],[1,11]]
    ,[[7,2],[10,3],[0,11]]
    ,[[1,19],[4,18],[6,1]]],
    [[[2,7],[15,3],[19,4]]
    ,[[0,7],[8,3],[11,4]]
    ,[[6,2],[2,15],[5,19]]],
    [[[3,8],[16,4],[15,5]]
    ,[[0,8],[9,4],[7,5]]
    ,[[6,3],[3,16],[1,15]]],
    [[[16,1],[4,9],[17,5]]
    ,[[8,1],[0,9],[10,5]]
    ,[[2,16],[6,4],[4,17]]]
    ])

def icosp():
    theta=np.hstack([0, np.tile(np.pi/2-np.arctan(1/2),5),np.pi,np.tile(np.pi/2+np.arctan(1/2),5)]);
    phi=np.hstack([0,np.arange(0,2*np.pi,2*np.pi/5),0,np.mod(np.arange(0,2*np.pi,2*np.pi/5)+np.pi,2*np.pi)]);
    return theta,phi;
    
def dodsp():
    ph=(1+np.sqrt(5))/2;
    alp=np.array([np.arcsin(ph/np.sqrt(3))-np.arccos(ph/np.sqrt(ph+2)),np.arctan(2*ph**2)]);
    eta=np.arange(np.pi/5,2*np.pi,2*np.pi/5);
    theta=np.reshape(np.transpose(np.reshape(np.transpose(np.tile(np.hstack([alp,np.pi-alp]),[5,1])),[20,1])),20);
    phi=np.reshape(np.hstack([np.tile(eta,[1,2]),np.tile(eta-np.pi,[1,2])]),20);
    return theta,phi;

def nbase(n):
    return (n%2)*4*int(2*(1/2-sh0))+n*16

def symmetrize(un):
    for n in np.arange(0,Nmax):
        npoly=20-(n+sh0)%2*8
        for l in np.arange(int(npoly/2),npoly):
            un[:,nbase(n)+l]=un[:,nbase(n)+l-int(npoly/2)].conj()
    return un;

def get_nls(n,l):
    if(n%2==sh0): # i.e. if zero or even
        nlp=np.zeros((2,9,2),dtype=int)
        nlp[0,0:3,0]=n-2
        nlp[0,0:3,1]=n-1
        nlp[0,3:6,0]=n-1
        nlp[0,3:6,1]=n+1
        nlp[0,6:9,0]=n+1
        nlp[0,6:9,1]=n+2
        nlp[1,:,:]=links_dod[l,:,:,:].reshape(nlp[1,:,:].shape)
    else:
        nlp=np.zeros((2,15,2),dtype=int)
        nlp[0,0:5,0]=n-2
        nlp[0,0:5,1]=n-1
        nlp[0,5:10,0]=n-1
        nlp[0,5:10,1]=n+1
        nlp[0,10:15,0]=n+1
        nlp[0,10:15,1]=n+2
        nlp[1,:,:]=links_ico[l,:,:,:].reshape(nlp[1,:,:].shape)
    return nlp

def get_interacting_nodes(n,l):
    npl=get_nls(n,l)
    np=nbase(npl[0,:,0])
    lp=npl[1,:,0]
    npp=nbase(npl[0,:,1])
    lpp=npl[1,:,1]
    nlp=np+lp
    nlpp=npp+lpp
    return nlp,nlpp

Nmax=40
NLmax=Nmax*16
np.floor((Nmax/2)).astype(int)*2*16+(Nmax%2)*(20-((Nmax+1+sh0)%2)*8)
#k0=2.0**(-4)
k0=1.0
g=np.sqrt((1+np.sqrt(5))/2)
lam=np.sqrt(np.sqrt(5)/3)
th,ph=dodsp()
khatd=np.array([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])
th,ph=icosp()
khati=np.array([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)])
k=np.zeros((3,NLmax));
links=[]
for n in np.arange(0,Nmax):
    if(n%2==sh0):
        k[:,nbase(n):nbase(n+1)]=khatd*g**n
    else:
        k[:,nbase(n):nbase(n+1)]=khati*g**n*lam

    for l in np.arange(0,20-(n+sh0)%2*8):
        nl=nbase(n)+l
        nlp,nlpp=get_interacting_nodes(n,l)
        inds=np.nonzero((nlp>=0) & (nlpp>=0) & (nlp<NLmax) & (nlpp<NLmax))
        links.append([nlp[inds],nlpp[inds],inds]);
unls=np.zeros((3,NLmax),dtype=complex)
kn=np.sqrt(np.sum(k**2,0))
unls[0,:]=0.0001+0.0001j;#1e-6*np.exp(-1e-8*kn**2+1j*2*np.pi*np.random.random((NLmax,)));
unls[1,:]=0.0001+0.0001j;#1e-6*np.exp(-1e-8*kn**2+1j*2*np.pi*np.random.random((NLmax,)));
unls[2,:]=0.0001+0.0001j;#1e-6*np.exp(-1e-8*kn**2+1j*2*np.pi*np.random.random((NLmax,)));
unls=unls-(np.einsum("ij,ij->j",k,unls)*k)/kn**2
unls=symmetrize(unls)
Mkpq=np.einsum("kl,ij->kijl",k,np.eye(3))+np.einsum("jl,ik->kijl",k,np.eye(3))-2*np.einsum("il,jl,kl,l->kijl",k,k,k,1/kn**2)

uu0=unls.reshape(unls.shape[0]*unls.shape[1]);
t0=0.0;
t1=100.0;
dt=0.1;
nu=1e-13;
dudt=np.zeros(np.shape(unls),dtype=unls.dtype);
Nt=int(t1/dt)+1;

Dn=nu*kn**4
Fn=np.zeros((3,NLmax),dtype=unls.dtype);
Fn[:,nbase(3):nbase(4)]=0.01+0.01j
Fn[:,nbase(4):nbase(5)]=0.01-0.01j
Fn=Fn-(np.einsum("ij,ij->j",k,Fn)*k)/kn**2
Fn=symmetrize(Fn)

def func(t,y,par):
    u=y.reshape((3,NLmax));
    for n in np.arange(0,NLmax):
        nlp=links[n][0]
        nlpp=links[n][1]
        up=u[:,nlp]
        upp=u[:,nlpp]
	dudt[:,n]=Fn[:,n]-Dn[n]*u[:,n]-1j*np.einsum("kij,kl,jl->i",Mkpq[:,:,:,n],up.conj(),upp.conj())
        dydt=np.reshape(dudt,NLmax*3);
    return dydt;

grp.create_dataset("k",data=k)
ures=grp.create_dataset("u",(1,3,NLmax),maxshape=(None,3,NLmax),dtype=complex)
r = ode(func, 0).set_integrator('zvode', method='bdf', with_jacobian=False,atol=1e-14,rtol=1e-8,nsteps=1E6);
r.set_initial_value(uu0, t0).set_f_params(0.0).set_jac_params(0.0);
i=0;
while r.successful() and r.t < t1:
    print("t=",r.t);
    r.integrate(r.t+dt)
    ures=fl["/fields/u"]
    ures.resize((i+1,3,NLmax))
    ures[i,:,:]=np.reshape(r.y,(3,NLmax))
    fl.flush()
    i=i+1;
