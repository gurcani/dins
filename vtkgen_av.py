import numpy as np
import matplotlib.pylab as plt
import h5py

fl=h5py.File('test.h5');
k=fl['/fields/k'].value[0]
u=fl['/fields/u'].value
u=u['real']+1j*u['imag']

fid = open ("out.vtk", "w");
fid.write("# vtk DataFile Version 3.0\n");
fid.write("patates\n");
fid.write("ASCII\nDATASET UNSTRUCTURED_GRID\n");
fid.write("POINTS %i FLOAT\n" % k.shape[1]);

N=k.shape[1]

for l in range(N):
    kx=k[0,l];
    ky=k[1,l];
    kz=k[2,l];
    kk=np.sqrt(kx**2+ky**2+kz**2);
    th=np.arccos(kz/kk);
    if(kx==0):
        ph=0;
    else:
        ph=np.arctan2(ky,kx);
    kkx=np.log10(kk)*np.sin(th)*np.cos(ph);
    kky=np.log10(kk)*np.sin(th)*np.sin(ph);
    kkz=np.log10(kk)*np.cos(th);
    fid.write("% 4.4f % 4.4f % 4.4f\n" %(kkx,kky,kkz));
    kn=np.sqrt(kkz**2+kky**2+kkz**2)

fid.write("CELLS %i %i\n"%(N,2*N));

for l in range(N):
    fid.write("1 %i\n"%l)

fid.write("CELL_TYPES %i\n"%N);
for l in range(N):
    fid.write("1\n");

fid.write("POINT_DATA %i\n"%N)
fid.write("SCALARS sample_scalars float 1\n");
fid.write("LOOKUP_TABLE default\n");
for l in range(N):
    fid.write("%f\n"%np.log10(np.mean(np.sum(np.abs(u[5000:6800,:,l])**2,1),0)/kn));
fid.close()
