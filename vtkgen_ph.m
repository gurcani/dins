idd=[[0,4,3,2,1];[5,17,9,4,0];[1,6,18,5,0];[7,19,6,1,2];[8,15,7,2,3];[16,8,3,4,9];[10,11,12,13,14];[15,10,14,19,7];[11,10,15,8,16];[17,12,11,16,9];[18,13,12,17,5];[6,19,14,13,18]];
dii=[[0,1,2];[0,2,3];[0,3,4];[0,4,5];[0,5,1];[2,1,10];[2,11,3];[3,7,4];[4,8,5];[5,1,9];[6,8,7];[6,9,8];[6,10,9];[6,11,10];[6,7,11];[8,4,7];[8,9,5];[9,10,1];[10,11,2];[11,7,3]];
iidd=[[0,1,0,4];[0,2,1,0];[0,3,2,1];[0,4,3,2];[0,5,4,3];[1,2,0,5];[1,5,9,4];[1,9,17,9];[1,10,5,17];[2,10,18,5];[2,11,6,18];[2,3,1,6];[3,4,2,7];[3,7,7,19];[3,11,19,6];[4,5,3,8];[4,8,8,15];[4,7,15,7];[5,8,16,8];[5,9,9,16];[6,7,14,10];[6,11,13,14];[6,10,12,13];[6,9,11,12];[6,8,10,11];[7,8,15,10];[7,11,14,19];[8,9,16,11];[9,10,17,12];[10,11,18,13]];

fl=load('test.h5');
k=fl.fields.k;
u=fl.fields.u;
fid = fopen ("out_ph.vtk", "w");
fprintf (fid,"# vtk DataFile Version 3.0\n");
fprintf (fid,"patates\n");
fprintf (fid,"ASCII\nDATASET UNSTRUCTURED_GRID\n");
fprintf (fid,"POINTS %i FLOAT\n",size(k,1));

for l=1:size(k,1)
  kx=k(l,1);
  ky=k(l,2);
  kz=k(l,3);
  kk=sqrt(kx.^2+ky.^2+kz.^2);
  th=acos(kz/kk);
  if(kx==0)
    ph=0;
  else
    ph=atan2(ky,kx);
  end
  kkx=log10(kk)*sin(th)*cos(ph);
  kky=log10(kk)*sin(th)*sin(ph);
  kkz=log10(kk)*cos(th);
  fprintf (fid, "% 4.4f % 4.4f % 4.4f\n",kkx,kky,kkz);
end
%first icosahedron:
Nico=12;
Ndod=20;
Nedge=30;
Ns=Nico+Nedge+Ndod;
Nsdat=Ndod*7+Nico*11+Nedge*7;
Nshells=9;
fprintf (fid,"CELLS %i %i\n",Ns*Nshells,Nsdat*Nshells);

for n=1:Nshells
  for l=1:Ndod
    d0=l-1+Nico+(n-1)*(Nico+Ndod);
    ii=dii(l,:)+n*(Nico+Ndod);
    i0=ii(1);i1=ii(2);i2=ii(3);
    fprintf (fid, "6 %i %i %i %i %i %i \n",d0,i0,i1,i2,d0,i0);
  end
  for l=1:Nico
    i0=l-1+n*(Nico+Ndod);
    dd=idd(l,:)+Nico+(n-1)*(Nico+Ndod);
    d0=dd(1);d1=dd(2);d2=dd(3);d3=dd(4);d4=dd(5);
    fprintf (fid, "10 %i %i %i %i %i %i %i %i %i %i\n",d1,i0,d2,d3,d4,i0,d0,d1,d4,d2);
  end
  for l=1:Nedge
    i0=iidd(l,1)+n*(Nico+Ndod);i1=iidd(l,2)+n*(Nico+Ndod);
    d0=iidd(l,3)+Nico+(n-1)*(Nico+Ndod);d1=iidd(l,4)+Nico+(n-1)*(Nico+Ndod);
    fprintf (fid, "6 %i %i %i %i %i %i\n",i0,i1,d0,d1,i0,i1);
  end
end
fprintf (fid,"CELL_TYPES %i\n",Ns*Nshells);
for l=1:Ns*Nshells
  fprintf (fid, "6\n");
end

fprintf(fid,"POINT_DATA %i\n",size(k,1));
fprintf(fid,"SCALARS sample_scalars float 1\n");
fprintf(fid,"LOOKUP_TABLE default\n");
fprintf(fid,"%f\n",arg(u(:,1,end)));
fclose(fid);
