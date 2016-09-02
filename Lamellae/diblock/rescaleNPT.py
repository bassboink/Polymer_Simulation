#Python code for rescaling molecular coordinates to fit in the averaged box size after NPT equilibration
#Made by Youngmi Seo
#060513
#python rescaleNPT.py < equil.lammpstrj
#-------------------------------------------

import sys,string
from numpy import *
from math import *
from random import *
from itertools import chain

#Original Sample
iskip = 50
Lz = 17.652658000000002
rho = 0.8422596736868155
npoly = 480
nbeads = 40
nmon = npoly*nbeads

#Output file name
file = 'rescaledinput.txt'

for i in range(0,iskip):
    print i
    sys.stdin.readline()
    sys.stdin.readline()
    sys.stdin.readline()
    line = sys.stdin.readline()
    fields = string.split(line)

    if i==0:
        natoms = int(fields[0])  #reads number of atoms from first configuration; don't support changing num atoms
        dim=natoms+1
        xc=zeros(dim,float32)
        yc=zeros(dim,float32)
        zc=zeros(dim,float32)
        ix=zeros(dim)
        iy=zeros(dim)
        iz=zeros(dim)
        typea=[0]*dim
        mol=[0]*dim
        conf = 0
        
    sys.stdin.readline()
    line = sys.stdin.readline()      
    [xm,xp] = map(float,line.split())
    line = sys.stdin.readline()      
    [ym,yp] = map(float,line.split())
    line = sys.stdin.readline()      
    [zm,zp] = map(float,line.split())
    line = sys.stdin.readline()
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)

        
for i in range(0,1):
    print i
    sys.stdin.readline()
    sys.stdin.readline()
    sys.stdin.readline()
    line = sys.stdin.readline()
    fields = string.split(line)        
    sys.stdin.readline()
    line = sys.stdin.readline()      
    [xm,xp] = map(float,line.split())
    line = sys.stdin.readline()      
    [ym,yp] = map(float,line.split())
    line = sys.stdin.readline()      
    [zm,zp] = map(float,line.split())
    line = sys.stdin.readline()
    nzbox = Lz*4                  #new box size at equilibrium state
    nxbox = sqrt(natoms/rho/nzbox)
    nybox = sqrt(natoms/rho/nzbox)
    xbox = xp - xm
    ybox = yp - ym
    zbox = zp - zm
    ndens = natoms/(xbox*ybox*zbox)
    nxbox2 = nxbox/2
    nybox2 = nybox/2
    nzbox2 = nzbox/2
    vol = xbox*ybox*zbox
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)
        k = int(ii)
        typea[k] = int(typej)           
        mol[k] = int(molj)
        ix[k] = int(n1)
        iy[k] = int(n2)
        iz[k] = int(n3)
        xc[k] = float(x1)
        yc[k] = float(x2)
        zc[k] = float(x3)

        xc[k] = nxbox*(float(x1) + ix[k])      #Rescale to new box size
        yc[k] = nybox*(float(x2) + iy[k])
        zc[k] = nzbox*(float(x3) + iz[k])

        if xc[k] > nxbox:                       #Image Flags
            ix[k] = int(xc[k]/nxbox)
            xc[k] = xc[k] - ix[k]*nxbox - nxbox2
        elif xc[k] < 0:
            ix[k] = -int((-xc[k] + nxbox)/nxbox)
            xc[k] = xc[k] - ix[k]*nxbox - nxbox2
        else:
            ix[k] = 0
            xc[k] = xc[k] - nxbox2
        if yc[k] > nybox:
            iy[k] = int(yc[k]/nybox)
            yc[k] = yc[k] - iy[k]*nybox - nybox2
        elif yc[k] < 0:
            iy[k] = -int((-yc[k] + nybox)/nybox)
            yc[k] = yc[k] - iy[k]*nybox - nybox2
        else:
            iy[k] = 0
            yc[k] = yc[k] - nybox2
        if zc[k] > nzbox:
            iz[k] = int(zc[k] /nzbox)
            zc[k] = zc[k] - iz[k]*nzbox - nzbox2
        elif zc[k] < 0:
            iz[k] = -int((-zc[k] + nzbox)/nzbox)
            zc[k] = zc[k] - iz[k]*nzbox - nzbox2
        else:
            iz[k] = 0
            zc[k] = zc[k] - nzbox2

        #scaled coords go from 0 to 1; want to go from -xbox/2 to xbox/2
    print(nxbox2,nybox2,nzbox2)

                
    conf = conf + 1


    #OUTPUT
    OUT = open(file,'w')
    OUT.write("configuration number\n")
    OUT.write("#%7i\n" % (conf))
    OUT.write("ii molj typej x y z xi yi zi\n")
    for k in range(1,natoms+1):
        OUT.write("%6d %6d %2d %9.4f %9.4f %9.4f %6d %6d %6d\n" % (k,mol[k],typea[k],xc[k],yc[k],zc[k],ix[k],iy[k],iz[k]))
    OUT.write("Bonds\n")
    molk = 0
    bond = 0
    for k in range(1,nmon+1):
        molk = molk + 1
        if molk % nbeads != 0:
            bond = bond + 1
            if typea[k] != typea[k+1]:
                OUT.write("%8d %2d %8d %8d\n" % (bond,2,molk,molk+1))
            else:
                OUT.write("%8d %2d %8d %8d\n" % (bond,1,molk,molk+1))
                
    OUT.write("Masses\n")
    OUT.write("%d %1.1f\n" % (1,1.0))
    OUT.write("%d %1.1f\n" % (2,1.0))
    OUT.close()
    
    

            
