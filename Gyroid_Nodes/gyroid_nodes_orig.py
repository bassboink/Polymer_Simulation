#Python coding for classifying node segnebt in double gyroid for different coloring in VMD 
#Made by Kevin Shen
#071416
#python gyroid_nodes.py < equil.lammpstrj
#-------------------------------------------

import sys,string
from numpy import *
from math import *
from random import *
from itertools import chain

#Original Sample
nconf = 10
#iskip = 0
npoly = 1000
nbeads = 40
accRange = 14.8
inOrNot = [0]*npoly

INPUT_LAMMPS = open('nodes.lammpstrj', 'w')

def gyroid(x, y, z):
    x = x/a;    
    y = y/a;
    z = z/a;
    cx = cos(2*x)
    cy = cos(2*y)
    cz = cos(2*z)
    g1 = 10.0*(cos(x) * sin(y) + cos(y) * sin(z) + cos(z) * sin(x)) - 0.5*(cx*cy + cy*cz + cz*cx)

    x = -x
    y = -y
    z = -z
    cx = cos(2*x)
    cy = cos(2*y)
    cz = cos(2*z)
    g2 = 10.0*(cos(x) * sin(y) + cos(y) * sin(z) + cos(z) * sin(x)) - 0.5*(cx*cy + cy*cz + cz*cx)

    return max(g1, g2)

for i in range(0,nconf):
    print i
    if i==0:
        sys.stdin.readline()
        line = sys.stdin.readline()
        fields = string.split(line)
        timestep = int(fields[0])
        sys.stdin.readline()
        line = sys.stdin.readline()
        fields = string.split(line)
        natoms = int(fields[0])  #reads number of atoms from first configuration; don't support changing num atoms
        dim=natoms+1
        xc=zeros(dim,float32)
        yc=zeros(dim,float32)
        zc=zeros(dim,float32)
        xu=zeros(dim,float32)
        yu=zeros(dim,float32)
        zu=zeros(dim,float32)
        ix=zeros(dim)
        iy=zeros(dim)
        iz=zeros(dim)
        typea=[0]*dim
        mol=[0]*dim
        sys.stdin.readline()
        line = sys.stdin.readline()      
        [xm,xp] = map(float,line.split())
        line = sys.stdin.readline()      
        [ym,yp] = map(float,line.split())
        line = sys.stdin.readline()      
        [zm,zp] = map(float,line.split())
        line = sys.stdin.readline()
        xbox = xp - xm
        ybox = yp - ym
        zbox = zp - zm
        xbox2 = xbox/2
        ybox2 = ybox/2
        zbox2 = zbox/2
        a = xbox/(2*pi)
        for j in range(1,dim):
            line = sys.stdin.readline()
            [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)
            k = int(ii)
            mol[k] = int(molj)
            typea[k] = int(typej)
            ix[k] = int(n1)
            iy[k] = int(n2)
            iz[k] = int(n3)
            xu[k] = float(x1)
            yu[k] = float(x2)
            zu[k] = float(x3)
            xc[k] = float(x1)*xbox
            yc[k] = float(x2)*ybox
            zc[k] = float(x3)*zbox 

        for i in range(npoly):
            for j in range(nbeads):
                na = i * nbeads + (j+1)
                print na, i
                if typea[na] == 1:
                    if inOrNot[i] == 0 and gyroid(xc[na], yc[na], zc[na]) >= accRange:
                        print "yes"
                        inOrNot[i] = 1
                        for k in range(nbeads):
                            na1 = i * nbeads + (k+1)
                            if typea[na1] == 1 or typea[na1] == 4:
                                typea[na1] = 3
                    else:
                        typea[na] = 4

    else:
        sys.stdin.readline()
        line = sys.stdin.readline()
        fields = string.split(line)
        timestep = int(fields[0])
        sys.stdin.readline()
        sys.stdin.readline()     
        sys.stdin.readline()
        line = sys.stdin.readline()      
        [xm,xp] = map(float,line.split())
        line = sys.stdin.readline()      
        [ym,yp] = map(float,line.split())
        line = sys.stdin.readline()      
        [zm,zp] = map(float,line.split())
        line = sys.stdin.readline()
        xbox = xp - xm
        ybox = yp - ym
        zbox = zp - zm
        xbox2 = xbox/2
        ybox2 = ybox/2
        zbox2 = zbox/2
        for j in range(1,dim):
            line = sys.stdin.readline()
            [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)
            k = int(ii)
            mol[k] = int(molj)
            typea[k] = int(typej)
            ix[k] = int(n1)
            iy[k] = int(n2)
            iz[k] = int(n3)
            xu[k] = float(x1)
            yu[k] = float(x2)
            zu[k] = float(x3)

           
    INPUT_LAMMPS.write("ITEM: TIMESTEP\n")
    INPUT_LAMMPS.write("%i\n" % (timestep))
    INPUT_LAMMPS.write("ITEM: NUMBER OF ATOMS\n")
    INPUT_LAMMPS.write("%i\n" % (natoms))
    INPUT_LAMMPS.write("ITEM: BOX BOUNDS pp pp pp\n")
    INPUT_LAMMPS.write("%.4f %.4f\n" % (-xbox2,xbox2))
    INPUT_LAMMPS.write("%.4f %.4f\n" % (-ybox2,ybox2))
    INPUT_LAMMPS.write("%.4f %.4f\n" % (-zbox2,zbox2))
    INPUT_LAMMPS.write("ITEM: ATOMS id mol type xs ys zs ix iy iz \n")
    for i in range(1,dim):
        INPUT_LAMMPS.write("%i %i %i %.8f %.8f %.8f %i %i %i\n" % (i, mol[i], typea[i], xu[i], yu[i], zu[i], ix[i], iy[i], iz[i]))

INPUT_LAMMPS.close()
print "LAMMPS output completed."
