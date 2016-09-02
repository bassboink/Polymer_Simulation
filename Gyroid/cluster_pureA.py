#Python coding for classifying two segments in double gyroid for different coloring in VMD 
#Made by Youngmi Seo
#060513
#python cluster.py < equil2.lammpstrj
#-------------------------------------------

import sys,string
from numpy import *
from math import *
from random import *
from itertools import chain

#Original Sample
nconf = 1
npoly = 1000
nbeads = 40
fA = 0.325
taper = 0
npureB = nbeads*(1-fA-taper/2)
npureA = nbeads*(fA - taper/2)
deltr = 1.5    #should not exceed the shortest distance btw the two gyroids
initAtom = 2   #random atom from pure A

INPUT_LAMMPS = open('clustered_pswap.lammpstrj', 'w')

for i in range(0,nconf):
    print i
    if i==0:
        conf = 0
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
            xc[k] = (float(x1)-0.5)*xbox
            yc[k] = (float(x2)-0.5)*ybox
            zc[k] = (float(x3)-0.5)*zbox

        init = initAtom
        collect = [init]
        collect1 = []
        next1 = [[0 for x in range(0)] for x in range(100)]   #set number of list at most
        next1[0].append(init)

        for j in [x for x in range(1,dim) if x not in collect]:
            if (j % nbeads) <= (npureA) and (j % nbeads) != 0:  #pure A blocks
                idist = sqrt(((xc[init]-xc[j])-xbox*round((xc[init]-xc[j])/xbox))**2+((yc[init]-yc[j])-ybox*round((yc[init]-yc[j])/ybox))**2+((zc[init]-zc[j])-zbox*round((zc[init]-zc[j])/zbox))**2)
                if idist < deltr:
                    collect.append(j)
                    next1[1].append(j)
            
        count = 0
        while len(next1[count+1]) != 0: #len(collect) <= npureB*npoly*0.5:
            count += 1
            for i in next1[count]:
                for j in [x for x in range(1,dim) if x not in collect]:
                    if (j % nbeads) <= (npureA) and (j % nbeads) != 0:  #pure A blocks
                        idist = sqrt(((xc[i]-xc[j])-xbox*round((xc[i]-xc[j])/xbox))**2+((yc[i]-yc[j])-ybox*round((yc[i]-yc[j])/ybox))**2+((zc[i]-zc[j])-zbox*round((zc[i]-zc[j])/zbox))**2)
                        if idist < deltr:
                            collect.append(j)
                            next1[count+1].append(j) 

            collect1[:] = [x-1 for x in collect]
            print '\n'
            print ' '.join(str(x) for x in collect1)        

        for ii in [x for x in collect]:
            typea[ii] = 3
        for jj in [x for x in range(1,dim) if x not in collect]:
            if (jj % nbeads) <= (npureA) and (jj % nbeads) != 0:  #pure A blocks
                typea[jj] = 4

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
            xc[k] = (float(x1)-0.5)*xbox
            yc[k] = (float(x2)-0.5)*ybox
            zc[k] = (float(x3)-0.5)*zbox

           
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
