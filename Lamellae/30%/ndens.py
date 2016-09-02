#Python coding for reading Lz, ndens
#Made by Youngmi Seo
#110415
#python ndens.py < equil.lammpstrj
#-------------------------------------------

import sys,string
from numpy import *
from math import *
from random import *
from itertools import chain

#Original Sample
nconf = 50
iskip = 449
nlayers = 4

for i in range(0,iskip+1):
    #print i
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
        Lz = 0
        ndens = 0
        avgLz = 0
        avgndens = 0
        
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
    Lz = zbox/nlayers
    ndens = natoms/(xbox*ybox*zbox)
    print i, Lz, ndens
    vol = xbox*ybox*zbox
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3,n1,n2,n3,v1,v2,v3] = string.split(line)

for i in range(iskip+1,nconf+iskip+1):
    #print i
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
    xbox = xp - xm
    ybox = yp - ym
    zbox = zp - zm
    vol = xbox*ybox*zbox
    hx2 = xbox/2
    hy2 = ybox/2
    hz2 = zbox/2
    Lz = zbox/nlayers
    ndens = natoms/(xbox*ybox*zbox)
    print i, Lz, ndens
    avgLz += Lz/nconf
    avgndens += ndens/nconf
    for ii in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3,n1,n2,n3,v1,v2,v3] = string.split(line)

print(avgLz,avgndens)

    
    

            
