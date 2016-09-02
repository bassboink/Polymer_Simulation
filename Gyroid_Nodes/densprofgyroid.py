#!/usr/bin/python

# Script:  densprofgyroid.py
# Purpose: To find the density profile of gyroid systems by g value
# Syntax:  densprofgyroid.py < nodes.lammpstrj (dump file with scaled coordinates and nodes and tubes separated by type with gyroid_nodes_chains.py)
# Author:  Mitchell Wendt, derived from Youngmi Seo, derived from Mark Stevens' g(r) code
# -------------------------------------------------------------------------

import sys,string
from numpy import *
from math import * 

## INPUT PARAMETERS
nconf =   9 # number of configurations to average over
iskip =   0 #number of configurations to skip before taking data
mg    =    120  #sets number of bins


file  = 'densprof.txt'

def gyroid(x, y, z):
    x = x/a   
    y = y/a
    z = z/a
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

    maxg = max(g1, g2)

    return maxg

def den():
    for ii in range(1,natoms+1):
        gr = gyroid(xc[ii], yc[ii], zc[ii])
        ig = int(gr / drg) + 1
        if typea[ii] == 4:
            d1[ig] = d1[ig] + 1.0
            if ig >= 117:
                print atomid[ii], " ", gr
        elif typea[ii] == 2:
            d2[ig] = d2[ig] + 1.0

# Initial Variables
natoms = 0
conf = 0

# Skip configurations
for i in range(0,iskip+1):  #skip configurations before getting data
    sys.stdin.readline()
    sys.stdin.readline()      #time step
    sys.stdin.readline()
    line = sys.stdin.readline()      #number of atoms
    fields = string.split(line)
    conf = conf + 1
    if i==0:
        natoms = int(fields[0])  #reads number of atoms from first configuration; don't support changing num atoms
        dim=natoms+1
        xc=zeros(dim,float32)
        yc=zeros(dim,float32)
        zc=zeros(dim,float32)
        cx=zeros(dim)
        cy=zeros(dim)
        cz=zeros(dim)
        atomid=[0]*dim
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
    vol = xbox*ybox*zbox
    for j in range(1,dim):
        line = sys.stdin.readline()
        #[ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)
        
    # den() prefactors
    drg = zbox/mg
    numg = mg + 1
 
print "Reading config file...."

# Read configurations to average out
d1 = zeros(numg,float32)  
d2 = zeros(numg,float32)
istart = iskip+1
for kconf in range(istart,nconf+istart): 
    sys.stdin.readline()
    sys.stdin.readline()             # time step
    sys.stdin.readline()
    line = sys.stdin.readline()      # number of atoms
    fields = string.split(line)
    num = fields[0]
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
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3,n1,n2,n3] = string.split(line)
        k=int(ii)
        atomid[k] = k
        typea[k] = int(typej)
        mol[k] = int(molj)
        xc[k] = xbox*float(x1)         #scaled coords go from 0 to 1 in dump file
        yc[k] = ybox*float(x2)
        zc[k] = zbox*float(x3)


    # den() prefactors
    drg = 15./mg
    numg = mg + 1
    lbox = vol**(1./3.)
    a = lbox/(2*pi)             # Parameter deciding the gyroid lattice spacing 

    # call den()
    den()
    
    conf = conf + 1
    dd1 = d1 / nconf
    dd2 = d2 / nconf

    ## OUTPUT
    OUT = open(file, 'w')
    OUT.write("#%7i\n" % (conf))
    OUT.write("z den(A) den(B)\n")
    for ig in range(1,numg):
        OUT.write("%8.4f %8.4f %8.4f\n" % (ig,dd1[ig],dd2[ig]))
    OUT.close()

print "Done."




