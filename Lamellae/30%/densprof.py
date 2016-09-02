#!/usr/bin/python

# Script:  density profile by densprof.py
# Purpose: First, find center of mass of type 1 molecules in lamellar structure by circling. The idea is to combine multiple layers into one layer before circling. With the center of mass as interface, calculate average density profiles averaging over some last configurations in dump file.
# Syntax:  densprof.py < test.lammpstrj (dump file with scaled coordinates)
# Author:  Youngmi Seo, derived from Mark Stevens' g(r) code
# -------------------------------------------------------------------------

import sys,string
from numpy import *
from math import * 

## INPUT PARAMETERS
nconf =   2 # number of configurations to average over
iskip =   0 #number of configurations to skip before taking data
mg    =    120  #sets number of bins


file  = 'densprof.txt'

def den():
    for ii in range(1,natoms+1):
        zr = zc[ii]
        if zr >= zcm:
            ig = int((zr - zcm) / drg) + 1
            if ig < len(d1_under):
                if typea[ii] == 1:
                    for ik in range(0, len(molsUnder50)):
                        if mol[ii] == molsUnder50[ik]:
                            d1_under[ig] = d1_under[ig] + 1.0
                            if (id[ii] % 40) <= 14:
                                d1_under_pure[ig] = d1_under_pure[ig] + 1.0
                    for ik in range(0, len(mols50)):
                        if mol[ii] == mols50[ik]:
                            d1_50[ig] = d1_50[ig] + 1.0
                            if (id[ii] % 40) <= 14:
                                d1_50_pure[ig] = d1_50_pure[ig] + 1.0
                    for ik in range(0, len(molsOver50)):
                        if mol[ii] == molsOver50[ik]:
                            d1_over[ig] = d1_over[ig] + 1.0
                            if (id[ii] % 40) <= 14:
                                d1_over_pure[ig] = d1_over_pure[ig] + 1.0
                elif typea[ii] == 2:
                    for ik in range(0, len(molsUnder50)):
                        if mol[ii] == molsUnder50[ik]:
                            d2_under[ig] = d2_under[ig] + 1.0
                            if (id[ii] % 40) >= 29:
                                d2_under_pure[ig] = d2_under_pure[ig] + 1.0
                    for ik in range(0, len(mols50)):
                        if mol[ii] == mols50[ik]:
                            d2_50[ig] = d2_50[ig] + 1.0
                            if (id[ii] % 40) >= 29:
                                d2_50_pure[ig] = d2_50_pure[ig] + 1.0
                    for ik in range(0, len(molsOver50)):
                        if mol[ii] == molsOver50[ik]:
                            d2_over[ig] = d2_over[ig] + 1.0
                            if (id[ii] % 40) >= 29:
                                d2_over_pure[ig] = d2_over_pure[ig] + 1.0
        elif zr < zcm:
            ig = int((zr+zbox - zcm) / drg) + 1
            if ig < len(d1_under):
                if typea[ii] == 1:
                    for ik in range(0, len(molsUnder50)):
                        if mol[ii] == molsUnder50[ik]:
                            d1_under[ig] = d1_under[ig] + 1.0
                            if (id[ii] % 40) <= 14:
                                d1_under_pure[ig] = d1_under_pure[ig] + 1.0
                    for ik in range(0, len(mols50)):
                        if mol[ii] == mols50[ik]:
                            d1_50[ig] = d1_50[ig] + 1.0
                            if (id[ii] % 40) <= 14:
                                d1_50_pure[ig] = d1_50_pure[ig] + 1.0
                    for ik in range(0, len(molsOver50)):
                        if mol[ii] == molsOver50[ik]:
                            d1_over[ig] = d1_over[ig] + 1.0
                            if (id[ii] % 40) <= 14:
                                d1_over_pure[ig] = d1_over_pure[ig] + 1.0
                elif typea[ii] == 2:
                    for ik in range(0, len(molsUnder50)):
                        if mol[ii] == molsUnder50[ik]:
                            d2_under[ig] = d2_under[ig] + 1.0
                            if (id[ii] % 40) >= 29:
                                d2_under_pure[ig] = d2_under_pure[ig] + 1.0
                    for ik in range(0, len(mols50)):
                        if mol[ii] == mols50[ik]:
                            d2_50[ig] = d2_50[ig] + 1.0
                            if (id[ii] % 40) >= 29:
                                d2_50_pure[ig] = d2_50_pure[ig] + 1.0
                    for ik in range(0, len(molsOver50)):
                        if mol[ii] == molsOver50[ik]:
                            d2_over[ig] = d2_over[ig] + 1.0
                            if (id[ii] % 40) >= 29:
                                d2_over_pure[ig] = d2_over_pure[ig] + 1.0

# Initial Variables
natoms = 0
num1 = 0
num2 = 0
cirxcm = 0
cirycm = 0
thetacm = 0
zcm = 0
conf = 0
numA = [0] * 480
numUnder50 = 0
num50 = 0
numOver50 = 0

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
        zc1=zeros(dim,float32)
        theta=zeros(dim,float32)
        cirx=zeros(dim,float32)
        ciry=zeros(dim,float32)
        cx=zeros(dim)
        cy=zeros(dim)
        cz=zeros(dim)
        typea=[0]*dim
        mol=[0]*dim
        id=[0]*dim
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
    zbox2 = zbox/2.0
    zbox4 = zbox/4.0
    cirxcm = 0    #initialized to zero
    cirycm = 0
    for j in range(1,dim):
        line = sys.stdin.readline()
        #[ii,molj,typej,x1,x2,x3,n1,n2,n3,v1,v2,v3] = string.split(line)
        
    # den() prefactors
    drg = zbox/mg
    numg = mg + 1
 
print "Reading config file...."

# Read configurations to average out
d1_under = zeros(numg,float32)  
d2_under = zeros(numg,float32)
d1_50 = zeros(numg,float32)  
d2_50 = zeros(numg,float32)
d1_over = zeros(numg,float32)  
d2_over = zeros(numg,float32)
d1_under_pure = zeros(numg,float32)  
d2_under_pure = zeros(numg,float32)
d1_50_pure = zeros(numg,float32)  
d2_50_pure = zeros(numg,float32)
d1_over_pure = zeros(numg,float32)  
d2_over_pure = zeros(numg,float32)
istart = iskip+1
for kconf in range(istart,nconf+iskip+1): 
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
    zbox2 = zbox/2.0
    zbox4 = zbox/4.0
    cirxcm = 0    #initialized to zero
    cirycm = 0
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,x1,x2,x3,n1,n2,n3,v1,v2,v3] = string.split(line)
        k=int(ii)
        id[k] = k
        typea[k] = int(typej)
        mol[k] = int(molj)
        xc[k] = xbox*float(x1)         #scaled coords go from 0 to 1 in dump file
        yc[k] = ybox*float(x2)
        zc[k] = zbox*float(x3)
        if zc[k] <= zbox2:             #combine the two layers of lamellae
            zc1[k] = zc[k]
        elif zc[k] > zbox4 and zc[k] <= 2*zbox4:
            zc1[k] = zc[k] - zbox4
        elif zc[k] > 2*zbox4 and zc[k] <= 3*zbox4:
            zc1[k] = zc[k] - 2*zbox4
        elif zc[k] > 3*zbox4 and zc[k] <= 4*zbox4:
            zc1[k] = zc[k] - 3*zbox4
        theta[k] = (2*pi/zbox4)*zc1[k] #circling z-axis for one layer (combined with others)
        cirx[k] = (zbox4/(2*pi))*sin(theta[k]) #x coordinate in cylindrical coordinates(xytheta), R=zbox/(2*pi)
        ciry[k] = (zbox4/(2*pi))*cos(theta[k]) #y coordinate in cylindrical coordinates(xytheta), R=zbox/(2*pi)        
        if typea[k] == 1:
            cirxcm += cirx[k]/(natoms*0.5) #get the tubular center of mass; mass is not included since it is 1.0
            cirycm += ciry[k]/(natoms*0.5)
            numA[mol[k]-1] += 1.0/nconf  #Count number of A atoms
        thetacm = atan2(-cirxcm,-cirycm) + pi #atan2(X,Y); return theta in positive way from (0,R), arctan(X/Y), range from -pi to pi, atan2(1,1)=pi/4 and atan2(-1,-1)=-3*pi/4; +pi to make the range from 0 to 2pi
        zcm = (zbox4/(2*pi))*thetacm

    #split into groups (<50%, 50%, >50%)
    if kconf == nconf:
        for x in range(0,len(numA)):
            if numA[x] < 20:
                numUnder50 += 1
            elif numA[x] == 20:
                num50 += 1
            else:
                numOver50 += 1   
        molsUnder50 = [0] * numUnder50
        mols50 = [0] * num50
        molsOver50 = [0] * numOver50
        indexUnder = 0
        index50 = 0
        indexOver = 0
        for y in range(0,len(numA)):
            if numA[y] < 20:
                molsUnder50[indexUnder] = y+1 #stores molecule number in list
                indexUnder += 1
            elif numA[y] == 20:
                mols50[index50] = y+1
                index50 += 1
            else:
                molsOver50[indexOver] = y+1
                indexOver += 1  
        print("Under 50%% A: %i (%2.4f%%)" % (numUnder50, (numUnder50*100/480.0)))
        print("50%% A: %i (%2.4f%%)" % (num50, (num50*100/480.0)))
        print("Over 50%% A: %i (%2.4f%%)" % (numOver50, (numOver50*100/480.0)))

    

        # den() prefactors
        drg = zbox/mg
        numg = mg + 1
    
        # call den()
        den()
    
        conf = conf + 1

        dd1_under = d1_under / nconf
        dd2_under = d2_under / nconf
        dd1_50 = d1_50 / nconf
        dd2_50 = d2_50 / nconf
        dd1_over = d1_over / nconf
        dd2_over = d2_over / nconf

        dd1_under_pure = d1_under_pure / nconf
        dd2_under_pure = d2_under_pure / nconf
        dd1_50_pure = d1_50_pure / nconf
        dd2_50_pure = d2_50_pure / nconf
        dd1_over_pure = d1_over_pure / nconf
        dd2_over_pure = d2_over_pure / nconf

        ## OUTPUT
        OUT = open(file, 'w')
        OUT.write("#%7i\n" % (conf))

        OUT.write("Under 50% A\n")
        OUT.write("z den(A) den(B)\n")
        for ig in range(1,numg):
            OUT.write("%8.4f %8.4f %8.4f\n" % (ig,dd1_under[ig],dd2_under[ig]))
        OUT.write("\n")

        OUT.write("50% A\n")
        OUT.write("z den(A) den(B)\n")
        for ig in range(1,numg):
            OUT.write("%8.4f %8.4f %8.4f\n" % (ig,dd1_50[ig],dd2_50[ig]))
        OUT.write("\n")

        OUT.write("Over 50% A\n")
        OUT.write("z den(A) den(B)\n")
        for ig in range(1,numg):
            OUT.write("%8.4f %8.4f %8.4f\n" % (ig,dd1_over[ig],dd2_over[ig]))
        OUT.write("\n")

        OUT.write("Under 50% A (Pure A only)\n")
        OUT.write("z den(A) den(B)\n")
        for ig in range(1,numg):
            OUT.write("%8.4f %8.4f %8.4f\n" % (ig,dd1_under_pure[ig],dd2_under_pure[ig]))
        OUT.write("\n")

        OUT.write("50% A (Pure A only)\n")
        OUT.write("z den(A) den(B)\n")
        for ig in range(1,numg):
            OUT.write("%8.4f %8.4f %8.4f\n" % (ig,dd1_50_pure[ig],dd2_50_pure[ig]))
        OUT.write("\n")

        OUT.write("Over 50% A (Pure A only)\n")
        OUT.write("z den(A) den(B)\n")
        for ig in range(1,numg):
            OUT.write("%8.4f %8.4f %8.4f\n" % (ig,dd1_over_pure[ig],dd2_over_pure[ig]))
        OUT.write("\n")

        OUT.close()




