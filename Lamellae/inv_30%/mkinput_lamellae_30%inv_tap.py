# Script:  mkinput_lamellae_30%inv_tap.py
# Purpose: Make input file of few layers lamellae of tapered diblock copolymers in Kremer-Grest approach (JCP, 105, 1996); graft half-polymers on two facing walls and then mirror them to the opposite directions to complete one layer)
# Syntax: "python mkinput_lamellae_30%inv_tap.py"
# Author: Modified Mark Stevens's mkinput.py by Youngmi Seo (10/2014)

# chains are random walk.

#  LAMMPS types (all neutral)
#  1         soft polymer beads (A)
#  2         hard polymer beads (B)

import sys, string, math, itertools
from numpy import *
from random import *
#-------------------------------------------------------------------------

## 40 bead/monomes, 480 poly = 19200 beads
## ADJUSTED DENSITY
## INPUT PARAMETERs

#iseed = 9328             # random number seed
#seed(iseed)              # to generate same set of random numbers every time
nbeads = 40               # number of beads in monomer
hnbeads = nbeads/2        # half number of beads in monomer
nmonomersperpoly = 1      # total (including x) number of monomers in polymer
npoly = 300               # number of polymers
hnpoly = npoly/2          # half number of polymers
ntot = npoly*nbeads       # number of total beads
nlayers = 3               # number of lamellar layers
ntypes = 2                # number of monomer types
nbtypes = 2               # number of bond types
minsep = 1.0              # allowed separation in overlap check
dens = 0.85               # bead number density
cdens = 0.1               # coverage density on a planar surface
bond = 0.97               # bond length. depends on bond potential, but close to 1 is good enough
percentTaper = 30./100.                       # tapered midblock length in percent
npurebeads = int(nbeads*(1-percentTaper))     # number of pure beads on two ends
ntapers = nbeads - npurebeads                 # number of beads on taper midblock

# Label schemes
# type names: conducting block(PEO), nonconducting block(PS)
#                       1                      2   

atype = ( '00', 'A' , 'B' )

# Define monomer sequence: inverse tapered diblock copolymers
pureA = [1]*(npurebeads/2)   #modify these lines to control fraction of monomers
pureB = [2]*(npurebeads/2)
taperAB = [0]*ntapers
j = 0
for ij in range(ntapers):
    j = j + 1
    taperAB[j-1] = [1]*int(round(npoly*(j-0.5)/ntapers)) + [2]*(npoly-int(round(npoly*(j-0.5)/ntapers)))
    shuffle(taperAB[j-1])
monomerTypeTaper = [[row[i] for row in taperAB] for i in range(npoly)]
    


# END INPUT Parameters ------------------------------------------------------
# files
INPUT_LAMMPS = open('input_30inv_tap.lammps', 'w')

# Simulation cell parameters for grafting polymers
gnpoly = npoly/nlayers   #for grafting polymers
gntot = gnpoly*hnbeads
gvol = gntot/dens
surf = gnpoly/2/cdens
ghx = surf**0.5
ghy = surf**0.5
ghz = gvol/surf
ghx2 = ghx/2.
ghy2 = ghy/2.
ghz2 = ghz/2.
gvol = ghx*ghy*ghz
gnbonds = gntot-gnpoly

# Simulation cell parameters
dim = ntot
vol = ntot/dens
hx = ghx
hy = ghy
hz = ghz*2*nlayers + (nlayers*2 - 1)
hx2 = hx/2.
hy2 = hy/2.
hz2 = hz/2.
vol = hx*hy*hz
nbonds = ntot-npoly

print "Total number of beads:",ntot
print "Number of chains =", npoly
print "Number of Beads in monomer =", nbeads
print "Taper percent =", percentTaper*100,"%"
print "Number of atoms types = ",ntypes
#print "seed = ", iseed
print " "
print "Geometry:"
print "dens = ", dens
print "vol = ", vol
print " "
print "metric: %10.4f %10.4f %10.4f\n\n" % (hx, hy, hz)

# Initial variables
xc=zeros(dim+1,float32)
yc=zeros(dim+1,float32)
zc=zeros(dim+1,float32)
cx=zeros(dim+1)
cy=zeros(dim+1)
cz=zeros(dim+1)

rg2ave=0.0
rgave=0.0
rend2ave = 0.0
typeb=[0]*(dim+1)
molnum=[0]*(dim+1)
k=0

## BUILD POLYMERS
for nl in range(0,nlayers):
    for ix in range((hnpoly/nlayers)*nl*2,(hnpoly/nlayers)*(nl*2+1)):      
        lengthcurrentpoly = 0
        sequence = pureA + monomerTypeTaper[ix] + pureB
        for iy in range(nmonomersperpoly):
            currentmonomer = ix*nmonomersperpoly + iy
            seqnum = 0
            seq = sequence
            for iz in seq:
                seqnum = seqnum + 1
                k = k + 1
                kk = k + 19 - ((k-1) % 20)*2
                lengthcurrentpoly = lengthcurrentpoly + 1
                if iy == 0 and seqnum == 1:
                    typeb[k] = iz
                    molnum[kk] = ix + 1 
                    k1 = k
                    xc[kk] = random()*ghx
                    yc[kk] = random()*ghy
                    zc[kk] = 0.0000 + (ghz*nl*2 + nl*2)
                elif seqnum <= hnbeads:          # pick random direction; scale to be bond length
                    typeb[k] = iz
                    molnum[kk] = ix + 1 
                    theta = random()*2*pi
                    dz = random()*2 - 1 
                    dx = sqrt(1-dz**2)*cos(theta)
                    dy = sqrt(1-dz**2)*sin(theta)                  
                    r = sqrt(dx*dx+dy*dy+dz*dz)
                    scale = bond/r
                    dx = scale*dx
                    dy = scale*dy
                    dz = scale*dz
                    xc[kk] = xc[kk+1] + dx
                    yc[kk] = yc[kk+1] + dy
                    zc[kk] = zc[kk+1] + dz        # don't need to shift, already done for the first bead
                    if zc[kk] >= (ghz + (ghz*nl*2 + nl*2)) or zc[kk] <= (ghz*nl*2 + nl*2):
                        zc[kk] = zc[kk+1] - dz
                else:
                    typeb[k] = iz
                    molnum[k] = ix + 1
                    xc[k] = xc[kk-20]
                    yc[k] = yc[kk-20]
                    zc[k] = zc[kk-20] - (zc[kk-20] - zc[kk-20 + ((k-21) % 20)] + 0.5)*2
                        
    for ix in range((hnpoly/nlayers)*(nl*2+1),(hnpoly/nlayers)*(nl*2+2)):
        sequence = pureA + monomerTypeTaper[ix] + pureB
        for iy in range(nmonomersperpoly):
            currentmonomer = ix*nmonomersperpoly + iy
            seqnum = 0
            seq = sequence
            for iz in seq:
                seqnum = seqnum + 1
                k = k + 1
                kk = k + 19 - ((k-1) % 20)*2
                lengthcurrentpoly = lengthcurrentpoly + 1                
                if iy == 0 and seqnum == 1:
                    typeb[k] = iz
                    molnum[kk] = ix + 1 
                    k1 = k
                    xc[kk] = random()*ghx
                    yc[kk] = random()*ghy
                    zc[kk] = ghz + (ghz*nl*2 + nl*2)
                elif seqnum <= hnbeads:        # pick random direction; scale to be bond length
                    typeb[k] = iz
                    molnum[kk] = ix + 1 
                    theta = random()*2*pi
                    dz = random()*2 - 1
                    dx = sqrt(1-dz**2)*cos(theta)
                    dy = sqrt(1-dz**2)*sin(theta)
                    r = sqrt(dx*dx+dy*dy+dz*dz)
                    scale = bond/r
                    dx = scale*dx
                    dy = scale*dy
                    dz = scale*dz   
                    xc[kk] = xc[kk+1] + dx
                    yc[kk] = yc[kk+1] + dy
                    zc[kk] = zc[kk+1] + dz
                    if zc[kk] >= (ghz + (ghz*nl*2 + nl*2)) or zc[kk] <= (ghz*nl*2 + nl*2):
                        zc[kk] = zc[kk+1] - dz
                else:
                    typeb[k] = iz
                    molnum[k] = ix + 1
                    xc[k] = xc[kk-20]
                    yc[k] = yc[kk-20]
                    zc[k] = zc[kk-20] - (zc[kk-20] - zc[kk-20 + ((k-21) % 20)] - 0.5)*2

# Calculate R and R_G
#    k2= k1 + lengthcurrentpoly
#    xcm = sum(xc[k1:k2+1])/lengthcurrentpoly
#    ycm = sum(yc[k1:k2+1])/lengthcurrentpoly
#    zcm = sum(zc[k1:k2+1])/lengthcurrentpoly
#    xg = xc[k1:k2+1]-xcm
#    yg = yc[k1:k2+1]-ycm
#    zg = zc[k1:k2+1]-zcm
#    rg2 = (dot(xg,xg) + dot(yg,yg) + dot(zg,zg))/lengthcurrentpoly  # radius of gyration
#    rend2 = (xc[k1]-xc[k2])**2 + (yc[k1]-yc[k2])**2 + (zc[k1]-zc[k2])**2  # end to end
#    rend2ave = rend2ave + rend2
#    rg2ave = rg2ave + rg2
#    rgave = rgave + sqrt(rg2)
#    print "current rg", rg2

## IMAGE FLAGS
for k in xrange(1,ntot+1):
    if (xc[k] > hx):
        cx[k] = int(xc[k]/hx)
        xc[k] = xc[k] - cx[k]*hx - hx2
    elif (xc[k] < 0.0):
        cx[k] = -int((-xc[k]+hx)/hx)
        xc[k] = xc[k] - cx[k]*hx - hx2
    else:
        cx[k] = 0
        xc[k] = xc[k] - hx2
        
    if (yc[k] > hy):
        cy[k] = int(yc[k]/hy)
        yc[k] = yc[k] - cy[k]*hy - hy2
    elif (yc[k] < 0.0):
        cy[k] = -int((-yc[k]+hy)/hy)
        yc[k] = yc[k] - cy[k]*hy - hy2
    else:
        cy[k] = 0
        yc[k] = yc[k] - hy2
        
    if (zc[k] > hz):
        cz[k] = int(zc[k]/hz)
        zc[k] = zc[k] - cz[k]*hz - hz2
    elif (zc[k] < 0.0):
        cz[k] = -int((-zc[k]+hz)/hz)
        zc[k] = zc[k] - cz[k]*hz - hz2
    else:
        cz[k] = 0
        zc[k] = zc[k] - hz2

print "Polymers built."

#rg2ave = rg2ave/npoly
#rgave = rgave/npoly
#rend2ave = rend2ave/npoly
#rave = rend2ave/rg2ave
#print "<R_G^2> <R_G> = ",rg2ave,rgave
#print "<R_end^2>= ",rend2ave


## OUTPUT Headers ---------------------------------------------------------------

# Input.lammps header 
INPUT_LAMMPS.write("# Kremer-Grest type Lamellae YS 10/2014\n")
INPUT_LAMMPS.write("# Number of lamellar layers = %1.0f\n" % nlayers)
INPUT_LAMMPS.write("# Tapered midblock Lengths = %1.1f\n" % percentTaper)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atoms\n" % dim)
INPUT_LAMMPS.write("%10i    bonds\n" % nbonds)
INPUT_LAMMPS.write("%10i    angles\n" % 0)
INPUT_LAMMPS.write("%10i    dihedrals\n" % 0)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atom types\n" % ntypes)
INPUT_LAMMPS.write("%10i    bond types\n" % nbtypes)
INPUT_LAMMPS.write("%10i    angle types\n" % 0)
INPUT_LAMMPS.write("%10i    dihedral types\n" % 0)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write(" %16.8f %16.8f   xlo xhi\n" % (-hx2,hx2))
INPUT_LAMMPS.write(" %16.8f %16.8f   ylo yhi\n" % (-hy2,hy2))
INPUT_LAMMPS.write(" %16.8f %16.8f   zlo zhi\n" % (-hz2,hz2))
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Atoms\n")
INPUT_LAMMPS.write("\n")

## END OUTPUT Headers -----------------------------------------------------------
mass = 1.0
i = 0
imol = 0

# Positions 
for i in xrange(1,dim+1):
    itype = typeb[i]
    #aname = atype[itype]
    if itype != 3:
        imol = molnum[i]
    elif itype == 3:       #for the case there are penetrants or ions added
        imol = i-ntot+npoly
    INPUT_LAMMPS.write("%6i %6i %2i %9.4f %9.4f %9.4f %6i %6i %6i\n" % (i, imol, typeb[i], xc[i], yc[i], zc[i], cx[i], cy[i], cz[i]))

# Bonds
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Bonds\n")
INPUT_LAMMPS.write("\n")
jbond1 = zeros(nbonds+1)
jbond2 = zeros(nbonds+1)
pbond1 = zeros(nbonds+1)
pbond2 = zeros(nbonds+1)
ibond=0
for i in xrange(1,ntot):        
        if molnum[i+1] == molnum[i]:       #if not at the end of the polymer
            ibond = ibond+1 #bond number
            j=i+1
            if typeb[i] != typeb[j]:
                INPUT_LAMMPS.write("%8i  2 %8i %8i\n" % (ibond,i,j))
            else: 
                INPUT_LAMMPS.write("%8i  1 %8i %8i\n" % (ibond,i,j))


# Masses
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Masses\n")
INPUT_LAMMPS.write("\n")

for ii in xrange(1,ntypes+1):
    INPUT_LAMMPS.write("%3i  1.0\n" % ii)

# Close files
INPUT_LAMMPS.close()

print "LAMMPS output complete."



