# Script:  mkinput_gyroid_diblock.py 
# Purpose: Make input file of gyroid phase of diblock copolymers
# Syntax: "python mkinput_gyroid_diblock.py"
# Author: Mitchell Wendt, based on Mark Stevens and YS 10/2014

#  LAMMPS types (all neutral)
#  1         soft polymer beads (A)
#  2         hard polymer beads (B)

import sys, string, math, itertools
from numpy import *
from random import *

def maxg(g1, g2):
    if (math.fabs(g1) > math.fabs(g2)):
        return math.fabs(g1)
    else:
        return math.fabs(g2)

def get_g1(a, k, xc = [], yc = [], zc = [], *args):
    val = math.fabs(10.0 * (cos(xc[k] / a) * sin(yc[k] / a) + cos(yc[k] / a) * sin(zc[k] / a) + cos(zc[k] / a) * sin(xc[k] / a)) - 0.5 * (cos(2 * (xc[k] / a)) * cos(2 * (yc[k] / a)) + cos(2 * (yc[k] / a)) * cos(2 * (zc[k] / a)) + cos(2 * (zc[k] / a)) * cos(2 * (xc[k] / a))))
    return val

def get_g2(a, k, xc = [], yc = [], zc = [], *args):
    val = math.fabs(10.0 * (cos(-1 * xc[k] / a) * sin(-1 * yc[k] / a) + cos(-1 * yc[k] / a) * sin(-1 * zc[k] / a) + cos(-1 * zc[k] / a) * sin(-1 * xc[k] / a)) - 0.5 * (cos(2 * (-1 * xc[k] / a)) * cos(2 * (-1 * yc[k] / a)) + cos(2 * (-1 * yc[k] / a)) * cos(2 * (-1 * zc[k] / a)) + cos(2 * (-1 * zc[k] / a)) * cos(2 * (-1 * xc[k] / a))))
    return val

# INPUT PARAMETERs
nbeads = 40                        # number of beads in monomer
nmonomersperpoly = 1               # total (including x) number of monomers in polymer
npoly = 1000                       # number of polymers
ntot = npoly * nbeads              # number of total beads
ntypes = 2                         # number of monomer types
dens = 0.85                        # bead number density
bond = 0.97                        # bond length
fA= 0.325                          # fraction of A monomers
g0 = 10.23                         # outer bound of generating gyroid phase
diblock_split = int(nbeads * fA)   # molecule at which chains switch from A to B
interval = 0.005                   # interval around g0

# type names: conducting block(PEO), nonconducting block(PS)
#           0     1     2     
atype = ( '00', 'CB', 'NB' )

# Define monomer sequence: diblock copolymers
sequence = [1] * int(nbeads * fA) + [2] * (int(nbeads * (1 - fA)))

# Files
INPUT_LAMMPS = open('input_gyroid.lammps', 'w')

# Simulation cell parameters
dim = ntot + 1
vol = ntot / dens
nbonds = ntot - npoly
nmonomers = nmonomersperpoly * npoly
hx = vol ** (1. / 3.)
hy = hx
hz = hx
hx2 = hx / 2.
hy2 = hy / 2.
hz2 = hz / 2.
nbonds = ntot - npoly

a = hx / (2. * math.pi)

print
print "Total number of beads:", ntot
print "Number of chains =", npoly
print "beads in monomer =", nbeads
print "monomers total =", nmonomers
print "Number of atoms types = ",ntypes
print "dens = ", dens
print "vol = ", vol
print "metric: %10.4f %10.4f %10.4f\n\n" % (hx, hy, hz)

# init position variables
xc = zeros(dim, float32)
yc = zeros(dim, float32)
zc = zeros(dim, float32)
g1 = zeros(dim, float32)
g2 = zeros(dim, float32)
gmax = zeros(dim, float32)
cx = zeros(dim)
cy = zeros(dim)
cz = zeros(dim)


# BUILD POLYMERS
typeb = [0] * dim
molnum = [0] * dim
k = 0

#g1 = 10.0 * (cos(x / a) * sin(y / a) + cos(y / a) * sin(z / a) + cos(z / a) * sin(x / a)) - 0.5 * (cos(2 * (x / a)) * cos(2 * (y / a)) + cos(2 * (y / a)) * cos(2 * (z / a)) + cos(2 * (z / a)) * cos(2 * (x / a))) 

#g2 = 10.0 * (cos(-x / a) * sin(-y / a) + cos(-y / a) * sin(-z / a) + cos(-z / a) * sin(-x / a)) - 0.5 * (cos(2 * (-x / a)) * cos(2 * (-y / a)) + cos(2 * (-y / a)) * cos(2 * (-z / a)) + cos(2 * (-z / a)) * cos(2 * (-x / a)))

k = 0

for ix in xrange(npoly):
    lengthcurrentpoly = 0
    for iy in range(nmonomersperpoly):
        currentmonomer = ix * nmonomersperpoly + iy
        seqnum = 0
        seq = sequence
        for iz in seq:
                seqnum = seqnum + 1
                k = k + 1
                kk = k + diblock_split - 1 - ((k - 1) % nbeads) * 2
                lengthcurrentpoly = lengthcurrentpoly + 1
                if iy == 0 and seqnum == 1:
                    typeb[kk] = iz
                    molnum[kk] = ix + 1
                    xc[kk] = random() * hx
                    yc[kk] = random() * hy
                    zc[kk] = random() * hz
                    g1[kk] = get_g1(a, kk, xc, yc, zc)
                    g2[kk] = get_g2(a, kk, xc, yc, zc)
                    gmax[kk] = maxg(g1[kk], g2[kk])
                    print k, " ", gmax[k], " ", g1[k], " ", g2[k], " ", kk
                    while (gmax[kk] > (g0 + interval)) or (gmax[kk] < (g0 - interval)):
                        xc[kk] = random() * hx
                        yc[kk] = random() * hy
                        zc[kk] = random() * hz
                        g1[kk] = get_g1(a, kk, xc, yc, zc)
                        g2[kk] = get_g2(a, kk, xc, yc, zc)
                        gmax[kk] = maxg(g1[kk], g2[kk])
                        print k, " ", gmax[kk], " ", g1[kk], " ", g2[kk], " ", kk
                elif seqnum <= diblock_split:
                    typeb[kk] = iz
                    molnum[kk] = ix + 1
                    theta = random() * 2 * pi
                    dz = random() * 2 - 1
                    dx = sqrt(1 - dz ** 2) * cos(theta)
                    dy = sqrt(1 - dz ** 2) * sin(theta)
                    r = sqrt(dx*dx + dy*dy + dz*dz)
                    scale = bond / r
                    dx = scale * dx
                    dy = scale * dy
                    dz = scale * dz
                    xc[kk] = xc[kk + 1] + dx
                    yc[kk] = yc[kk + 1] + dy
                    zc[kk] = zc[kk + 1] + dz
                    g1[kk] = get_g1(a, kk, xc, yc, zc)
                    g2[kk] = get_g2(a, kk, xc, yc, zc)
                    gmax[kk] = maxg(g1[kk], g2[kk])
                    print k, " ", gmax[kk], " ", g1[kk], " ", g2[kk], " ", kk
                    while gmax[kk] < g0:
                        theta = random() * 2 * pi
                        dz = random() * 2 - 1
                        dx = sqrt(1 - dz ** 2) * cos(theta)
                        dy = sqrt(1 - dz ** 2) * sin(theta)
                        r = sqrt(dx*dx + dy*dy + dz*dz)
                        scale = bond / r
                        dx = scale * dx
                        dy = scale * dy
                        dz = scale * dz
                        xc[kk] = xc[kk + 1] + dx
                        yc[kk] = yc[kk + 1] + dy
                        zc[kk] = zc[kk + 1] + dz
                        g1[kk] = get_g1(a, kk, xc, yc, zc)
                        g2[kk] = get_g2(a, kk, xc, yc, zc)
                        gmax[kk] = maxg(g1[kk], g2[kk])
                        print k, " ", gmax[kk], " ", g1[kk], " ", g2[kk], " ", kk
                else:
                    typeb[k] = iz
                    molnum[k] = ix + 1
                    theta = random() * 2 * pi
                    dz = random() * 2 - 1
                    dx = sqrt(1 - dz ** 2) * cos(theta)
                    dy = sqrt(1 - dz ** 2) * sin(theta)
                    r = sqrt(dx*dx + dy*dy + dz*dz)
                    scale = bond / r
                    dx = scale * dx
                    dy = scale * dy
                    dz = scale * dz
                    xc[k] = xc[k - 1] + dx
                    yc[k] = yc[k - 1] + dy
                    zc[k] = zc[k - 1] + dz
                    g1[k] = get_g1(a, k, xc, yc, zc)
                    g2[k] = get_g2(a, k, xc, yc, zc)
                    gmax[k] = maxg(g1[k], g2[k])
                    print k, " ", gmax[k], " ", g1[k], " ", g2[k], " ", kk
                    while gmax[k] > g0:
                        theta = random() * 2 * pi
                        dz = random() * 2 - 1
                        dx = sqrt(1 - dz ** 2) * cos(theta)
                        dy = sqrt(1 - dz ** 2) * sin(theta)
                        r = sqrt(dx*dx + dy*dy + dz*dz)
                        scale = bond / r
                        dx = scale * dx
                        dy = scale * dy
                        dz = scale * dz
                        xc[k] = xc[k - 1] + dx
                        yc[k] = yc[k - 1] + dy
                        zc[k] = zc[k - 1] + dz
                        g1[k] = get_g1(a, k, xc, yc, zc)
                        g2[k] = get_g2(a, k, xc, yc, zc)
                        gmax[k] = maxg(g1[k], g2[k])
                        print k, " ", gmax[k], " ", g1[k], " ", g2[k], " ", kk

print "Gyroid Box Generated"

# IMAGE FLAGS
for k in xrange(1, ntot + 1):
    if (xc[k] > hx):
        cx[k] = int(xc[k] / hx)
        xc[k] = xc[k] - cx[k] * hx - hx2
    elif (xc[k] < 0.0):
        cx[k] = -int((-xc[k] + hx) / hx)
        xc[k] = xc[k] - cx[k] * hx - hx2
    else:
        cx[k] = 0
        xc[k] = xc[k] - hx2
        
    if (yc[k] > hy):
        cy[k] = int(yc[k] / hy)
        yc[k] = yc[k] - cy[k] * hy - hy2
    elif (yc[k] < 0.0):
        cy[k] = -int((-yc[k] + hy) / hy)
        yc[k] = yc[k] - cy[k] * hy - hy2
    else:
        cy[k] = 0
        yc[k] = yc[k] - hy2
        
    if (zc[k] > hz):
        cz[k] = int(zc[k] / hz)
        zc[k] = zc[k] - cz[k] * hz - hz2
    elif (zc[k] < 0.0):
        cz[k] = -int((-zc[k] + hz)/hz)
        zc[k] = zc[k] - cz[k] * hz - hz2
    else:
        cz[k] = 0
        zc[k] = zc[k] - hz2

print "Polymers built."


# OUTPUT headers ---------------------------------------------------------------

# input.lammps header 
INPUT_LAMMPS.write("#Diblock Copolymer KS 1/2016\n")
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atoms\n" %     ntot)
INPUT_LAMMPS.write("%10i    bonds\n" %     nbonds)
INPUT_LAMMPS.write("%10i    angles\n" %    0)
INPUT_LAMMPS.write("%10i    dihedrals\n" % 0)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atom types\n" % ntypes)
INPUT_LAMMPS.write("%10i    bond types\n" % 2)
INPUT_LAMMPS.write("%10i    angle types\n" % 0)
INPUT_LAMMPS.write("%10i    dihedral types\n" % 0)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write(" %16.8f %16.8f   xlo xhi\n" % (-hx2,hx2))
INPUT_LAMMPS.write(" %16.8f %16.8f   ylo yhi\n" % (-hy2,hy2))
INPUT_LAMMPS.write(" %16.8f %16.8f   zlo zhi\n" % (-hz2,hz2))
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Atoms\n")
INPUT_LAMMPS.write("\n")

# END OUTPUT headers -----------------------------------------------------------
# Atoms output
mass = 1.0

# Polymers
i = 0
imol = 0

# positions 
for i in xrange(1, dim):
    itype = typeb[i]
    aname = atype[itype]
    if itype != 3:
        imol = molnum[i]
    elif itype == 3:
        imol = i - ntot + npoly
    INPUT_LAMMPS.write("%6i %6i %2i %9.4f %9.4f %9.4f %6i %6i %6i\n" % (i, imol, typeb[i], xc[i], yc[i], zc[i], cx[i], cy[i], cz[i]))


# Bonds
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Bonds\n")
INPUT_LAMMPS.write("\n")
ibond = 0

for i in xrange(1, ntot):
        #if not at the end of the polymer
        if molnum[i + 1] == molnum[i]:
            ibond = ibond + 1 #the bond number
            j = i + 1
            if i % nbeads == int(nbeads * fA):
                INPUT_LAMMPS.write("%8i  2 %8i %8i\n" % (ibond, i, j))
            else: 
                INPUT_LAMMPS.write("%8i  1 %8i %8i\n" % (ibond, i, j))


# Masses
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Masses\n")
INPUT_LAMMPS.write("\n")
for ii in xrange(1, ntypes + 1):
    INPUT_LAMMPS.write("%3i  1.0\n" % ii)


#Close files
INPUT_LAMMPS.close()

print "LAMMPS output complete."
