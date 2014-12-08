#!/usr/bin/env python
"""\
 ising.py - Play with the Ising model
"""

from math import exp
from random import randrange,choice,random
from numpy import zeros
import numpy as np

def init_ising_lattice(n):
    lattice = zeros((n,n),dtype=int)
    options = [-1,1]
    for i in range(n):
        for j in range(n):
            lattice[i,j] = choice(options)
    return lattice

def energydiff(S0,Sn,J,H): return 2*S0*(H+J*Sn)

def ising(n=200,nsteps=500000,H=0,J=1,T=1):
    lattice = init_ising_lattice(n)
    energy = 0
    energies = []
    lattices = []
    for step in range(nsteps):
        i = randrange(n)
        j = randrange(n)
        Sn = lattice[(i-1)%n,j]+lattice[(i+1)%n,j]+\
             lattice[i,(j-1)%n]+lattice[i,(j+1)%n]
        dE = energydiff(lattice[i,j],Sn,J,H)
        if dE < 0 or random() < exp(-dE/T):
            lattice[i,j] = -lattice[i,j]
            energy += dE
            lattices.append(lattice)
            energies.append(energy)
    return lattices,energies

def pil_image(lattice,fname="ising.png"):
    # creates a snapshot image of the current state of the lattice
    import Image, ImageDraw
    n,m = lattice.shape
    img = Image.new("RGB",(m,n),(255,255,255))
    draw = ImageDraw.Draw(img)

    for i in range(n):
        for j in range(m):
            if lattice[i,j] > 0: draw.point((i,j),(0,0,0))
    img.save(fname,"PNG")
    return
    

def main():
    import sys,getopt
    import matplotlib.pyplot as plt
    opts,args = getopt.getopt(sys.argv[1:],'n:s:h:j:t:')
    n = 200
    nsteps = 500000
    H = 0
    J = 1
    T = 20
    for key,val in opts:
        if key == '-n': n = int(val)
        elif key == '-s': nsteps = int(val)
        elif key == '-h': H = float(val)
        elif key == '-j': J = float(val)
	elif key == '-t': T = float(val)
    lattices,energies = ising(n,nsteps,H,J,T)
    S = []
    for lattice in lattices:
        S.append(np.sum(lattice)/float(n**2))
    plt.plot(range(len(S)), S)
    plt.show()
    #pil_image(lattice)
    #print len(energies)/float(nsteps)
    #plot(energies)
    return

def plot(energies):
    # if you have Gnuplot module, you can get the plot automatically:
    #
    # import Gnuplot
    # p = Gnuplot.Gnuplot() 
    # d = Gnuplot.Data(range(len(energies)),energies)
    # p.xlabel('MC cycle')
    # p.ylabel('Energy')
 
    # p.plot(d)
    # ans = raw_input('Enter f to create .png file, or Enter to quit ')
    # if ans == 'f':
    #   p.hardcopy('energies.png',terminal = 'png')
    # p.reset()

    # but we'll just write the energies out to a file:

    outFile = open('ising.out', 'w') 

    for i in range(len(energies)):
       outFile.write("%d\t%f\n" % ( i, energies[i]))

    outFile.close()

    return

if __name__ == '__main__': main()
