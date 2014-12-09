
#!/usr/bin/env python
import os
from ising import ising
import sys,getopt
opts,args = getopt.getopt(sys.argv[1:],'n:s:i:t:w')

for key, val in opts:

    if key == '-n': n = int(val)
    elif key == '-s': nsteps = int(val)
    elif key == '-t': T = float(val)
    elif key == '-i': index = int(val)
    elif key == '-w': wd = str(val)
lattice, energies, spins = ising(n=n, nsteps=nsteps, T=T)
    
with open(os.path.join(wd,'temp-{1}.out'.format(wd, index)), 'w') as f:
    for i, spin in enumerate(spins):
        if i % 1000 == 0:
            f.write("{0}\t{1}\n".format(i, spin))   
