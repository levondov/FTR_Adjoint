# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 12:27:33 2021

@author: levon

Example script of running the moments with the new MomentSolver Class
"""

from Moments import MomentSolver


# create momentSolver objects
mom = MomentSolver(energy=5e3, current=3.0e-3)


### FTR transformer example ###

# define a quadrupole lattice
dB = [-0.20157, 0.24166, -0.19439] # gradients [T/m]
qstart = [0.00444, 0.12462, 0.2232] # start positions [m]
qend = [0.01444, 0.13462, 0.2332] # end positions [m]
quadrot = [0.883077, 0.8406295, 0.8419267] # rotation angle [rad]

# create lattice
mom.CreateLatticeProfile(dB, qstart, qend, quadrot, zstart=0.0, zend=0.275, repeat=1, verbose=True)

# run moment solver and grab output
z,y,motion,k = mom.RunMoments(verbose=True)

# plot beam size moments
mom.PlotBeamSize()


### Periodic 

# let's make the previous lattice periodicity equal to 3 and repeat the integration

# create lattice , note the change to the repeat variable since we want a 3x periodicity
mom.CreateLatticeProfile(dB, qstart, qend, quadrot, zstart=0.0, zend=0.275, repeat=3, verbose=True)

# run moment solver and grab output
z,y,motion,k = mom.RunMoments(verbose=True)

# plot beam size moments
mom.PlotBeamSize()






