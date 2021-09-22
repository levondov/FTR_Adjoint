# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 12:27:33 2021

@author: levon

Example script of running the moments with the new MomentSolver Class
"""

from Moments import MomentSolver


# create momentSolver objects
mom = MomentSolver(energy=10e3, current=0.0)

# change some properties
mom.initialConditions = [1e-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
mom.h = 0.0001

# define a quadrupole lattice
dB = [0.02661, -0.02661, 0.02661] # gradients [T/m]
qstart = [0.0, 0.052125, 0.131375] # start positions [m]
qend = [0.027125, 0.106375, 0.1585] # end positions [m]
quadrot = [0.0, 0.0, 0.0] # rotation angle [rad]

# create lattice
mom.CreateLatticeProfile(dB, qstart, qend, quadrot, zstart=0.0, zend=0.1585, repeat=5, verbose=False)

# run moment solver and grab output
z,y,motion,k = mom.RunMoments(verbose=True)

# plot beam size moments
mom.PlotBeamSize()
