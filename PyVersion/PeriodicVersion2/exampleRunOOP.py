# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 12:27:33 2021

@author: levon

Example script of running the moments with the new MomentSolver Class
"""

from Moments import MomentSolver


# create momentSolver objects
mom = MomentSolver(energy=5e3, current=3.0e-3)

#mom.rigidity = 10


if 0:
    ### FTR transformer example ###
    
    # define a quadrupole lattice
    dB = [-0.20157, 0.24166, -0.19439] # gradients [T/m]
    qstart = [0.00444, 0.12462, 0.2232] # start positions [m]
    qend = [0.01444, 0.13462, 0.2332] # end positions [m]
    quadrot = [0.883077, 0.8406295, 0.8419267] # rotation angle [rad]
    
    
    
    # create lattice
    mom.CreateLatticeProfile(dB, qstart, qend, quadrot, zstart=0.0, zend=0.322, repeat=1, verbose=False)
    
    
    # run moment solver and grab output
    z,y,motion,k = mom.RunMoments(verbose=True)
    z,y,ym,mot,k = mom.RunMomentsAdjoint(verbose=True)
    # plot beam size moments
    mom.PlotBeamSize()
    

if 0:    
    ### Periodic 
    
    # let's make the previous lattice periodicity equal to 3 and repeat the integration
    
    # create lattice , note the change to the repeat variable since we want a 3x periodicity
    mom.CreateLatticeProfile(dB, qstart, qend, quadrot, zstart=0.0, zend=0.275, repeat=3, verbose=True)
    # run moment solver and grab output
    z,y,motion,k = mom.RunMoments(verbose=True)
    # plot beam size moments
    mom.PlotBeamSize()


if 1:
    ### Periodic 2
    
    # define a quadrupole lattice
    qlength = 0.02
    dB = [0.15, -0.15, 0.15] # gradients [T/m]
    qstart = [0.0, 0.03, 0.07] # start positions [m]
    qend = [qstart[0]+qlength/2.0, qstart[1]+qlength, qstart[2]+qlength/2.0] # end positions [m]
    quadrot = [0.0, 0.0, 0.0] # rotation angle [rad]
    
    # create lattice
    mom.CreateLatticeProfile(dB, qstart, qend, quadrot, zstart=0.0, zend=0.08, repeat=10, verbose=True)
    # adjust beam properties
    mom.energy = 10e3
    mom.current = 0.0
    mom.initialConditions = [
       0.005305303507350*1e-4,
      -0.000704766300903*1e-4,
      -0.000003802517788*1e-4,
      -0.000013223559645*1e-4,
       0.000005829207635*1e-4,
      -0.000010686777115*1e-4,
       0.138523562476793*1e-4,
       0.018317745992450*1e-4,
       0.000011963896800*1e-4,
      -0.000003463399474*1e-4,
      0.0
     ]
    mom.h = 100
    mom.CalculateBeamProperties()
    # run moment solver and grab output
    z,y,motion,k = mom.RunMoments(verbose=True)
    # plot beam size moments
    mom.PlotBeamSize(scale=0.2)
    
    import scipy.io as sio
    sio.savemat('run1py.mat',{'data': [mom.z,mom.y]})








