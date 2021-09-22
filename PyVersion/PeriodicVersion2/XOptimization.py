# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 20:54:20 2021

@author: levon
"""

import numpy as np

from Moments import MomentSolver

def a2params(a):
    #dB = [0.02661, -0.02661, 0.02661] # gradients [T/m]
    dB = [-0.20157, 0.24166, -0.19439]
    return [dB[0]*a[0], dB[1]*a[1], dB[2]*a[2]]

def UpdateLattice(mom,a):
    # define a quadrupole lattice
    dB = a2params(a) # gradients [T/m]
    #qstart = [0.0, 0.052125, 0.131375] # start positions [m]
    #qend = [0.027125, 0.106375, 0.1585] # end positions [m]
    #quadrot = [0.0, 0.0, 0.0] # rotation angle [rad]
    qstart = [0.00444, 0.12462, 0.2232] # start positions [m]
    qend = [0.01444, 0.13462, 0.2332] # end positions [m]
    quadrot = [0.883077, 0.8406295, 0.8419267] # rotation angle [rad]
    
    # create lattice
    #mom.CreateLatticeProfile(dB, qstart, qend, quadrot, zstart=0.0, zend=0.1585, repeat=1, verbose=False)
    mom.CreateLatticeProfile(dB, qstart, qend, quadrot, zstart=0.0, zend=0.25, repeat=1, verbose=False)
    
    return mom

def GetdF(mom, a):

    # grab base O and N matrices
    OmatBase = np.copy(mom.Omat)
    NmatBase = np.copy(mom.Nmat)
    
    # copy the params and dF creation
    aCopy = np.copy(a)
    dF = np.zeros((1,len(a)))

    perturb = 0.01    
    for i,ai in enumerate(aCopy):
        
        # perturb a parameter
        a[i] = aCopy[i] + aCopy[i]*perturb
        
        # reGenerate lattice
        mom = UpdateLattice(mom, a)
        
        # reCalculate new O,N matrix
        Onew, Nnew = mom.RecalculateON()

        # calculate perturbation in O and N matrix
        for j in range(len(mom.z)):
            Onew[0,j] = Onew[0,j] - OmatBase[0,j]
            Nnew[0,j] = Nnew[0,j] - NmatBase[0,j]
    
        # calculate the adjoint integral via O and N matrix perturbation
        dF[0,i] = mom.GetAdjointIntegral(Onew, Nnew)
    
    return dF

#############################
# Starting setup
#############################

# tuning parameters
an = np.array([1.2, 0.9, 1.1])

# initial moments
X0 = [1e-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] # initial conditions for the beam

# create momentSolver objects
mom = MomentSolver(energy=5e3, current=3.0e-3)#, initialMoments=X0)

# create lattice
mom = UpdateLattice(mom, an)

#############################
# Optimization
#############################

### compute adjoint equations
mom.RunMoments(verbose=True)
mom.RunMomentsAdjoint(verbose=True)

### compute initial FoM values
f0,fp0,_,_ = mom.GetFoMAndDFoM()
df0 = GetdF(mom, an)
gamma0 = f0 / np.sum(df0**2)

### Define data arrays to store everthing
gammah = np.array([gamma0]).reshape(1,1)
anh = np.copy(an).reshape(1,len(an))
fh = np.array([f0]).reshape(1,1)
fph = np.copy(fp0).reshape(1,len(fp0))
dfh = np.copy(df0)

### Solve for the starting value of gamma
# iterate gradient descent
tmp = (anh[0,:] - gammah[0,0] * dfh[0,:]).reshape(1,len(an))
anh = np.concatenate((anh, tmp))
# update the lattice based on the new parameters 
mom = UpdateLattice(mom,anh[-1,:])
mom.RunMoments()
# grab the value of the FoM
tmp1,tmp2,_,_ = mom.GetFoMAndDFoM()
fh = np.concatenate((fh,tmp1.reshape(1,1)))
fph = np.concatenate((fph, tmp2.reshape(1,len(tmp2))))
print('FoM: ',fh[-1,0])
# Now repeat while decreasing gamma until we find a value that actually reduces the FoM
jj = 1
while ( fh[-1,0] >= f0 ):
    gammatmp = gammah[-1,0] / 2.0
    gammah = np.concatenate((gammah, gammatmp.reshape(1,1)))
    
    tmp = (anh[0,:] - gammah[-1,0] * dfh[0,:]).reshape(1,len(an))
    anh = np.concatenate((anh, tmp))
    
    mom = UpdateLattice(mom,anh[-1,:])
    mom.RunMoments()
    tmp1,tmp2,_,_ = mom.GetFoMAndDFoM()
    fh = np.concatenate((fh,tmp1.reshape(1,1)))
    fph = np.concatenate((fph, tmp2.reshape(1,len(tmp2))))
    print('FoM: ',fh[-1,0])
    
    jj +=1
    if jj > 20: break

### Start the optimization
while True:
    ii = 1 # counter to keep track of how many successful gradient descent steps we have taken
    
    while fh[-1,0] < fh[-2,0]:
        print('Iterating: ',str(ii))

        # iterate
        tmp = (anh[-1,:] - gammah[-1,0] * dfh[-1,:]).reshape(1,len(an))
        anh = np.concatenate((anh, tmp))
        
        # Update lattice and compute FoM
        mom = UpdateLattice(mom,anh[-1,:])
        mom.RunMoments()
        tmp1,tmp2,_,_ = mom.GetFoMAndDFoM()
        fh = np.concatenate((fh,tmp1.reshape(1,1)))
        fph = np.concatenate((fph, tmp2.reshape(1,len(tmp2))))    
        print('FoM: ',fh[-1,0])
        
        # update counter
        ii += 1
        
        # if we have taken too many (small) steps, increase gamma size
        if ( ii > 20 ):
            gammatmp = gammah[-1,0] * 2.0
            gammah = np.concatenate((gammah, gammatmp.reshape(1,1)))
            
    # if the gradient descent stop reducing the FoM, recalculate the adjoint equations
    print('Recomputing adjoint equations')
    
    # grab the last good setting
    tmp = anh[-2,:]
    tmp1 = fh[-2,0]
    tmp2 = fph[-2,:]
    anh = np.concatenate((anh, tmp.reshape(1,len(an))))
    fh = np.concatenate((fh,tmp1.reshape(1,1)))
    fph = np.concatenate((fph, tmp2.reshape(1,len(tmp2))))
    
    # calculate adjoint equations
    mom = UpdateLattice(mom,anh[-1,:])
    mom.RunMoments()
    mom.RunMomentsAdjoint()
    
    # calculate the new gradient
    tmp = GetdF(mom, anh[-1,:])
    dfh = np.concatenate((dfh, tmp))
    
    # update gamma if we recalculated adjoint and there was no improvement
    iii = 1
    if (ii == 2):
        print('Updating gamma')
        f0n = fh[-1,0]
        ann = anh[-1,:]
        
        while fh[-1,0] >= f0n:
            if iii > 20: break
            # update gamma
            gammatmp = gammah[-1,0] / 2.0
            gammah = np.concatenate((gammah, gammatmp.reshape(1,1)))
            
            # iterate
            tmp = (ann - gammah[-1,0] * dfh[-1,:]).reshape(1,len(an))
            anh = np.concatenate((anh, tmp))
            
            # grab FoM
            mom = UpdateLattice(mom,anh[-1,:])
            mom.RunMoments()
            tmp1,tmp2,_,_ = mom.GetFoMAndDFoM()
            fh = np.concatenate((fh,tmp1.reshape(1,1)))
            fph = np.concatenate((fph, tmp2.reshape(1,len(tmp2))))    
            print('FoM: ',fh[-1,0])
            
            iii+=1
    
    # no more reductions
    if iii > 20:
        break
    # max iterations
    if len(fh) > 5000:
        break
    

mom = UpdateLattice(mom,anh[-1,:])
mom.RunMoments()
mom.PlotBeamSize()
















