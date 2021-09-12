# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 12:21:15 2021

@author: levon
"""

import numpy as np
import matplotlib.pyplot as plt

## The Moments class solves the moment equations and adjoint equations
#
#
class MomentSolver:
    
    ## instance attributes
    # @param energy Beam energy
    # @param current Beam current
    # @param pipeRadius vacuum pipe radius used in image force calculations. A value of 0.0 means the image force calculation is ignored
    # @param initialMoments Initial moment values for the integration
    def __init__(self, energy=10e3, current=0.0, pipeRadius=0.0, initialMoments=None):
        
        # grab physics parameters
        ## Particle beam energy [eV]
        self.energy = energy
        ## Particle beam current [A]
        self.current = current
        ## Vacuum pipe radius [m]
        self.pipeRadius = pipeRadius
        ## Beam Rigidity
        self.rigidity = None
        ## Beam perveance, a measure of space charge
        self.perveance = 0.0
        
        # generate some more physics parameters
        self.rigidity, self.perveance = self.__Get_beamridg_and_perv(self.energy, self.current) 
        
        ## Accelerator lattice variable. See CreateLatticeProfile()
        self.lattice = None
        ## The starting z position of the lattice [m]
        self.zstart = None
        ## The ending z position of the lattice [m]
        self.zend = None
        
        ## integration step size.
        self.h = 0.0001
        
        ## z position of the moment solution
        self.z = None
        ## Solution to all the 10 moment variables at each z position
        self.y = None
        ## quadrupole strength k [1/m^2] at each z position
        self.k = None
        
        ## Initial moment values for the integration
        if initialMoments == None:
            self.initialConditions = [
                0.5*(2.2581**2*1e-6 + 0.2258**2*1e-6), # Q+
                0.5*(2.2581**2*1e-6 - 0.2258**2*1e-6), # Q-
                0.0, # Qx
                0.0, # P+
                0.0, # P-
                0.0, # Px
                (7.0855**2*1e-6 + 0.70855**2*1e-6), # E+
                (7.0855**2*1e-6 - 0.70855**2*1e-6), # E-
                0.0, # Ex
                0.0, # L
                0.0  # phi , used for solenoid larmor rotations, is always 0 unless we include a solenoid         
                ]          
            
        
    ## Creates the lattice magnet profile that will be integrated. The lengths of dB, qstart, and qend should be equal.
    # @param dB The quadrupole magnet field gradient dB/dm [T/m]
    # @param qstart The z start position for the quadrupole magnet [m]
    # @param qend The z end position for the quadrupole magnet [m]
    # @param qrot The quad rotation angle [radians]
    # @param zstart The z starting position of the lattice [m]
    # @param zend The z ending position of the lattice [m]
    # @param repeat The periodicity of the lattice, e.g. repeat X times
    # @param verbose Verbose diagnostics output
    def CreateLatticeProfile(self, dB=[1.0], qstart=[0.0], qend=[0.1], qrot=[0.0], zstart=0, zend=1, repeat=1, verbose=False):
        
        Nquads = len(dB) # number of quads
        numCols = 4 # q start, q end, k strength, q rot
        
        tmpLattice = np.zeros((1, numCols))
        
        # drift from z start to q1 start
        tmpLattice[0,:] = [zstart, qstart[0], 0.0, 0.0]
        
        # loop through each quad and add to the lattice
        for i in range(Nquads):
            
            # add quad to lattice
            tmpQuadPos = [[qstart[i], qend[i], dB[i], qrot[i]]]
            tmpLattice = np.concatenate((tmpLattice, tmpQuadPos), axis=0)
            
            # if this is not the last quad yet, add in the drift to the next quad
            if (i+1) < Nquads:
                tmpDrftPos = [[qend[i], qstart[i+1], 0.0, 0.0]]
                tmpLattice = np.concatenate((tmpLattice, tmpDrftPos), axis=0)
                
        # drfit from last q end to z end
        tmpDrftPos = [[qend[-1], zend, 0.0, 0.0]]
        tmpLattice = np.concatenate((tmpLattice, tmpDrftPos), axis=0)
        
        # are we repeating the lattice?
        if repeat > 1:
            # make a copy of the lattice
            copyTmpLattice = np.copy(tmpLattice)
            for ii in range(repeat-1):
                # grab the ending location of the last element
                zendPos = tmpLattice[-1,1]
                
                # create a copy of the lattice
                newLatticePiece = np.copy(copyTmpLattice)
                
                # adjust the starting and ending positions of this new piece
                newLatticePiece[:,0:2] += zendPos
                
                # append to the tmp lattce
                tmpLattice = np.concatenate((tmpLattice, newLatticePiece), axis=0)
                
                
        # assign to class parameter
        self.lattice = tmpLattice
        self.zstart = zstart
        self.zend = zend
        
        if verbose:
            print('Start | End | k strength | quad rotation')
            print(self.lattice)

    ## Run the moment equations through the lattice
    #
    # @param verbose verbose diagnostics print statements
    def RunMoments(self, verbose=False):
        
        odefunc = lambda z,Y,dbdx,rot : self.__ode_moments(z,Y,dbdx,rot)
        z,y,k = self.__ode3(odefunc, verbose=verbose)
        
        # assign variables
        self.z = z
        self.y = y
        self.k = k
    
        # Constant of Motion %
        motion = self.GetCOM(y)
        
        return z,y,motion,k  
    
    
    ## Calculates and returns the beam ridgidity and perveance
    #
    # @param energy The beam energy [eV]
    # @param current The beam current [A]
    def __Get_beamridg_and_perv(self, energy=5e3, current=0.0):
        
        # Parameters
        e         = 1.60217733E-19 #C
        m         = 9.1093897E-31 #kg
        Energy    = energy # eV
        c         = 2.997924E8 # m/s
    
        gamma     = 1+((Energy)/(510998.9461));
        beta      = np.sqrt((gamma*gamma)-1)/gamma
        v         = beta*c
        bg        = beta*gamma
        rho       = bg*c*(m/e) 
        
        k_perv = (1.0/(4.0*np.pi))*(c*377.0*current) / (m*v**3*gamma**3/e);   
        
        return rho,k_perv
    
    ## Return the constant of motion (COM) from the moment equations. See Tom's notes
    #
    # @param y All 10 moment values at a particular z to calculate the constant of motion
    def GetCOM(self, y):
        L = y[9,:]
        EQ = y[6,:]*y[0,:] + y[7,:]*y[1,:] + y[8,:]*y[2,:]
        PP = y[3,:]**2 + y[4,:]**2 + y[5,:]**2
        motion = EQ + (0.5)*L**2 - (0.5)*PP
        
        return motion
    
    ## The main integrator that solves the moment and adjoint equations. Based on a third order Runge-Kutta method
    # This functions goes through each element and integrates the moment equations through a given quadrupole field profile
    #
    # @param F Reference to the function being solved, in this case __ode_moments()
    # @param verbose Print out verbose info while running the integration
    def __ode3(self, F, verbose=False):
        
        y0 = self.initialConditions
        h = self.h
        
        yout = np.array([]).reshape(len(y0),0)
        tout = np.array([])
        kval = np.array([])
        
        # integrate through each element one at a time, making sure to match the initial conditions at the boundary points of elements
        for j,elem in enumerate(self.lattice):
        
            if verbose:
                print("Integrating element: "+str(j+1))
                
            # start, stop, and k value
            t0 = elem[0]
            t1 = elem[1]
            k = elem[2]
            tsteps = np.arange(t0,t1,h)
            tout = np.concatenate((tout,tsteps))
            kval = np.concatenate((kval,[k]*len(tsteps)))
            
            # initialize output values for this element
            ytmp = np.zeros((len(y0),len(tsteps)+1)) # +1 for initial value
            # initial conditions are the very last set of points integrated in the previous element (except for the starting element)
            ytmp[:,0] = yout[:,-1] if (j > 0) else y0
            
            # run rk3 ode solver algorithm through the element
            y = ytmp[:,0]        
            for i,t in enumerate(tsteps):
                t1 = F(t,y,elem[2], elem[3])
                s1 = h*t1
                t2 = F(t+h/2.0, y+s1/2.0, elem[2], elem[3])
                s2 = h*t2
                t3 = F(t+h, y-s1+2.0*s2, elem[2], elem[3])
                s3 = h*t3
                y = y + (s1 + 4.0*s2 + s3)/6.0
                ytmp[:,i+1] = y
            
            # append to main output before moving onto next element        
            yout = np.concatenate((yout,ytmp),1) if (j==0) else np.concatenate((yout,ytmp[:,1:]),1)
    
        tout = np.concatenate((tout,np.array([tout[-1]+h])))
        kval = np.concatenate((kval,np.array([k])))  
        return tout,yout,kval
    
    ## Main function that solves the moment equations for all 10 moments.
    #
    # This should not be called directly. Only the integrator ode3() should be calling this as it is integrating through the equations
    #
    # @param z The z step to use [m]
    # @param Y The value of the 10 moments at the particular z value
    # @param dB The quadrupole gradient [T/m] at the particular z value
    # @param rot The quadrupole rotation angle in radians
    def __ode_moments(self, z, Y, dB, rot):
        '''
        Main function that solves Tom's moment equations
        
        Y(1) - Q+
        Y(2) - Q-
        Y(3) - Qx
        Y(4) - P+
        Y(5) - P-
        Y(6) - Px
        Y(7) - E+
        Y(8) - E-
        Y(9) - Ex
        Y(10) - L
        
        '''
    
        # quads
        psi = rot
        k_sol = 0.0
        k_quad = dB / self.rigidity
                
        # cq, sq quad calculations         
        cq = np.cos(2.0*Y[10]-2.0*psi)
        sq = np.sin(2.0*Y[10]-2.0*psi)
        
        # space charge calculation
        Q_delta = np.sqrt( Y[0]**2 - Y[1]**2 - Y[2]**2 ) if Y[0]**2 - Y[1]**2 - Y[2]**2 >= 0 else 2**(-64.0)
        ab4 = 1.0 / Q_delta; # the 4/ab term in equation
        ca_ab4 = -Y[1] / ( (Y[0]+Q_delta)*Q_delta ) # 4c_alpha/ab
        sa_ab4 = -Y[2] / ( (Y[0]+Q_delta)*Q_delta ) # 4s_alpha/ab
        
        # calculate O,N matrix
        O_mat,N_mat = self.__calcON(Y, k_sol, k_quad, cq, sq, ab4, ca_ab4, sa_ab4)
        
        # System of 10 equations
        dydt = np.array([
        # dQ/dz
        Y[3],
        Y[4],
        Y[5],
        # dP/dz
        Y[6] + np.matmul( O_mat[0,:], np.reshape(Y[0:3],(3,1)) )[0],
        Y[7] + np.matmul( O_mat[1,:], np.reshape(Y[0:3],(3,1)) )[0],
        Y[8] + np.matmul( O_mat[2,:], np.reshape(Y[0:3],(3,1)) )[0],
        # dE/dz
        np.matmul( O_mat[0,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[0,0]*Y[9],
        np.matmul( O_mat[1,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[1,0]*Y[9],
        np.matmul( O_mat[2,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[2,0]*Y[9],
        # dL/dz
        -1.0*np.matmul( N_mat.T, Y[0:3] )[0],
        # dphi/dz
        -1.0*k_sol/2.0
        ])
        
        return dydt    
    
    ## Helper function for __ode_moments() that calculates the O and N matrices in the moment equations
    # See Tom's notes for the definitions of these matrices
    #
    def __calcON(self, Y, k_sol, k_quad, cq, sq, ab4, ca_ab4, sa_ab4): 
        # pipe radius calculation
        pipe_constant = (8*self.perveance/self.pipeRadius**4) if self.pipeRadius != 0 else 0 # zero pipe radius leads to infinity
        
        # O and N matrix calculations
        O_mat = np.array([[-k_sol**2/2.0 + ab4*self.perveance, 2.0*k_quad*cq + ca_ab4*self.perveance, -2.0*k_quad*sq + sa_ab4*self.perveance],
            [2.0*k_quad*cq + ca_ab4*self.perveance, -k_sol**2/2.0 + ab4*self.perveance, 0],
            [-2.0*k_quad*sq + sa_ab4*self.perveance, 0, -k_sol**2/2.0 + ab4*self.perveance]]) + pipe_constant*np.array([ [0, -Y[1], -Y[2]],[-Y[1], 0, 0],[-Y[2], 0, 0] ])
    
        N_mat = np.array([[0], [2.0*k_quad*sq - sa_ab4*self.perveance], [2.0*k_quad*cq + ca_ab4*self.perveance]])
        - pipe_constant*np.array([[0], [-Y[2]], [Y[1]]]) # pipe radius image forces addition  
        
        return O_mat,N_mat    
    
    ## Plot the beam size after running a moment solution
    #
    #
    def PlotBeamSize(self):
        
        plt.figure()
        
        plt.plot(self.z, self.y[0,:] + self.y[1,:], label='x^2')
        plt.plot(self.z, self.y[0,:] - self.y[1,:], label='y^2')
        
        plt.plot(self.z, self.k / 1e2*1e-6*0.2 / self.rigidity, color='k', label='quad strength')
        
        plt.grid(True)
        plt.legend()
        #plt.show()
        

        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    