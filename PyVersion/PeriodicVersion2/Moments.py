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
        ## Reversed Accelerator lattice variable for adjoint equation integration
        self.latticeReversed = None
        ## The starting z position of the lattice [m]
        self.zstart = None
        ## The ending z position of the lattice [m]
        self.zend = None
        
        ## integration step size
        self.h = 0.0001
        ## reversed integration step size
        self.hReversed = -0.0001
        
        ## z position of the moment solution
        self.z = None
        ## Solution to all the 10 moment variables at each z position
        self.y = None
        ## quadrupole strength k [1/m^2] at each z position
        self.k = None
        ## A reverse integration of the moment equations; this should be equal to self.y
        self.yOri = None
        ## Solution to the adjoint moment variables at each z position.
        self.yAdj = None
        ## z position of the adjoint moment solution; this should be equal to self.z
        self.zAdj = None
        ## quadrupole strength k [1/m^2] at each z position; this should match self.k
        self.kAdj = None     
        
        ## normalization factor in the FoM
        self.k0 = 7
        ## O matrix calculated during the moment equation integration at each z value
        self.Omat = None
        ## N matrix calculated during the moment equation integration at each z value
        self.Nmat = None
        
        ## Initial moment values for the integration
        self.initialConditions = initialMoments
        if initialMoments == None:
            self.initialConditions = [
                0.5*(2.2581**2 + 0.2258**2)*1e-6, # Q+
                0.5*(2.2581**2 - 0.2258**2)*1e-6, # Q-
                0.0, # Qx
                0.0, # P+
                0.0, # P-
                0.0, # Px
                (7.0855**2 + 0.70855**2)*1e-6, # E+
                (7.0855**2 - 0.70855**2)*1e-6, # E-
                0.0, # Ex
                0.0, # L
                0.0  # phi , used for solenoid larmor rotations, is always 0 unless we include a solenoid         
                ]
            
        ## Initial moment values for the reverse adjoint integration
        self.initialMomentsAdjoint = None
            
        
    ## Creates the lattice magnet profile that will be integrated. The lengths of dB, qstart, and qend should be equal.
    # @param[in] dB The quadrupole magnet field gradient dB/dm [T/m]
    # @param[in] qstart The z start position for the quadrupole magnet [m]
    # @param[in] qend The z end position for the quadrupole magnet [m]
    # @param[in] qrot The quad rotation angle [radians]
    # @param[in] zstart The z starting position of the lattice [m]
    # @param[in] zend The z ending position of the lattice [m]
    # @param[in] repeat The periodicity of the lattice, e.g. repeat X times
    # @param[in] verbose Verbose diagnostics output
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
                
        
        # create a reversed lattice for backwards integration
        latticeReversed = np.copy(tmpLattice)
        latticeReversed[:,0] = np.flip(tmpLattice[:,1])
        latticeReversed[:,1] = np.flip(tmpLattice[:,0])
        latticeReversed[:,2] = np.flip(tmpLattice[:,2])
        latticeReversed[:,3] = np.flip(tmpLattice[:,3])
        
        # assign to class parameter
        self.lattice = tmpLattice
        self.latticeReversed = latticeReversed
        self.zstart = zstart
        self.zend = zend
        
        if verbose:
            print('Start | End | k strength | quad rotation')
            print(self.lattice)

    ## Run the moment equations through the lattice
    #
    # @param[in] verbose verbose diagnostics print statements
    # @param[out] z longitudinal Z position of the integration. The independent variable.
    # @param[out] y Solutions to the 11 moment variables at each Z position.
    # @param[out] motion Constant of motion value at each Z position.
    # @param[out] k The quadrupole strength value at each Z position.
    def RunMoments(self, verbose=False):
        
        if verbose:
            print("Forward")
        odefunc = lambda z,Y,dbdx,rot : self.__ode_moments(z,Y,dbdx,rot)
        z,y,k,Omat,Nmat = self.__ode3(odefunc, self.h, self.initialConditions, self.lattice, verbose=verbose)
        
        # assign variables
        self.z = z
        self.y = y
        self.k = k
        self.Omat = Omat
        self.Nmat = Nmat
    
        # Constant of Motion 
        motion = self.GetCOM(y)
        
        # calculate initial conditions for adjoint equations
        _,_,_,dFoMpFinal = self.GetFoMAndDFoM(self.y)
        self.initialMomentsAdjoint = np.concatenate((dFoMpFinal,y[:,-1]))
        
        return z,y,motion,k
    
    ## Run the adjoint moment equations backwards through the lattice
    #
    # @param[in] verbose verbose diagnostics print statements
    # @param[out] z longitudinal Z position of the integration. The independent variable.
    # @param[out] y_mom Solutions to the 11 moment variables at each Z position.
    # @param[out] y_adj Solutions to the 11 adjoint moment variables at each Z position
    # @param[out] motion Constant of motion value at each Z position.
    # @param[out] k The quadrupole strength value at each Z position.    
    def RunMomentsAdjoint(self, verbose=False):

        if verbose:
            print("Reverse")
        odefunc = lambda z,Y,dbdx,rot : self.__ode_moments_adjoint(z,Y,dbdx,rot)
        z,y,k,_,_ = self.__ode3(odefunc, self.hReversed, self.initialMomentsAdjoint, self.latticeReversed, verbose=verbose)
        
        y_mom = np.flip(y[11:,:],1)
        y_adj = np.flip(y[0:11,:],1)
        
        # assign variables
        self.yOri = y_mom
        self.yAdj = y_adj
        self.zAdj = np.flip(z)
        self.kAdj = np.flip(k)
        
        # Constant of Motion
        motion = self.GetCOM(y)
        
        return np.flip(z),y_mom,y_adj,motion, np.flip(k)    
        
    ## Recalculate beam properties such as perveance and rigidity
    #
    #
    def CalculateBeamProperties(self):
        # generate some more physics parameters
        self.rigidity, self.perveance = self.__Get_beamridg_and_perv(self.energy, self.current)         
    
    ## Return the constant of motion (COM) from the moment equations. See Tom's notes
    #
    # @param y All 10 moment values at a particular z to calculate the constant of motion
    def GetCOM(self, y):
        L = y[9,:]
        EQ = y[6,:]*y[0,:] + y[7,:]*y[1,:] + y[8,:]*y[2,:]
        PP = y[3,:]**2 + y[4,:]**2 + y[5,:]**2
        motion = EQ + (0.5)*L**2 - (0.5)*PP
        
        return motion
    
    ## Calculate and return the FoM and partial derivatives of the FoM
    #
    # @param[in] Y The moment solutions as a function of z. The output from RunMoments()
    #
    # @param[out] FoM Calculated FoM
    # @param[out] FoMp Calculated FoM pieces
    # @param[out] dFoMp Calculated variation of the FoM pieces
    # @param[out] dFoMpFinal Calculated initial conditions for adjoint integration
    def GetFoMAndDFoM(self, Y=None):
        
        # if no user input, use the set of last calculated Y values from moments run
        if Y is None:
            Y = self.y
            
        # new FoM matching condition
        if 0:
            # figure of merit broken into pieces for ease of reading    
            FoM1 = 0.5*np.sum( (Y[0:3,-1]-Y[0:3,0])**2 )*self.k0**2.
            FoM2 = 0.5*np.sum( (Y[3:6,-1]-Y[3:6,0])**2 )
            FoM3 = 0.5*np.sum( (Y[6:9,-1]-Y[6:9,0])**2 )*self.k0**(-2.)
            FoM4 = 0.5*np.sum( (Y[9,-1]-Y[9,0])**2 )
                
            FoM = FoM1 + FoM2 + FoM3 + FoM4
            FoMp = np.array([FoM1,FoM2,FoM3,FoM4])
            
            # derivative of FoM with respect to Q,P,E,L
            # dF/dQ
            dQ_p = ( Y[0,-1]-Y[0,0] )*self.k0**2.
            dQ_m = ( Y[1,-1]-Y[1,0] )*self.k0**2.
            dQ_x = ( Y[2,-1]-Y[2,0] )*self.k0**2.
            # dF/dP
            dP_p = ( Y[3,-1]-Y[3,0] )
            dP_m = ( Y[4,-1]-Y[4,0] )
            dP_x = ( Y[5,-1]-Y[5,0] )
            # dF/dE
            dE_p = ( Y[6,-1]-Y[6,0] )*self.k0**(-2.)
            dE_m = ( Y[7,-1]-Y[7,0] )*self.k0**(-2.)
            dE_x = ( Y[8,-1]-Y[8,0] )*self.k0**(-2.)
            # dF/dL
            dL = np.abs( Y[9,-1]-Y[9,0] )
            
            dFoMp = np.array([dQ_p,dQ_m,dQ_x,dP_p,dP_m,dP_x,dE_p,dE_m,dE_x,dL,0.0])
            dFoMpFinal = np.array([-dE_p, -dE_m, -dE_x, dP_p, dP_m, dP_x, -dQ_p, -dQ_m, -dQ_x, -dL, 0.0])
            
        # old FoM from FTR
        if 1:
            e1 = 0
            e2 = 1
            # figure of merit kiss 2 no force balance (FoM4)
            f5_tmp = ( Y[6,-1] )
            f4_tmp = ( Y[6,-1] + self.perveance )
            FoM1 = 0.5*( Y[4,-1]**2 + Y[5,-1]**2 ) + 0.5*Y[3,-1]**2
            FoM2 = 0.5*(self.k0**2)*(Y[1,-1]**2 + Y[2,-1]**2)
            FoM3 = 0.5*(self.k0**(-2))*(Y[7,-1]**2+Y[8,-1]**2)
            FoM4 = 0.5*(self.k0**(-2))*e1*(f4_tmp)**2
            FoM5 = (self.k0**(-2))*0.5*e2*(f5_tmp)**2
            
            FoM = FoM1 + FoM2 + FoM3 + FoM4 + FoM5
            FoMp = np.array([FoM1,FoM2,FoM3,FoM4,FoM5])
            
            dP_p = Y[3,-1]
            dP_m = Y[4,-1]
            dP_x = Y[5,-1]
            dE_p = 0.0
            dE_m = -self.k0**(2)*Y[1,-1]
            dE_x = -self.k0**(2)*Y[2,-1] 
            dQ_p = (self.k0**(-2))*e2*(f5_tmp)*(-1) + (self.k0**(-2))*e1*(f4_tmp)*(-1)
            dQ_m = -self.k0**(-2)*Y[7,-1]
            dQ_x = -self.k0**(-2)*Y[8,-1]
            dL = 0.0
            
            dFoMp = np.array([dQ_p,dQ_m,dQ_x,dP_p,dP_m,dP_x,dE_p,dE_m,dE_x,dL,0.0])
            dFoMpFinal = np.array([dQ_p, dQ_m, dQ_x, dP_p, dP_m, dP_x, dE_p, dE_m, dE_x, -dL, 0.0])
        
        return FoM, FoMp, dFoMp, dFoMpFinal
        
    ## Calculate the adjoint variation integral
    #
    # @param[in] Omat The O matrix to use in the calculation
    # @param[in] Nmat The N matrix to use in the calculation
    # @param[out] intVal integral value, see Tom's notes
    def GetAdjointIntegral(self, Omat, Nmat):
        int1 = np.zeros(len(self.z))
        int2 = np.zeros(len(self.z))
        int3 = np.zeros(len(self.z))
        int4 = np.zeros(len(self.z))
        
        for i in range(len(self.z)):
            int1[i] = np.dot( self.yAdj[3:6,i], np.matmul( Omat[0,i], self.y[0:3,i] ) )
            int2[i] = np.dot( self.y[0:3,i], Nmat[0,i] ) * self.yAdj[9,i]
            int3[i] = -np.dot( self.yAdj[0:3,i], np.matmul( Omat[0,i], self.y[3:6,i] ) )
            int4[i] = -np.dot( self.yAdj[0:3,i], Nmat[0,i] ) * self.y[9,i]
            
        intVal = np.trapz(int1+int2+int3+int4, self.z)
        
        return intVal
        
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
    
    ## The main integrator that solves the moment and adjoint equations. Based on a third order Runge-Kutta method
    # This functions goes through each element and integrates the moment equations through a given quadrupole field profile
    #
    # @param[in] F Reference to the function being solved, in this case __ode_moments()
    # @param[in] h The integration step size, see self.h and self.hReversed
    # @paran[in] y0 The initial conditions for the integration.
    # @param[in] The lattice used for the integration
    # @param[in] verbose Print out verbose info while running the integration
    #
    # @param[out] tout The independent variable of the integration
    # @param[out] yout The dependent variable of the integration as a function of tout
    # @param[out] kval The quadrupole strength as a function of tout
    # @param[out] Oval The O matrix values as a function of tout
    # @param[out] Nval The N matrix values as a function of tout
    def __ode3(self, F, h, y0, lattice, verbose=False):
        
        yout = np.array([]).reshape(len(y0),0)
        tout = np.array([])
        kval = np.array([])
        Oval = np.array([], dtype=object).reshape(1,0)
        Nval = np.array([], dtype=object).reshape(1,0)
        
        if verbose:
            print("Integrating elements ... ")
        
        # integrate through each element one at a time, making sure to match the initial conditions at the boundary points of elements
        for j,elem in enumerate(lattice):
        
            if verbose:
                print( str(j+1).zfill(3) + " / " + str(len(lattice[:,0])).zfill(3) )
                
            # start, stop, and k value
            t0 = elem[0]
            t1 = elem[1]
            k = elem[2]
            tsteps = np.arange(t0,t1,h)
            tout = np.concatenate((tout,tsteps))
            kval = np.concatenate((kval,[k]*len(tsteps)))
            
            # initialize output values for this element
            ytmp = np.zeros((len(y0),len(tsteps)+1)) # +1 for initial value
            Otmp = np.empty((1,len(tsteps)), dtype=object)
            Ntmp = np.empty((1,len(tsteps)), dtype=object)
            
            # initial conditions are the very last set of points integrated in the previous element (except for the starting element)
            ytmp[:,0] = yout[:,-1] if (j > 0) else y0
            
            # run rk3 ode solver algorithm through the element
            y = ytmp[:,0]        
            for i,t in enumerate(tsteps):
                t1, Otmp[0,i], Ntmp[0,i] = F(t,y,elem[2], elem[3])
                s1 = h*t1
                t2,_,_ = F(t+h/2.0, y+s1/2.0, elem[2], elem[3])
                s2 = h*t2
                t3,_,_ = F(t+h, y-s1+2.0*s2, elem[2], elem[3])
                s3 = h*t3
                y = y + (s1 + 4.0*s2 + s3)/6.0
                ytmp[:,i+1] = y
                
            # append to main output before moving onto next element        
            yout = np.concatenate((yout,ytmp),1) if (j==0) else np.concatenate((yout,ytmp[:,1:]),1)
            Oval = np.concatenate((Oval,Otmp),axis=1)
            Nval = np.concatenate((Nval,Ntmp),axis=1)
        
        tout = np.concatenate((tout,np.array([tout[-1]+h])))
        kval = np.concatenate((kval,np.array([k])))
        
        # final eval for O matrices
        _, Otmp, Ntmp = F(tout[-1], yout[:,-1], lattice[-1,2], lattice[-1,3])
        Otmp2, Ntmp2 = np.empty((1,1),dtype=object), np.empty((1,1),dtype=object)
        Otmp2[0,0], Ntmp2[0,0] = Otmp, Ntmp
        Oval = np.concatenate((Oval,Otmp2),axis=1)
        Nval = np.concatenate((Nval,Ntmp2),axis=1)
        
        return tout,yout,kval,Oval,Nval
    
    ## Main function that solves the moment equations for all 11 moments.
    #
    # This should not be called directly. Only the integrator ode3() should be calling this as it is integrating through the equations
    #
    # @param z The z step to use [m]
    # @param Y The value of the 11 moments at the particular z value
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
        Y(11) - phi (solenoid related variable)
        
        '''
    
        # quads
        psi = rot
        k_sol = 0.0
        k_quad = dB / self.rigidity
                
        # cq, sq quad calculations         
        cq = np.cos(2.0*Y[10]-2.0*psi)
        sq = np.sin(2.0*Y[10]-2.0*psi)
        
        # space charge calculation
        Q_delta = np.sqrt( Y[0]**2 - Y[1]**2 - Y[2]**2 ) if Y[0]**2 - Y[1]**2 - Y[2]**2 > 0 else 2**(-64.0)
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
        
        return dydt, O_mat, N_mat    

    ## Main function that solves the adjoint moment equations for all 10 moments.
    #
    # This should not be called directly. Only the integrator ode3() should be calling this as it is integrating through the equations
    #
    # @param z The z step to use [m]
    # @param Yt The value of the 22 moments at the particular z value. 
    # This contains the original 11 moments + the 11 adjoint moments for a total of a 22 variable to solve for
    # @param dB The quadrupole gradient [T/m] at the particular z value
    # @param rot The quadrupole rotation angle in radians    
    def __ode_moments_adjoint(self, z, Yt, dB, rot):  
        
        Y = Yt[0:11] # adjoint variables
        Y2 = Yt[11:] # moment variables
        
        # quads
        psi = rot
        k_sol = 0.0
        k_quad = dB / self.rigidity 
                
        # calculations            
        cq = np.cos(2.0*Y2[10]-2.0*psi)
        sq = np.sin(2.0*Y2[10]-2.0*psi)
        
        # space charge stuff
        Q_delta = np.sqrt( Y2[0]**2 - Y2[1]**2 - Y2[2]**2 ) if Y2[0]**2 - Y2[1]**2 - Y2[2]**2 > 0 else 2**(-64.0)
        ab4 = 1.0 / Q_delta; # the 4/ab term in equation
        ca_ab4 = -Y2[1] / ( (Y2[0]+Q_delta)*Q_delta ) # 4c_alpha/ab
        sa_ab4 = -Y2[2] / ( (Y2[0]+Q_delta)*Q_delta ) # 4s_alpha/ab
        
        # calculate O,N matrix
        O_mat,N_mat = self.__calcON(Y2, k_sol, k_quad, cq, sq, ab4, ca_ab4, sa_ab4)
        
        # Calculate special matrix due to space charge variations
        [Mq,Mp,Mn] = self.__GetSCVM(Y2)
        
        # System of 20 equations
        
        # Solve regular envelope equations
        dydt = np.array([
        # dQ/dz
        Y2[3],
        Y2[4],
        Y2[5],
        # dP/dz
        Y2[6] + np.matmul( O_mat[0,:], np.reshape(Y2[0:3],(3,1)) )[0],
        Y2[7] + np.matmul( O_mat[1,:], np.reshape(Y2[0:3],(3,1)) )[0],
        Y2[8] + np.matmul( O_mat[2,:], np.reshape(Y2[0:3],(3,1)) )[0],
        # dE/dz
        np.matmul( O_mat[0,:], np.reshape(Y2[3:6],(3,1)) )[0] + N_mat[0,0]*Y2[9],
        np.matmul( O_mat[1,:], np.reshape(Y2[3:6],(3,1)) )[0] + N_mat[1,0]*Y2[9],
        np.matmul( O_mat[2,:], np.reshape(Y2[3:6],(3,1)) )[0] + N_mat[2,0]*Y2[9],
        # dL/dz
        -1.0*np.matmul( N_mat.T, Y2[0:3] )[0],
        # dphi/dz
        -1.0*k_sol/2.0
        ])
        
        # dE dot 1 (y)
        dedot1 = np.matmul( np.reshape(Y[3:6], (1,3)), Mq )[0] + \
            Y[9]*np.matmul( np.reshape(Y2[0:3], (1,3)), Mn )[0] - \
            np.matmul( np.reshape(Y[0:3], (1,3)), Mp )[0] - \
            Y2[9]*np.matmul( np.reshape(Y[0:3], (1,3)), Mn )[0] 
        # dE dot 2 (y)
        dedot2 = np.array([0,0,0])
        # dE dot total (y)
        dedot = dedot1 + dedot2
        # dQ dot (y)
        dqdot = np.array([0,0,0])
        if z == 0.0:
            dqdot = np.array([2*Y2[0],0,0])
        # dP dot (y)
        dpdot = np.array([0,0,0])
        # dL dot (y)
        dldot = 0 
        
        # Solve adjoint equations
        dydt2 = np.array([
        #dQ/dz (y)
        Y[3] + dqdot[0],
        Y[4] + dqdot[1],
        Y[5] + dqdot[2],
        # dP/dz (y)
        Y[6] + np.matmul( O_mat[0,:], np.reshape(Y[0:3],(3,1)) )[0] + dpdot[0],
        Y[7] + np.matmul( O_mat[1,:], np.reshape(Y[0:3],(3,1)) )[0] + dpdot[1],
        Y[8] + np.matmul( O_mat[2,:], np.reshape(Y[0:3],(3,1)) )[0] + dpdot[2],
        # dE/dz (y) 
        np.matmul( O_mat[0,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[0,0]*Y[9] + dedot[0],
        np.matmul( O_mat[1,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[1,0]*Y[9] + dedot[1],
        np.matmul( O_mat[2,:], np.reshape(Y[3:6],(3,1)) )[0] + N_mat[2,0]*Y[9] + dedot[2],
        # dL/dz    
        -1.0*np.matmul( N_mat.T, Y[0:3] )[0] + dldot,   
        # dphi/dz
        -1.0*k_sol/2.0
        ])
        
        return np.concatenate((dydt2,dydt)), O_mat, N_mat
    
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
    
    ## Recalculate the O and N matrices using the current lattice
    #
    # @param[out] Onew The new O matrix
    # @param[out] Nnew The new N matrix
    def RecalculateON(self):
        
        # initialize new O matrix arrays
        Onew = np.empty((1,len(self.z)), dtype=object)
        Nnew = np.empty((1,len(self.z)), dtype=object)

        for i, (zi,Y) in enumerate(zip(self.z, self.y.T)):

            # find the lattice element we are in
            for elem in self.lattice:
                if zi <= elem[1]:
                    # quads
                    psi = elem[3]
                    k_sol = 0.0
                    k_quad = elem[2] / self.rigidity
                    break;

            # cq, sq quad calculations         
            cq = np.cos(2.0*Y[10]-2.0*psi)
            sq = np.sin(2.0*Y[10]-2.0*psi)
            
            # space charge calculation
            Q_delta = np.sqrt( Y[0]**2 - Y[1]**2 - Y[2]**2 ) if Y[0]**2 - Y[1]**2 - Y[2]**2 > 0 else 2.0**(-64.0)
            ab4 = 1.0 / Q_delta; # the 4/ab term in equation
            ca_ab4 = -Y[1] / ( (Y[0]+Q_delta)*Q_delta ) # 4c_alpha/ab
            sa_ab4 = -Y[2] / ( (Y[0]+Q_delta)*Q_delta ) # 4s_alpha/ab
            
            # calculate O,N matrix
            O_mat,N_mat = self.__calcON(Y, k_sol, k_quad, cq, sq, ab4, ca_ab4, sa_ab4)

            Onew[0,i] = O_mat
            Nnew[0,i] = N_mat
            
        return Onew, Nnew    
    
    ## Calculate the space charge variation matrices used in the adjoint equations
    #
    # @param Y A vector of the 10 moment values
    def __GetSCVM(self, Y):
        k_perv = self.perveance
        
        Q_delta = np.sqrt( Y[0]**2 - Y[1]**2 - Y[2]**2 )
        Q_deltaplus = Q_delta + Y[0]
    
        V1 = -(k_perv/(Q_delta**2))*np.array([Y[3]-Y[1]*Y[4]/Q_deltaplus-Y[2]*Y[5]/Q_deltaplus,
            -Y[1]*Y[3]/Q_deltaplus+Y[4],
            -Y[2]*Y[3]/Q_deltaplus+Y[5]]).reshape(3,1)
        U1t = (1./Q_delta)*np.array([Y[0], -Y[1], -Y[2]]).reshape(1,3)
    
        V2 = (k_perv/(Q_delta*(Q_deltaplus**2)))*np.array([Y[1]*Y[4]+Y[2]*Y[5],Y[1]*Y[3],Y[2]*Y[3]]).reshape(3,1)
        U2t = U1t + np.array([1.,0,0]).reshape(1,3)
    
        V3 = -(k_perv/(Q_delta*Q_deltaplus))*np.array([Y[4], Y[3], 0]).reshape(3,1)
        V4 = -(k_perv/(Q_delta*Q_deltaplus))*np.array([Y[5], 0, Y[3]]).reshape(3,1)
        U3t = np.array([0, 1., 0]).reshape(1,3)
        U4t = np.array([0,0,1.]).reshape(1,3)
    
        Mp = V1*U1t + V2*U2t + V3*U3t + V4*U4t
        ####
    
        W1 = -(k_perv/Q_delta)*np.array([1., Y[1]/Q_deltaplus, Y[2]/Q_deltaplus]).reshape(3,1)
        W2 = (k_perv/(Q_delta*(Q_deltaplus**2)))*np.array([Y[1]**2+Y[2]**2,Y[1]*Y[0],Y[2]*Y[0]]).reshape(3,1)
        W3 = -(k_perv/(Q_delta*Q_deltaplus))*np.array([Y[1],Y[0],0]).reshape(3,1)
        W4 = -(k_perv/(Q_delta*Q_deltaplus))*np.array([Y[2],0,Y[0]]).reshape(3,1)
    
        Mq = W1*U1t + W2*U2t + W3*U3t + W4*U4t
        ####
    
        X1 = -k_perv*np.array([0, Y[2]/(Q_deltaplus*(Q_delta**2)), -Y[1]/(Q_deltaplus*(Q_delta**2))]).reshape(3,1)
        X2 = -k_perv*np.array([0, Y[2]/(Q_delta*(Q_deltaplus**2)), -Y[1]/(Q_delta*(Q_deltaplus**2))]).reshape(3,1)
        X3 = k_perv*np.array([0,0,-1./(Q_delta*Q_deltaplus)]).reshape(3,1)
        X4 = k_perv*np.array([0,1./(Q_delta*Q_deltaplus),0]).reshape(3,1)
    
        Mn = X1*U1t + X2*U2t + X3*U3t + X4*U4t
        
        return Mq,Mp,Mn    
    
    ## Convert cartesian coordinates to the larmor frame
    #
    # @param[in] y a set of moments in the cartesian frame
    # @param[in] dphi solenoid/larmor rotation variable
    # @param[out] yl moments in the larmor frame
    def Cart2Larmor(self, y, dphi):
        yl = np.zeros(len(y))
        
        lc = np.cos(2*y[10])
        ls = np.sin(2*y[10])
        
        # Q+, Q-, Qx
        yl[0] = y[0]
        yl[1] = y[1]*lc + y[2]*ls
        yl[2] = y[2]*lc - y[1]*ls
        
        # P+, P-, Px
        yl[3] = y[3]
        yl[4] = y[4]*lc + y[5]*ls + 2*dphi*yl[2]
        yl[5] = -y[4]*ls + y[5]*lc - 2*dphi*yl[1]
        
        # L
        yl[9] = y[9] - 2*dphi*yl[0]
        
        # E+, E-, Ex
        # need to finish this
    
    ## Plot the beam size after running a moment solution
    #
    #
    def PlotBeamSize(self, useAdjointVariables=False):
        
        plt.figure()
        
        if useAdjointVariables:
            plt.plot(self.zAdj, self.yAdj[0,:] + self.yAdj[1,:], label='x^2')
            plt.plot(self.zAdj, self.yAdj[0,:] - self.yAdj[1,:], label='y^2')        
            plt.plot(self.zAdj, self.kAdj / 1e2*1e-6*0.2 / self.rigidity, color='k', label='quad strength')
        else:
            plt.plot(self.z, self.y[0,:] + self.y[1,:], label='x^2')
            plt.plot(self.z, self.y[0,:] - self.y[1,:], label='y^2')        
            plt.plot(self.z, self.k / 1e2*1e-6*0.2 / self.rigidity, color='k', label='quad strength')
        
        plt.grid(True)
        plt.legend()
        

        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    