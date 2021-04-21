import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from tabulate import tabulate
    
    
def Get_beamridg_and_perv(energy=5e3,current=0.0):
    '''
    Grab beam ridgidity
    '''
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


def ode3(F,h,y0,lattice,verbose_f):
    '''
    3rd order modified classical Runge-Kutta ODE solver
    
    This functions goes through each element and integrates the moment equations through a given quadrupole field profile
    '''
        
    yout = np.array([]).reshape(len(y0),0)
    tout = np.array([])
    kval = np.array([])
    
    # integrate through each element one at a time, making sure to match the initial conditions at the boundary points of elements
    for j,elem in enumerate(lattice):
    
        if verbose_f:
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
            t1 = F(t,y,elem[2])
            s1 = h*t1
            t2 = F(t+h/2.0, y+s1/2.0, elem[2])
            s2 = h*t2
            t3 = F(t+h, y-s1+2.0*s2, elem[2])
            s3 = h*t3
            y = y + (s1 + 4.0*s2 + s3)/6.0
            ytmp[:,i+1] = y
        
        # append to main output before moving onto next element        
        yout = np.concatenate((yout,ytmp),1) if (j==0) else np.concatenate((yout,ytmp[:,1:]),1)

    tout = np.concatenate((tout,np.array([tout[-1]+h])))
    kval = np.concatenate((kval,np.array([k])))  
    return tout,yout,kval
    
def getLatticeKvsZ(lattice,h):
    '''
    returns K as a function of z
    '''   
    
    kval = np.array([])
    
    for j,elem in enumerate(lattice):
        # start, stop, and k value
        t0 = elem[0]
        t1 = elem[1]
        k = elem[2]
        tsteps = np.arange(t0,t1,h)
        kval = np.concatenate((kval,[k]*len(tsteps)))   
    kval = np.concatenate((kval,np.array([k])))
    
    return kval          
    
def CreateLatticeProfile(amplitude=1.0,qlength=0.1,dlength=0.1,polarity=[1,-1],repeat=1,verbose=False):
    '''
    amplitude - quadrupole strength amplitude
    qlength - quadrupole length
    dlength - drift length inbetween quads
    repeat - number of times to repeat
    polarity - quad polarity 
    
    the length of qlength, dlength, amplitude, polarity should be equal
    
    The arrays should follow the lattice pattern: drift length, quad length, drift length, quad length, etc...
    note that either drift length or quad length can be zero. This way you can make patterns like drift-quad-quad-drift etc.
    '''

    # how many quads
    amplitude = np.tile(amplitude,repeat)
    qlength = np.tile(qlength,repeat)
    dlength = np.tile(dlength,repeat)
    polarity = np.tile(polarity,repeat)        
    
    # organize elements into three arrays, element start position, end position, and k value
    # organize into two arrays for the quad starting and ending positions
    elemLocations = np.zeros(len(qlength) + len(dlength))
             
    ii = 0
    for i in range(len(qlength)):
        elemLocations[ii] = dlength[i]
        ii+=1
        elemLocations[ii] = qlength[i]
        ii+=1
    elemLocations = np.cumsum(elemLocations) # this list should always be even length, since we need a dlength for every qlength
    
    qStartLocations = elemLocations[0::2]
    qEndLocations = elemLocations[1::2]
    dStartLocations = np.concatenate(([0],elemLocations[1::2][0:-1]))
    dEndLocations = elemLocations[0::2]
    
    ##########################
    if verbose:
        print("\n")
        print(amplitude)
        print(qlength)
        print(dlength)
        print(polarity)
        print(qStartLocations)
        print(qEndLocations)
        print(dStartLocations)
        print(dEndLocations)        
        print("\n")
    ##########################
    
    elemLatticeInfo = np.zeros((len(elemLocations),3))
    
    #
    # element # | start location | end location | k value
    ii=0
    for i in range(len(qStartLocations)):
        elemLatticeInfo[ii,0] = dStartLocations[i]
        elemLatticeInfo[ii,1] = dEndLocations[i]
        elemLatticeInfo[ii,2] = 0.0 # drift has no k value
        ii+=1
        elemLatticeInfo[ii,0] = qStartLocations[i]
        elemLatticeInfo[ii,1] = qEndLocations[i]
        elemLatticeInfo[ii,2] = amplitude[i]*polarity[i]
        ii+=1                  
       
    return elemLatticeInfo
    
def PlotLatticeProfile(lattice):                     
    N,_ = lattice.shape
    for ii,elem in enumerate(lattice):
        if (elem[0] != elem[1]):
            plt.plot([elem[0],elem[1]],[elem[2],elem[2]],color='k')
            # connect previous element to current element
            if (ii > 0 and ii < N): # no connecting first or last element
                if (lattice[ii-1,0] != lattice[ii-1,1]):
                    plt.plot([lattice[ii-1,1],lattice[ii,0]],[lattice[ii-1,2],lattice[ii,2]],color='k')





    
