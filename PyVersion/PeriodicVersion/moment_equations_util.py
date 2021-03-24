import numpy as np
import scipy as sp



def Initial_conditions():
    '''
    return starting conditions for the moment equations (based off Santiago's values)
    can be updated whenever.
    '''
    Q_plus = 0.5*(2.2581**2*1e-6 + 0.2258**2*1e-6)
    Q_minus = 0.0
    Q_x = 0.0
    P_plus = 0
    P_minus = 0
    P_x = 0
    E_plus = 0.0
    E_minus = 0.0
    E_x = 0
    L = 0
    phi = 0
   
    return np.array([Q_plus,Q_minus,Q_x,P_plus,P_minus,P_x,E_plus,E_minus,E_x,L,phi])
    
    
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


def ode3(F,t0,h,tfinal,y0,k,verbose_f):
    '''
    third order modified classical Runge-Kutta ODE solver
    '''
    y = y0
    tsteps = np.arange(t0,tfinal,h)
    ii=0 # counter for quadprofile
        
    yout = np.zeros((len(y0),len(tsteps)+1))
    yout[:,0] = y
    
    
    # if verbose, same algorithm but with integration progress print outs
    if verbose_f:
        N = len(tsteps)
        NN = int(N/10.0)
        for i,t in enumerate(tsteps):
            t1 = F(t,y,k[ii])
            s1 = h*t1
            t2 = F(t+h/2.0, y+s1/2.0, k[ii+1])
            s2 = h*t2
            t3 = F(t+h, y-s1+2.0*s2, k[ii+2])
            s3 = h*t3
            y = y + (s1 + 4.0*s2 + s3)/6.0
            yout[:,i+1] = y 
            ii+=2        
            if i % NN == 0:
                print(str(int(1.0*i/NN*10)) + '%')
    else:
        for i,t in enumerate(tsteps):
            t1 = F(t,y,k[ii])
            s1 = h*t1
            t2 = F(t+h/2.0, y+s1/2.0, k[ii+1])
            s2 = h*t2
            t3 = F(t+h, y-s1+2.0*s2, k[ii+2])
            s3 = h*t3
            y = y + (s1 + 4.0*s2 + s3)/6.0
            yout[:,i+1] = y
            ii+=2
            
    return yout
    
    
def CreateQuadProfile(amplitude=1.0,qlength=0.1,dlength=0.1,numQuads=2,polarity=[1,-1],stepsize=0.01,verbose=False):
    '''
    amplitude - quadrupole strength amplitude
    qlength - quadrupole length
    dlength - drift length inbetween quads
    numQuads - number of quads
    polarity - quad polarity
    stepsize -    
    '''

    # how many quads
    numQuads = (int)(numQuads)
    
    ########################## setup arrays if needed
    # setup amplitude array based on total quads
    if (not isinstance(amplitude,(list,np.ndarray))):
        amplitude = np.array([amplitude]*numQuads)
    else:
        if (len(amplitude) != numQuads):
            # if the length of the array is shorter than # of quads, just repeat the pattern
            if (isinstance(amplitude,list)):
                amplitude = np.array(amplitude)
            amplitude = np.tile(amplitude,(1+(int)(numQuads/len(amplitude))))            
    # setup qlength array based on total quads
    if (not isinstance(qlength,(list,np.ndarray))):
        qlength = np.array([qlength]*numQuads)
    else:
        if (len(qlength) != numQuads):
            # if the length of the array is shorter than # of quads, just repeat the pattern
            if (isinstance(qlength,list)):
                qlength = np.array(qlength)
            qlength = np.tile(qlength,(1+(int)(numQuads/len(qlength))))          
    # setup dlength array based on total quads 
    if (not isinstance(dlength,(list,np.ndarray))):
        dlength = np.array([dlength]*numQuads) 
    else:
        if (len(dlength) != numQuads):
            # if the length of the array is shorter than # of quads, just repeat the pattern
            if (isinstance(dlength,list)):
                dlength = np.array(dlength)
            dlength = np.tile(dlength,(1+(int)(numQuads/len(dlength))))         
    # setup polarity array based on total quads
    if (not isinstance(polarity,(list,np.ndarray))):
        polarity = np.array([polarity]*numQuads)   
    else:
        if (len(polarity) != numQuads):
            # if the length of the array is shorter than # of quads, just repeat the pattern
            if (isinstance(polarity,list)):
                polarity = np.array(polarity)
            polarity = np.tile(polarity,(1+(int)(numQuads/len(polarity))))           
                
    # organize into two arrays for the quad starting and ending positions
    qLocations = np.zeros(numQuads*2)           
    ii = 0
    for i in range(numQuads):
        qLocations[ii] = dlength[i]
        ii+=1
        qLocations[ii] = qlength[i]
        ii+=1
    print(qLocations)
    qLocations = np.cumsum(qLocations)
    qStartLocations = qLocations[0::2]
    qEndLocations = qLocations[1::2]
    
    ##########################
    if verbose:
        print("\n")
        print(amplitude)
        print(qlength)
        print(dlength)
        print(polarity)
        print(qStartLocations)
        print(qEndLocations)
        print("\n")
    ##########################        

    # setup the length of our quadrupole pattern
    zi = 0.0
    ze = np.sum(qlength[0:numQuads]) + np.sum(dlength[0:numQuads])
    z = np.arange(zi,round(ze,4)+stepsize,stepsize)
    N = len(z)
    k = np.zeros(N)
    
    ## iterate through z and fill in the quad values    
    ii=0 # current quad index values
    for i,zv in enumerate(z):
        if (zv >= qStartLocations[ii] and zv <= qEndLocations[ii]):
            k[i] = amplitude[ii]*polarity[ii]
        elif (zv >= qEndLocations[ii]):
            ii+=1      
            if  (ii < numQuads and zv >= qStartLocations[ii] and zv <= qEndLocations[ii]):
                k[i] = amplitude[ii]*polarity[ii]
    # last value is the same as the first value (periodic)
    k[-1] = k[0]
       
    return z,k                    
            
    
    







    
