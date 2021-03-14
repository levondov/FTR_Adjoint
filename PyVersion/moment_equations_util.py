import numpy as np
import scipy as sp



def initial_conditions():
    '''
    return starting conditions for the moment equations (based off Santiago's values)
    can be updated whenever.
    '''
    Q_plus = 0.5*(2.2581**2*1e-6 + 0.2258**2*1e-6)
    Q_minus = 0.5*(2.2581**2*1e-6 - 0.2258**2*1e-6)
    Q_x = 0.0
    P_plus = 0
    P_minus = 0
    P_x = 0
    E_plus = (7.0855**2*1e-6 + 0.70855**2*1e-6)
    E_minus = (7.0855**2*1e-6 - 0.70855**2*1e-6)
    E_x = 0
    L = 0
    phi = 0
   
    return np.array([Q_plus,Q_minus,Q_x,P_plus,P_minus,P_x,E_plus,E_minus,E_x,L,phi])
    
    
def initial_magnet_profiles():
    '''
    Default starting magnet profiles, based on a thin lens model. If you want to change the profiles DO NOT do it here. Instead use the "param_scale" parameter as input to the many functions
    ''' 
    params = np.array(
        [0.213,0.863,15e-4, # solenoid start, length, strength
        .00425,.0001,-18.236, # quad 1 start, length, strength
        0.10655,0.0001,21.3640, # quad 2 start, length, strength
        0.20895,0.0001,-18.236, # quad 3 start, length, strength
        45*np.pi/180,45*np.pi/180,45*np.pi/180]) # quad 1,2,3 rotations
        
    return params
    
def get_beamridg_and_perv(energy=5e3,current=0.0):
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


def ode3(F,t0,h,tfinal,y0,verbose=False):
    '''
    third order classical Runge-Kutta ODE solver
    '''
    y = y0
    tsteps = np.arange(t0,tfinal-h,h)
    
    # for extra params
    ksol = np.zeros(len(tsteps)+1)
    kquad = np.zeros(len(tsteps)+1)
        
    yout = np.zeros((len(y0),len(tsteps)+1))
    yout[:,0] = y
    if verbose:
        N = len(tsteps)
        NN = int(N/10.0)
        for i,t in enumerate(tsteps):
            t1,sol,quad = F(t,y)
            s1 = h*t1
            t2,_,_ = F(t+h/2.0, y+s1/2.0)
            s2 = h*t2
            t3,_,_ = F(t+h, y-s1+2.0*s2)
            s3 = h*t3
            y = y + (s1 + 4.0*s2 + s3)/6.0
            yout[:,i+1] = y
            ksol[i+1] = sol
            kquad[i+1] = quad            
            if i % NN == 0:
                print(str(int(1.0*i/NN*10)) + '%')
    else:
        for i,t in enumerate(tsteps):
            t1,sol,quad = F(t,y)
            s1 = h*t1
            t2,_,_ = F(t+h/2.0, y+s1/2.0)
            s2 = h*t2
            t3,_,_ = F(t+h, y-s1+2.0*s2)
            s3 = h*t3
            y = y + (s1 + 4.0*s2 + s3)/6.0
            yout[:,i+1] = y   
            ksol[i+1] = sol
            kquad[i+1] = quad             
            
    return yout,ksol,kquad
    
    
    
def get_params(a, param_scale={'ql': [1.0,1.0,1.0], 'qs': [1.0,1.0,1.0], 'qs_off': [1.0,1.0,1.0], 'qa': [1.0,1.0,1.0], 'sl': 1.0, 'ss': 1.0, 'ss_off': 1.0}):
    '''
    Given the vector 'a' , return the magnet param setting for the moment equations
    a = (s start, s strength, q1 start, q1 strength, q2 start, q2 strength, q3 start, q3 strength, q1 angle, q2 angle, q3 angle

    a = [  0     1        2       3      4       5       6       7      8      9    10 
    a = (s_sta, s_str, q1_sta, q1_sta, q2_sta, q2_sta, q3_sta, q3_sta, q1_a, q2_a, q3_a)    
    '''
    # scale starting magnet profiles if needed
    #
    # e.g. you want to start with a different effective length for magnets, or different starting positions, strengths, etc.
    #
    params_default = initial_magnet_profiles()
    params_scaled = np.array(
    [params_default[0]*param_scale['ss_off'],params_default[1]*param_scale['sl'],params_default[2]*param_scale['ss'], # solenoid start, length, strength
    params_default[3]*param_scale['qs_off'][0],params_default[4]*param_scale['ql'][0],params_default[5]*param_scale['qs'][0], # quad 1 start, length, strength
    params_default[6]*param_scale['qs_off'][1],params_default[7]*param_scale['ql'][1],params_default[8]*param_scale['qs'][1], # quad 2 start, length, strength
    params_default[9]*param_scale['qs_off'][2],params_default[10]*param_scale['ql'][2],params_default[11]*param_scale['qs'][2], # quad 3 start, length, strength
    params_default[12]*param_scale['qa'][0],params_default[13]*param_scale['qa'][1],params_default[14]*param_scale['qa'][2] # quad 1,2,3 angle
    ])
    
    # scaling profiles by 'a' vector
    #
    # e.g. the vector a is the term being 'optimized' in this case we scale magnet starting locations, strength, and angles as those are the terms being tuned.
    #
    params = np.array(
    [params_scaled[0]*a[0], params_scaled[1], params_scaled[2]*a[1],
    params_scaled[3]*a[2], params_scaled[4], params_scaled[5]*a[3],
    params_scaled[6]*a[4], params_scaled[7], params_scaled[8]*a[5],
    params_scaled[9]*a[6], params_scaled[10], params_scaled[11]*a[7],
    params_scaled[12]*a[8], params_scaled[13]*a[9], params_scaled[14]*a[10]
    ])
    
    return params
    
    
def createQuadProfile(amplitude=1.0,qlength=0.1,dlength=0.1,numQuads=2,polarity=[1,-1],stepsize=0.01):
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
    qLocations = np.cumsum(qLocations)
    qStartLocations = qLocations[0::2]
    qEndLocations = qLocations[1::2]
    
    ##########################
    print("\n")
    print(amplitude)
    print(qlength)
    print(dlength)
    print(polarity)
    print(qStartLocations)
    print(qEndLocations)
    print("\n")

    # setup the length of our quadrupole pattern
    zi = 0.0
    ze = np.sum(qlength[0:numQuads]) + np.sum(dlength[0:numQuads])
    z = np.arange(zi,ze,stepsize)
    N = len(z)
    k = np.zeros(N)
    print(zi,ze)
    
    ## iterate through z and fill in the quad values    
    ii=0 # current quad index values
    for i,zv in enumerate(z):
        if (zv >= qStartLocations[ii] and zv < qEndLocations[ii]):
            k[i] = amplitude[ii]*polarity[ii]
        elif (zv >= qEndLocations[ii]):
            ii+=1       
            
    return z,k                    
            
    
    







    
