import numpy as np
import scipy as sp
from moment_equations_util import *

# Python version of Tom's Moment equations

def run_moments(params, h, z_interval, energy, current=0.0, pipe_radius=0.0, hardedge_flag=1):
    '''
    Given a set of params, run the moment equations
    params are: 
    [solenoid start position, 
    solenoid length, 
    solenoid strength (T/rad), 
    quad 1 start position, 
    quad 1 length,
    quad 1 strength (T/m),
    quad 2 start position,
    quad 2 length,
    quad 2 strength,
    quad 3 start position,
    quad 3 length,
    quad 3 strength,
    quad 1 angle (rad),
    quad 2 angle,
    quad 3 angle]
    '''
    
    # grab initial conditions for QPEL
    init_conditions = initial_conditions()
    
    # grab some other params
    rho, k_perv = get_beamridg_and_perv(energy,current)    

    # run forward integration of moment equations
    z = np.arange(z_interval[0],z_interval[1],h) # all steps
    odefunc = lambda z,Y : ode_moments(z,Y,params,rho,k_perv,pipe_radius,hardedge_flag)
    y,ksol,kquad = ode3(odefunc,z_interval[0], h, z_interval[1], init_conditions,verbose=True)

    # Constant of Motion %
    motion = get_COM(y)
    
    return z,y,motion,ksol,kquad
    
def ode_moments(z,Y,params,rho,k_perv,r_pipe,hardedge_flag=1):
    '''
    Main function that solves Tom's moment equations
    
    # params = [solenoid start, solenoid length, solenoid strength,q1 start, q1 length, q1
    # strength, q2 start, q2 length, q2 strength, q3 start, q3 length, q3 strength]
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
    
    # solenoid
    k_sol = 0.0;
    k_solvv = params[2]/rho # just need max solenoid value for GD optimizations
    if z >= params[0] and z <= (params[0]+params[1]): # if we are inside solenoid
        k_sol = params[2]/rho # Ks = Bs / rho = 0.041 (T) / 0.0667
        if not hardedge_flag:
            # cos^2 profile, not yet implimented
            pass  

    # quads
    psi = 0
    k_quad = 0.0
    if z >= params[3] and z <= (params[3]+params[4]): # if inside quad 1
        k_quad = params[5]/rho
        psi = params[12]
        if not hardedge_flag:
            # not yet implemented
            pass
    if z >= params[6] and z <= (params[6]+params[7]): # if inside quad 2
        k_quad = params[8]/rho
        psi = params[13]
        if not hardedge_flag:
            # not yet implemented
            pass
    if z >= params[9] and z <= (params[9]+params[10]): # if inside quad 3
        k_quad = params[11]/rho
        psi = params[14]
        if not hardedge_flag:
            # not yet implemented
            pass
            
    # calculations            
    cq = np.cos(2.0*Y[10]-2.0*psi)
    sq = np.sin(2.0*Y[10]-2.0*psi)
    
    # space charge stuff
    Q_delta = np.sqrt( Y[0]**2 - Y[1]**2 - Y[2]**2 )
    ab4 = 1.0 / Q_delta; # the 4/ab term in equation
    ca_ab4 = -Y[1] / ( (Y[0]+Q_delta)*Q_delta ) # 4c_alpha/ab
    sa_ab4 = -Y[2] / ( (Y[0]+Q_delta)*Q_delta ) # 4s_alpha/ab
    
    # pipe radius calculation
    pipe_constant = (8*k_perv/r_pipe**4) if r_pipe != 0 else 0 # zero pipe radius leads to infinity
    
    # Calculate O and N matrix stuff
    O_mat = np.array([[-k_sol**2/2.0 + ab4*k_perv, 2.0*k_quad*cq + ca_ab4*k_perv, -2.0*k_quad*sq + sa_ab4*k_perv],
        [2.0*k_quad*cq + ca_ab4*k_perv, -k_sol**2/2.0 + ab4*k_perv, 0],
        [-2.0*k_quad*sq + sa_ab4*k_perv, 0, -k_sol**2/2.0 + ab4*k_perv]]) + pipe_constant*np.array([ [0, -Y[1], -Y[2]],[-Y[1], 0, 0],[-Y[2], 0, 0] ])

    N_mat = np.array([[0], [2.0*k_quad*sq - sa_ab4*k_perv], [2.0*k_quad*cq + ca_ab4*k_perv]])
    - pipe_constant*np.array([[0], [-Y[2]], [Y[1]]]) # pipe radius image forces addition
    
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
    
    return dydt,k_sol,k_quad
    
def ode_moments_adjoint(z,Yt,params,rho,k_perv,hardedge_flag=1):  
    '''
    Main function to solve the adjoint equations (backwards)
    
    Here we have the original 11 moments + the 11 adjoint moments for a total of a 22 variable ode solve
    '''
    Y = Yt[0:11] # adjoint variables
    Y2 = Yt[11:] # moment variables
    
    # solenoid
    k_sol = 0.0;
    k_solvv = params[2]/rho # just need max solenoid value for GD optimizations
    if z >= params[0] and z <= (params[0]+params[1]): # if we are inside solenoid
        k_sol = params[2]/rho # Ks = Bs / rho = 0.041 (T) / 0.0667
        if not hardedge_flag:
            # cos^2 profile, not yet implimented
            pass  

    # quads
    psi = 0
    k_quad = 0.0
    if z >= params[3] and z <= (params[3]+params[4]): # if inside quad 1
        k_quad = params[5]/rho
        psi = params[12]
        if not hardedge_flag:
            # not yet implemented
            pass
    if z >= params[6] and z <= (params[6]+params[7]): # if inside quad 2
        k_quad = params[8]/rho
        psi = params[13]
        if not hardedge_flag:
            # not yet implemented
            pass
    if z >= params[9] and z <= (params[9]+params[10]): # if inside quad 3
        k_quad = params[11]/rho
        psi = params[14]
        if not hardedge_flag:
            # not yet implemented
            pass    
            
    # calculations            
    cq = np.cos(2.0*Y2[10]-2.0*psi)
    sq = np.sin(2.0*Y2[10]-2.0*psi)
    
    # space charge stuff
    Q_delta = np.sqrt( Y2[0]**2 - Y2[1]**2 - Y2[2]**2 )
    ab4 = 1.0 / Q_delta; # the 4/ab term in equation
    ca_ab4 = -Y2[1] / ( (Y2[0]+Q_delta)*Q_delta ) # 4c_alpha/ab
    sa_ab4 = -Y2[2] / ( (Y2[0]+Q_delta)*Q_delta ) # 4s_alpha/ab
    
    # Calculate O and N matrix stuff
    O_mat = np.array([[-k_sol**2/2.0 + ab4*k_perv, 2.0*k_quad*cq + ca_ab4*k_perv, -2.0*k_quad*sq + sa_ab4*k_perv],
        [2.0*k_quad*cq + ca_ab4*k_perv, -k_sol**2/2.0 + ab4*k_perv, 0],
        [-2.0*k_quad*sq + sa_ab4*k_perv, 0, -k_sol**2/2.0 + ab4*k_perv]])

    N_mat = np.array([[0], [2.0*k_quad*sq - sa_ab4*k_perv], [2.0*k_quad*cq + ca_ab4*k_perv]]) 
    
    # Calculate special matrix due to space charge variations
    [Mq,Mp,Mn] = get_SCVM(Y2,k_perv);            
    
def get_COM(y):
    '''
    Calculate the constant of motion
    '''
    L = y[9,:]
    EQ = y[6,:]*y[0,:] + y[7,:]*y[1,:] + y[8,:]*y[2,:]
    PP = y[3,:]**2 + y[4,:]**2 + y[5,:]**2
    motion = EQ + (0.5)*L**2 - (0.5)*PP
    
    return motion
    
def get_FOM(y,k0,komega,k_perv,ynorm=None):
    '''
    calculates the Figure of Merit + adjoint equation initial conditions given a set of moment values
    
    i.e.
    Given Q+,Q-,Qx,P+,P-,Px,E+,E-,Ex,L , calculate FOM:
    
    ynorm is if you want to normalize the FOM by something
    '''        
    
    # figure of merit broken into pieces for ease of reading
    FoM1 = 0.5*np.sum(y[3:6]**2)
    FoM2 = 0.5*(k0**2)*(y[1]**2 + y[2]**2)
    FoM3 = 0.5*(k0**(-2))*(y[7]**2+y[8]**2)
    FoM4 = 0.5*(k0**(-2))*(y[6] - 0.5*(komega**2)*y[0] + k_perv)**2
    FoM5 = 0.5*(2*y[6]*y[0] - y[9]**2)**2
    
    if ynorm is not None:
        # normalize
        FoM1 = FoM1 / ((k0**2) * (ynorm[0]**2))
        FoM2 = FoM2 / ((k0**2) * (ynorm[0]**2))
        FoM3 = FoM3 / ((k0**2) * (ynorm[0]**2))
        FoM4 = FoM4 / ((k0**2) * (ynorm[0]**2))
        FoM5 = FoM5 / ((k0**2) * (ynorm[0]**2))
        
    FoM = FoM1 + FoM2 + FoM3 + FoM4 + FoM5
    FoMp = np.array([FoM1,FoM2,FoM3,FoM4,FoM5])  
          
    return FoM,FoMp        
        
def get_dFOM(y,k0,komega,k_perv):
    '''
    Calculate derivatives of FOM for initial conditions of the adjoint variables
    '''        
    # now calculate derivatives for adjoint variables
    # adjoint variables calculated from FoM
    dP_p = y[3]
    dP_m = y[4]
    dP_x = y[5]
    
    dE_p = k0**(-2)*(y[6]-0.5*komega**2*y[0]+k_perv)*0.5*komega**2 - 2*y[6]*(2*[6]*y[0]-y[9]**2)
    dE_m = -k0**(2)*y[1]
    dE_x = -k0**(2)*y[2]
    
    dQ_p = -k0**(-2)*(y[6]-0.5*komega**2*y[0]+k_perv) - 2*y[0]*(2*y[6]*y[0]-y[9]**2)
    dQ_m = -k0**(-2)*y[7]
    dQ_x = -k0**(-2)*y[8]
    
    dL = -2*y[9]*(2*y[6]*y[0]-y[9]**2)
    
    return np.array([dQ_p,dQ_m,dQ_x,dP_p,dP_m,dP_x,dE_p,dE_m,dE_x,dL])
    
    
def get_SCVM(Y,k_perv):
    '''
    find space charge variation matrices needed for adjoint calculation
    '''

    Q_delta = np.sqrt( Y[0]**2 - Y[1]**2 - Y[2]**2 )
    H = k_perv*np.log(Y[0] + Q_delta)
    Q_deltaplus = Q_delta + Y[0]

    V1 = -(k_perv/(Q_delta^2))*[Y(4)-Y(2)*Y(5)/Q_deltaplus-Y(3)*Y(6)/Q_deltaplus,
        -Y(2)*Y(4)/Q_deltaplus+Y(5),
        -Y(3)*Y(4)/Q_deltaplus+Y(6)];
    U1t = (1/Q_delta)*[Y(1), -Y(2), -Y(3)];

    V2 = (k_perv/(Q_delta*(Q_deltaplus^2)))*[Y(2)*Y(5)+Y(3)*Y(6),Y(2)*Y(4),Y(3)*Y(4)]#';
    U2t = U1t + [1,0,0];

    V3 = -(k_perv/(Q_delta*Q_deltaplus))*[Y(5), Y(4), 0]#';
    V4 = -(k_perv/(Q_delta*Q_deltaplus))*[Y(6), 0, Y(4)]#';
    U3t = [0, 1, 0];
    U4t = [0,0,1];

    Mp = V1*U1t + V2*U2t + V3*U3t + V4*U4t;
    ####

    W1 = -(k_perv/Q_delta)*[1, Y(2)/Q_deltaplus, Y(3)/Q_deltaplus]#';
    W2 = (k_perv/(Q_delta*(Q_deltaplus^2)))*[Y(2)^2+Y(3)^2,Y(2)*Y(1),Y(3)*Y(1)]#';
    W3 = -(k_perv/(Q_delta*Q_deltaplus))*[Y(2),Y(1),0]#';
    W4 = -(k_perv/(Q_delta*Q_deltaplus))*[Y(3),0,Y(1)]#';

    Mq = W1*U1t + W2*U2t + W3*U3t + W4*U4t;
    ####

    X1 = -k_perv*[0, Y(3)/(Q_deltaplus*(Q_delta^2)), -Y(2)/(Q_deltaplus*(Q_delta^2))]#';
    X2 = -k_perv*[0, Y(3)/(Q_delta*(Q_deltaplus^2)), -Y(2)/(Q_delta*(Q_deltaplus^2))]#';
    X3 = k_perv*[0,0,-1/(Q_delta*Q_deltaplus)]#';
    X4 = k_perv*[0,1/(Q_delta*Q_deltaplus),0]#';

    Mn = X1*U1t + X2*U2t + X3*U3t + X4*U4t;    
    
    
    
    
    
    
    
    
    

