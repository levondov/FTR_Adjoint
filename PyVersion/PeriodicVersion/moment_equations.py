import numpy as np
import scipy as sp
from moment_equations_util import *

# Python version of Tom's Moment equations

def run_moments(init_conditions, lattice, h, physics_params,verbose=False):
    '''
        INPUTS:
            init_conditions - initial conditions for the integration
            lattice - magnet profile (from CreateLatticeProfile()
            h - integration step size
            z_interval - integration interval
            physics_params - list of physics params
            
        OUTPUTS:
            z - z position of solved ode (independent variable)
            y - variables solved in ode
            motion - constant of motion in the system
    '''
    
    # grab some more physics params
    physics_params["rigidity"], physics_params["perveance"] = Get_beamridg_and_perv(physics_params["energy"],physics_params["current"])  
    
    odefunc = lambda z,Y,dbdx : ode_moments(z,Y,dbdx,physics_params)
    z,y = ode3(odefunc, h, init_conditions, lattice, verbose_f=verbose)

    # Constant of Motion %
    motion = get_COM(y)
    
    return z,y,motion
    
def run_moments_adjoint(init_conditions, lattice, h, physics_params,verbose=False):
    # grab some more physics params
    physics_params["rigidity"], physics_params["perveance"] = Get_beamridg_and_perv(physics_params["energy"],physics_params["current"])  
    
    odefunc = lambda z,Y,dbdx : ode_moments_adjoint(z,Y,dbdx,physics_params)
    z,y = ode3(odefunc, h, init_conditions, lattice, verbose_f=verbose)
    
    y_mom = np.flip(y[11:,:],1)
    y_adj = np.flip(y[0:11,:],1)
    
    # Constant of Motion %
    motion = get_COM(y_mom)
    
    return np.flip(z),y_mom,y_adj,motion    
    
def ode_moments(z,Y,quad_dbdx,physics_params):
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
    
    physics_params should have values for: rigidity, perveance, pipe_radius
    '''

    # quads
    psi = 0.0
    k_sol = 0.0
    k_quad = quad_dbdx / physics_params["rigidity"]
            
    # cq, sq quad calculations         
    cq = np.cos(2.0*Y[10]-2.0*psi)
    sq = np.sin(2.0*Y[10]-2.0*psi)
    
    # space charge calculation
    Q_delta = np.sqrt( Y[0]**2 - Y[1]**2 - Y[2]**2 ) if Y[0]**2 - Y[1]**2 - Y[2]**2 >= 0 else 2**(-64.0)
    ab4 = 1.0 / Q_delta; # the 4/ab term in equation
    ca_ab4 = -Y[1] / ( (Y[0]+Q_delta)*Q_delta ) # 4c_alpha/ab
    sa_ab4 = -Y[2] / ( (Y[0]+Q_delta)*Q_delta ) # 4s_alpha/ab
    
    # pipe radius calculation
    pipe_constant = (8*physics_params["perveance"]/physics_params["pipe_radius"]**4) if physics_params["pipe_radius"] != 0 else 0 # zero pipe radius leads to infinity
    
    # O and N matrix calculations
    O_mat = np.array([[-k_sol**2/2.0 + ab4*physics_params["perveance"], 2.0*k_quad*cq + ca_ab4*physics_params["perveance"], -2.0*k_quad*sq + sa_ab4*physics_params["perveance"]],
        [2.0*k_quad*cq + ca_ab4*physics_params["perveance"], -k_sol**2/2.0 + ab4*physics_params["perveance"], 0],
        [-2.0*k_quad*sq + sa_ab4*physics_params["perveance"], 0, -k_sol**2/2.0 + ab4*physics_params["perveance"]]]) + pipe_constant*np.array([ [0, -Y[1], -Y[2]],[-Y[1], 0, 0],[-Y[2], 0, 0] ])

    N_mat = np.array([[0], [2.0*k_quad*sq - sa_ab4*physics_params["perveance"]], [2.0*k_quad*cq + ca_ab4*physics_params["perveance"]]])
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
    
    return dydt
    
def ode_moments_adjoint(z,Yt,quad_dbdx,physics_params):  
    '''
    Main function to solve the adjoint equations (backwards)
    
    Here we have the original 11 moments + the 11 adjoint moments for a total of a 22 variable ode solve
    '''
    Y = Yt[0:11] # adjoint variables
    Y2 = Yt[11:] # moment variables
    
    # quads
    psi = 0.0
    k_sol = 0.0
    k_quad = quad_dbdx / physics_params["rigidity"]  
            
    # calculations            
    cq = np.cos(2.0*Y2[10]-2.0*psi)
    sq = np.sin(2.0*Y2[10]-2.0*psi)
    
    # space charge stuff
    Q_delta = np.sqrt( Y2[0]**2 - Y2[1]**2 - Y2[2]**2 ) if Y2[0]**2 - Y2[1]**2 - Y2[2]**2 >= 0 else 2**(-64.0)
    ab4 = 1.0 / Q_delta; # the 4/ab term in equation
    ca_ab4 = -Y2[1] / ( (Y2[0]+Q_delta)*Q_delta ) # 4c_alpha/ab
    sa_ab4 = -Y2[2] / ( (Y2[0]+Q_delta)*Q_delta ) # 4s_alpha/ab
    
    # pipe radius calculation
    pipe_constant = (8*physics_params["perveance"]/physics_params["pipe_radius"]**4) if physics_params["pipe_radius"] != 0 else 0 # zero pipe radius leads to infinity
        
    # O and N matrix calculations
    O_mat = np.array([[-k_sol**2/2.0 + ab4*physics_params["perveance"], 2.0*k_quad*cq + ca_ab4*physics_params["perveance"], -2.0*k_quad*sq + sa_ab4*physics_params["perveance"]],
        [2.0*k_quad*cq + ca_ab4*physics_params["perveance"], -k_sol**2/2.0 + ab4*physics_params["perveance"], 0],
        [-2.0*k_quad*sq + sa_ab4*physics_params["perveance"], 0, -k_sol**2/2.0 + ab4*physics_params["perveance"]]]) + pipe_constant*np.array([ [0, -Y2[1], -Y2[2]],[-Y2[1], 0, 0],[-Y2[2], 0, 0] ])

    N_mat = np.array([[0], [2.0*k_quad*sq - sa_ab4*physics_params["perveance"]], [2.0*k_quad*cq + ca_ab4*physics_params["perveance"]]])
    - pipe_constant*np.array([[0], [-Y2[2]], [Y2[1]]]) # pipe radius image forces addition
    
    # Calculate special matrix due to space charge variations
    [Mq,Mp,Mn] = get_SCVM(Y2,physics_params)
    
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
    
    # dE dot (y)
    dedot = np.matmul( np.reshape(Y[3:6], (1,3)), Mq )[0] + \
        Y[9]*np.matmul( np.reshape(Y2[0:3], (1,3)), Mn )[0] - \
        np.matmul( np.reshape(Y[0:3], (1,3)), Mp )[0] - \
        Y2[9]*np.matmul( np.reshape(Y[0:3], (1,3)), Mn )[0] 
    # dQ dot (y)
    dqdot = np.array([0,0,0])
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
    
    return np.concatenate((dydt2,dydt))
    
def get_COM(y):
    '''
    Calculate the constant of motion
    '''
    L = y[9,:]
    EQ = y[6,:]*y[0,:] + y[7,:]*y[1,:] + y[8,:]*y[2,:]
    PP = y[3,:]**2 + y[4,:]**2 + y[5,:]**2
    motion = EQ + (0.5)*L**2 - (0.5)*PP
    
    return motion
    
def get_FOM1(y,k0,komega,physics_params,ynorm=None):
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
    FoM4 = 0.5*(k0**(-2))*(y[6] - 0.5*(komega**2)*y[0] + physics_params["perveance"])**2
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
        
def get_dFOM1(y,k0,komega,physics_params):
    '''
    Calculate derivatives of FOM1 for initial conditions of the adjoint variables
    '''        
    k_perv = physics_params["perveance"]
    # now calculate derivatives for adjoint variables
    # adjoint variables calculated from FoM
    dP_p = y[3]
    dP_m = y[4]
    dP_x = y[5]
    
    dE_p = k0**(-2)*(y[6]-0.5*komega**2*y[0]+k_perv)*0.5*komega**2 - 2*y[6]*(2*y[6]*y[0]-y[9]**2)
    dE_m = -k0**(2)*y[1]
    dE_x = -k0**(2)*y[2]
    
    dQ_p = -k0**(-2)*(y[6]-0.5*komega**2*y[0]+k_perv) - 2*y[0]*(2*y[6]*y[0]-y[9]**2)
    dQ_m = -k0**(-2)*y[7]
    dQ_x = -k0**(-2)*y[8]
    
    dL = -2*y[9]*( 2*y[6]*y[0]-y[9]**2 )
    
    return np.array([dQ_p,dQ_m,dQ_x,dP_p,dP_m,dP_x,dE_p,dE_m,dE_x,dL])
    
    
def get_SCVM(Y,physics_params):
    '''
    find space charge variation matrices needed for adjoint calculation
    '''
    k_perv = physics_params["perveance"]
    
    Q_delta = np.sqrt( Y[0]**2 - Y[1]**2 - Y[2]**2 )
    H = k_perv*np.log(Y[0] + Q_delta)
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
    
    
    
    
    
    
    
    
    

