from moment_equations_util import *
from moment_equations import *
import matplotlib.pyplot as plt

# constants across functions
k0 = 10
h = 0.0001 # ode step size, make sure h << magnet thicknesses

repeat = 2 # periodicity
qstr = 0.02661
ql = 0.05425

amplitude = [qstr]*3 # quadrupole amplitude
qlength = [ql/2.0,ql,ql/2.0] # quadrupole length
dlength = [0.0, 0.025, 0.025] # drift space length
polarity = [1,-1,1] # quadrupole polarity

lattice = CreateLatticeProfile(amplitude,qlength,dlength,polarity,repeat,verbose=False)
# reverse lattice for adj
lattice_r = np.copy(lattice)
lattice_r[:,[0,1]] = lattice_r[:,[1,0]]
lattice_r = np.flip(lattice_r,0)

#physics settings
physics_params = {'energy': 10e3, #eV
                    'current': 0.0, #Amps
                    'pipe_radius': 0.0 #meters
                    }

#############################################################

def Initial_conditions():
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

# Pick a FoM and define functions
def get_FOM(y):
    
    # figure of merit broken into pieces for ease of reading    
    FoM1 = 0.5*np.sum( (y[0:3,-1]-y[0:3,0])**2 )
    FoM2 = 0.5*np.sum( (y[3:6,-1]-y[3:6,0])**2 )
    FoM3 = 0.5*np.sum( (y[6:9,-1]-y[6:9,0])**2 )
    FoM4 = 0.5*np.sum( (y[9,-1]-y[9,0])**2 )
        
    FoM = FoM1 + FoM2 + FoM3 + FoM4
    FoMp = np.array([FoM1,FoM2,FoM3,FoM4])  
          
    return FoM,FoMp        
        
def get_dFOM(y):
    # derivative of FoM
    # dF/dQ
    dQ_p = np.abs( y[0,-1]-y[0,0] )
    dQ_m = np.abs( y[1,-1]-y[1,0] )
    dQ_x = np.abs( y[2,-1]-y[2,0] )
    # dF/dP
    dP_p = np.abs( y[3,-1]-y[3,0] )
    dP_m = np.abs( y[4,-1]-y[4,0] )
    dP_x = np.abs( y[5,-1]-y[5,0] )
    # dF/dE
    dE_p = np.abs( y[6,-1]-y[6,0] )
    dE_m = np.abs( y[7,-1]-y[7,0] )
    dE_x = np.abs( y[8,-1]-y[8,0] )
    # dF/dL
    dL = np.abs( y[9,-1]-y[9,0] )
    
    return np.array([dQ_p,dQ_m,dQ_x,dP_p,dP_m,dP_x,dE_p,dE_m,dE_x,dL])
    
def gd_FOM(an):
    init_cond = Initial_conditions()*an
    z,y,motion = run_moments(init_cond, lattice, h, physics_params, verbose=True)
    return get_FOM(y)
    
def gd_adj(an):
    # run moment equations to get values at z=z_final
    init_cond = Initial_conditions()*an
    z,y_m,motion = run_moments(init_cond, lattice, h, physics_params, verbose=True)
    
    # setup adjoint initial conditions @ z=z_final
    dF = get_dFOM(y_m)
    init_cond_adj = np.array([ -dF[6],-dF[7],-dF[8],dF[3],dF[4],dF[5],-dF[0],-dF[1],-dF[2],dF[9],0 ])
    init_cond_adj = np.concatenate(( init_cond_adj,np.reshape(y_m[:,-1].T,(11)) ))
    
    # Run adjoint equations backwards starting from z=z_final to z=0
    z,y_mom,y_adj,motion = run_moments_adjoint(init_cond_adj, lattice_r, -h, physics_params, verbose=True)
    
    return z,y_mom,y_adj,motion

z,y,yadj,_ = gd_adj(np.ones(11))


   

