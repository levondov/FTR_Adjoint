from moment_equations_util import *
from moment_equations import *
import matplotlib.pyplot as plt

# constants across functions
k0 = 10
h = 0.0001 # ode step size, make sure h << magnet thicknesses

#physics settings
physics_params = {'energy': 10e3, #eV
                    'current': 0.0, #Amps
                    'pipe_radius': 0.0 #meters
                    }                   
#############################################################

def getLattice(an):

    qstr1,qstr2,qstr3 = an[0],an[1],an[2]
    qlen1,qlen2,qlen3 = an[3],an[4],an[5]
    dlen1,dlen2,dlen3 = an[6],an[7],an[8]
    repeat = 2 # periodicity

    amplitude = [qstr1,qstr2,qstr3] # quadrupole amplitude
    qlength = [qlen1,qlen2,qlen3] # quadrupole length
    dlength = [dlen1,dlen2,dlen3] # drift space length
    polarity = [1,-1,1] # quadrupole polarity

    lattice = CreateLatticeProfile(amplitude,qlength,dlength,polarity,repeat,verbose=False)
    # reverse lattice for adj
    lattice_r = np.copy(lattice)
    lattice_r[:,[0,1]] = lattice_r[:,[1,0]]
    lattice_r = np.flip(lattice_r,0)
    
    return lattice,lattice_r

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
    FoM1 = 0.5*np.sum( (y[0:3,-1]-y[0:3,0])**2 )*k0**2.
    FoM2 = 0.5*np.sum( (y[3:6,-1]-y[3:6,0])**2 )
    FoM3 = 0.5*np.sum( (y[6:9,-1]-y[6:9,0])**2 )*k0**(-2.)
    FoM4 = 0.5*np.sum( (y[9,-1]-y[9,0])**2 )
        
    FoM = FoM1 + FoM2 + FoM3 + FoM4
    FoMp = np.array([FoM1,FoM2,FoM3,FoM4])  
          
    return FoM,FoMp        
        
def get_dFOM(y):
    # derivative of FoM with respect to Q,P,E,L
    # dF/dQ
    dQ_p = ( y[0,-1]-y[0,0] )*k0**2.
    dQ_m = ( y[1,-1]-y[1,0] )*k0**2.
    dQ_x = ( y[2,-1]-y[2,0] )*k0**2.
    # dF/dP
    dP_p = ( y[3,-1]-y[3,0] )
    dP_m = ( y[4,-1]-y[4,0] )
    dP_x = ( y[5,-1]-y[5,0] )
    # dF/dE
    dE_p = ( y[6,-1]-y[6,0] )*k0**(-2.)
    dE_m = ( y[7,-1]-y[7,0] )*k0**(-2.)
    dE_x = ( y[8,-1]-y[8,0] )*k0**(-2.)
    # dF/dL
    dL = np.abs( y[9,-1]-y[9,0] )
    
    return np.array([dQ_p,dQ_m,dQ_x,dP_p,dP_m,dP_x,dE_p,dE_m,dE_x,dL])

def get_dFOM_X(y_adj):
    # gradient of FOM_p with respect to X
    # dF/dX for dQ
    dQ_p = -( y_adj[6,0]-y_adj[6,-1] )*k0**(-1.)
    dQ_m = -( y_adj[7,0]-y_adj[7,-1] )*k0
    dQ_x = -( y_adj[8,0]-y_adj[8,-1] )*k0
    # dF/dX for dP
    dP_p = ( y_adj[3,0]-y_adj[3,-1] )
    dP_m = ( y_adj[4,0]-y_adj[4,-1] )
    dP_x = ( y_adj[5,0]-y_adj[5,-1] )
    # dF/dX for dE
    dE_p = -( y_adj[0,0]-y_adj[0,-1] )*k0**(1.)
    dE_m = -( y_adj[1,0]-y_adj[1,-1] )*k0**(1.)
    dE_x = -( y_adj[2,0]-y_adj[2,-1] )*k0**(1.)
    # dF/dX for dL
    dL = -( y_adj[9,0]-y_adj[9,-1] )
    
    return np.array([dQ_p,dQ_m,dQ_x,dP_p,dP_m,dP_x,dE_p,dE_m,dE_x,dL])  
    
def get_dFOM_a(z,Y,Y_adj,O_per,N_per,ACT_per):
    # gradient of FOM_p with respect to "a" parameter
    
    # two integrals to calculate with 4 pieces in each integral
    int1_1 = np.zeros(len(z))
    int1_2 = np.copy(int1_1)
    int1_3 = np.copy(int1_1)
    int1_4 = np.copy(int1_1)
    
    int2_1 = np.copy(int1_1)
    int2_2 = np.copy(int1_1)
    int2_3 = np.copy(int1_1)
    int2_4 = np.copy(int1_1)
    
    for i in range(len(z)):       
        
        # calculate integral 1
        int1_1[i] = np.dot( ACT_per[3:6,i], Y[3:6,i] )
        int1_2[i] = -ACT_per[12,i]*Y[9,i]
        int1_3[i] = -np.dot( ACT_per[9:12,i], Y[0:3,i] ) 
        int1_4[i] = -np.dot( ACT_per[0:3,i], Y[6:9,i] )
        
        # calculate integral 2
        int2_1[i] = np.dot( Y_adj[3:6,i], np.matmul( O_per[i],Y[0:3,i] ) )
        int2_2[i] = np.dot( Y[0:3,i], N_per[i] )*Y_adj[9,i]
        int2_3[i] = -np.dot( Y_adj[0:3,i], np.matmul( O_per[i],Y[3:6,i] ) )
        int2_4[i] = -np.dot( Y_adj[0:3,i], N_per[i] )*Y[9,i]
    
    #int_val = np.trapz(int1_1+int1_2+int1_3+int1_4,z) + np.trapz(int2_1+int2_2+int2_3+int2_4,z)
    int_val = np.trapz(int2_1+int2_2+int2_3+int2_4,z)
    return int_val      
    
def gd_FOM(a,X):
    # grab parameters
    init_cond = X
    lattice,_ = getLattice(a)
    
    # run moment equations
    z,y,motion,k = run_moments(init_cond, lattice, h, physics_params, verbose=True)
    return get_FOM(y)
    
def gd_dFOM(a,X,z,y,y_adj,O_nope,N_nope,ACT_nope):    
    ## get dF_p / dX , respect to initial conditions
    dFOM_X = get_dFOM_X(y_adj)
    
    ## get dF_p / da , respect to parameters    
    # init some variables
    dFOM_a = np.zeros(len(a))
    a_copy =  a.copy()
    perturb = 0.001
    import time
    for i in range(len(a_copy)):
        # perturb
        a[i] = a_copy[i] + a_copy[i]*perturb
        # grab updated lattice
        lattice,_ = getLattice(a)
        k = getLatticeKvsZ(lattice,h)
        # calculate O,N perturbed matrices 
        O,N,ACT = get_ON_and_ACT(z,y,yadj,k,physics_params)
        for j in range(len(z)):
            O[j] = O[j] - O_nope[j]
            N[j] = N[j] - N_nope[j]
            
        # calculate integral
        dFOM_a[i] = get_dFOM_a(z,y,y_adj,O,N,ACT)
        
    return dFOM_a,dFOM_X        

def gd_adj(a,X):
    # grab parameters
    init_cond = X
    lattice,lattice_r = getLattice(a)
    
    # run moment equations to get values at z=z_final
    z,y_m,motion,_ = run_moments(init_cond, lattice, h, physics_params, verbose=True)
    
    # setup adjoint initial conditions @ z=z_final
    dF = get_dFOM(y_m)
    init_cond_adj = np.array([ -dF[6],-dF[7],-dF[8],dF[3],dF[4],dF[5],-dF[0],-dF[1],-dF[2],-dF[9],0 ])
    init_cond_adj = np.concatenate(( init_cond_adj,np.reshape(y_m[:,-1].T,(11)) ))
    
    # Run adjoint equations backwards starting from z=z_final to z=0
    z,y_mom,y_adj,motion,k = run_moments_adjoint(init_cond_adj, lattice_r, -h, physics_params, verbose=True)
    
    return z,y_mom,y_adj,motion,k

#############################################################

# a and X parameters for lattice optimization                    
qstr = 0.02661
ql = 0.05425
a0 = [qstr,qstr,qstr,ql/2.0,ql,ql/2.0,0.0,0.025,0.025]
X0 = Initial_conditions()
    
# compute adjoint equations
z,y,yadj,_,k = gd_adj(a0,X0)
O0,N0,ACT0 = get_ON_and_ACT(z,y,yadj,k,physics_params)

# compute initial values for FoM
[f0,f0p] = gd_FOM(a0,X0)
f00 = f0

# compute initial values for FoM derivatives w/ respect to X and a
dfX0,dfa0 = gd_dFOM(a0,X0,z,y,yadj,O0,N0,ACT0)

plt.figure()
plt.plot(z,k)
plt.show()

























   

