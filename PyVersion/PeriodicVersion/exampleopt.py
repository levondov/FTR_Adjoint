from moment_equations import *
from moment_equations_util import *
import matplotlib.pyplot as plt
from scipy.optimize import minimize

########################## ode solver settings
h = 0.0001 # ode step size, make sure h << magnet thicknesses


##########################  create magnet profile
#number of quads
numQuads = 3

amplitude = 1.0*1e-3 # quadrupole amplitude
qlength = [0.05,0.1,0.05] # quadrupole length
dlength = [0.0, 0.2, 0.2] # drift space length
polarity = [1,-1,1] # quadrupole polarity

stepsize=h/2.0 # step size - this needs to be half the integration step size
zq,quadprofile = CreateQuadProfile(amplitude,qlength,dlength,numQuads,polarity,stepsize)

z_interval=[0.0,zq[-1]] # length of integration (meters)

########################## physics settings
physics_params = {'energy': 10e3, #eV
                    'current': 0.0, #Amps
                    'pipe_radius': 0.0 #meters
                    }


########################## Initial conditions (you can make your own)
Q_plus = 0.5*(2.2581**2*1e-6 + 0.2258**2*1e-6)
Q_minus = 0.0
Q_x = 0.0
P_plus = 0
P_minus = 0
P_x = 0
E_plus = 0.0
E_minus = 0.0
E_x = 0
init_cond = np.array([Q_plus,Q_minus,Q_x,P_plus,P_minus,P_x,E_plus,E_minus,E_x])     
         
# opt function
def OptFunc(init_values):
    _,y,_ = run_moments(np.append(init_values,[0.0,0.0]), quadprofile, h, z_interval, physics_params, verbose=False)
    FoM_match = 0.5*np.sum((y[:,-1]-y[:,0])**2.0) + y[1,0]**2
    return FoM_match
    
Nfeval=1    
def callbackF(Xi):
    global Nfeval
    print('{0:4d}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}   {4: 3.6f}   {5: 3.6f}   {6: 3.6f}   {7: 3.6f}   {8: 3.6f}   {9: 3.6f}'.format(Nfeval, Xi[0]*1e6, Xi[1]*1e6, Xi[2]*1e6, Xi[3]*1e6, Xi[4]*1e6, Xi[5]*1e6, Xi[6]*1e6, Xi[7]*1e6, Xi[8]*1e6)) 
    Nfeval+=1
    
# run opt
res = minimize(OptFunc, init_cond, method='Nelder-Mead', tol=1e-8, callback=callbackF)


# run integration with opt results
numQuads = 12
amplitude = 1.0*1e-3 # quadrupole amplitude
qlength = [0.05,0.1,0.05] # quadrupole length
dlength = [0.0, 0.2, 0.2] # drift space length
polarity = [1,-1,1] # quadrupole polarity
stepsize=h/2.0 # step size - this needs to be half the integration step size
zq,quadprofile = CreateQuadProfile(amplitude,qlength,dlength,numQuads,polarity,stepsize)
z_interval=[0.0,zq[-1]] # length of integration (meters)
z,y,opt = run_moments(np.append(res.x,[0.0,0.0]), quadprofile, h, z_interval, physics_params, verbose=False)


##############################################################################
# plot x^2 y^2 moments
plt.figure()
xrms=1e6*(y[0,:]+y[1,:])
yrms=1e6*(y[0,:]-y[1,:])
plt.plot(z,xrms,color='C0')
plt.plot(z,yrms,color='C1')

kscale = (np.max(xrms) - np.min(xrms))/np.max(quadprofile)*0.25
koffset = np.max(xrms)*0.75
plt.plot(zq,quadprofile*kscale+koffset,color='k')
plt.show()

# save results
if 0:
    z=np.reshape(z,(1,len(z)))
    np.savetxt(output_file,np.r_[z,y].T[0::10,:],delimiter=',')










