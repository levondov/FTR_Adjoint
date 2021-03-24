from moment_equations import *
from moment_equations_util import *
import matplotlib.pyplot as plt

########################## ode solver settings
h = 0.0001 # ode step size, make sure h << magnet thicknesses


##########################  create magnet profile
#number of quads
numQuads = 12

amplitude = 1.0*1e-3 # quadrupole amplitude
qlength = [0.05,0.05] # quadrupole length
dlength = [0.0, 0.2, 0.0, 0.2] # drift space length
polarity = [1,-1,-1,1] # quadrupole polarity

stepsize=h/2.0 # step size - this needs to be half the integration step size
zq,quadprofile = CreateQuadProfile(amplitude,qlength,dlength,numQuads,polarity,stepsize)

z_interval=[0.0,zq[-1]] # length of integration (meters)

########################## physics settings
physics_params = {'energy': 10e3, #eV
                    'current': 0.0, #Amps
                    'pipe_radius': 0.0 #meters
                    }


########################## Initial conditions (you can make your own)
init_cond = Initial_conditions()     
         

# run integration
z,y,motion = run_moments(init_cond, quadprofile, h, z_interval, physics_params, verbose=True)


##############################################################################
# plot x^2 y^2 moments
plt.figure()
xrms=1e6*(y[0,:]+y[1,:])
yrms=1e6*(y[0,:]-y[1,:])
plt.plot(z,xrms,color='C0',label='$\langle x^2 \\rangle$')
plt.plot(z,yrms,color='C1',label='$\langle y^2 \\rangle$')

kscale = (np.max(xrms) - np.min(xrms))/np.max(quadprofile)*0.25
koffset = np.max(xrms)*0.75

plt.plot(zq,quadprofile*kscale+koffset,color='k',label='$K_q$')
plt.legend()
plt.grid(True)
plt.show()

# save results
if 0:
    z=np.reshape(z,(1,len(z)))
    np.savetxt(output_file,np.r_[z,y].T[0::10,:],delimiter=',')










