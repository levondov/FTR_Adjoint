from moment_equations import *
from moment_equations_util import *
import matplotlib.pyplot as plt

##########################  create magnet profile
repeat = 2 # periodicity

qstr = 0.02661
ql = 0.05425

amplitude = [qstr]*3 # quadrupole amplitude
qlength = [ql/2.0,ql,ql/2.0] # quadrupole length
dlength = [0.0, 0.025, 0.025] # drift space length
polarity = [1,-1,1] # quadrupole polarity

lattice = CreateLatticeProfile(amplitude,qlength,dlength,polarity,repeat,verbose=False)

########################## physics settings
physics_params = {'energy': 10e3, #eV
                    'current': 0.0, #Amps
                    'pipe_radius': 0.0 #meters
                    }


########################## Initial conditions (you can make your own)
init_cond = Initial_conditions()     
         

# run integration
h = 0.0001 # ode step size, make sure h << magnet thicknesses
z,y,motion = run_moments(init_cond, lattice, h, physics_params, verbose=True)


##############################################################################
# plot x^2 y^2 moments

plt.figure()

# plot moments
xrms=(y[0,:]+y[1,:])
yrms=(y[0,:]-y[1,:])
plt.plot(z,xrms,color='C0',label='$\langle x^2 \\rangle$')
plt.plot(z,yrms,color='C1',label='$\langle y^2 \\rangle$')

# plot quadrupole profile for reference
kscale = (np.max(xrms) - np.min(xrms))/np.max(lattice[:,-1])*0.25
koffset = np.max(xrms)*0.75
lattice_scaled = np.copy(lattice)
lattice_scaled[:,-1] = lattice_scaled[:,-1]*kscale+koffset
plt.plot([],color='k',label='$K_q$')
PlotLatticeProfile(lattice_scaled)

plt.legend()
plt.grid(True)
plt.ylabel('moments $[m^2]$')
plt.xlabel('Z position [m]')
plt.show()

# save results
if 0:
    z=np.reshape(z,(1,len(z)))
    np.savetxt(output_file,np.r_[z,y].T[0::10,:],delimiter=',')










