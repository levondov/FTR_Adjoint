from moment_equations import *
import matplotlib.pyplot as plt
import sys

# code for cmd line arguments, skip ahead
#############################################################################
output_file_flag = False  
if len(sys.argv) == 1:
# no inputs, need an a vector to run
    a = np.ones(11)
if len(sys.argv) == 2:
# dont save file
    a_file = sys.argv[1]
    a = np.genfromtxt(a_file)      
if len(sys.argv) == 3:
    a_file = sys.argv[1]
    a = np.genfromtxt(a_file) 
    output_file = sys.argv[2]
    output_file_flag = True    
############################################################################

# magnet field scaling parameters (all multiplicative scaling)
scale_fields={
'ql': [100.0,100.0,100.0], # quad1,2,3 length scaling
'qs': [0.01,0.01,0.01], # quad1,2,3 strength scaling
'qs_off': [1.0,1.0,1.0], # quad1,2,3 starting offset scaling
'qa': [1.0,1.0,1.0], # quad1,2,3 angle scaling
'sl': 1.0, # solenoid length scaling
'ss': 1.0, # solenoid strength scaling
'ss_off': 1.0} # solenoid starting offset scaling

# ode solver settings
h = 0.00001 # ode step size, make sure h << magnet thicknesses
z_interval=[0.0,0.313] # length of integration (meters)

# physics settings
current = 1.0e-3 # Amps
energy = 5e3 # eV
pipe_radius = 0.0 # meters

# grab magnet params scaled by vector a and scaled_fields
params = get_params(a, param_scale=scale_fields)

# run integration
z,y,motion,ksol,kquad = run_moments(params, h, z_interval, energy, current, pipe_radius, hardedge_flag=1)



# plot x^2 y^2 moments
plt.figure()
plt.plot(z,1e6*(y[0,:]+y[1,:]))
plt.plot(z,1e6*(y[0,:]-y[1,:]))
plt.plot(z,motion)
plt.plot(z,ksol/6e0,color='m')
plt.plot(z,kquad/1e3,color='k')
plt.show()

# save results
if output_file_flag:
    z=np.reshape(z,(1,len(z)))
    np.savetxt(output_file,np.r_[z,y].T[0::10,:],delimiter=',')










