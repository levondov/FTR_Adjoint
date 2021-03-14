from moment_equations import *
import matplotlib.pyplot as plt

a = np.array([
   1.129928869681529,
   1.042878646468118,
   0.996927320742096,
   1.040384838649907,
   1.070238743182279,
   1.075813979350096,
   1.026419246981094,
   1.047129360787238,
   1.054061572281601,
   1.035232396154389,
   1.031077426545717])

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
z = np.arange(z_interval[0],z_interval[1],h) # all steps

# physics settings
current = 1e-3 # Amps
energy = 5e3 # eV

# grab magnet params scaled by vector a and scaled_fields
params = get_params(a, param_scale=scale_fields)

radi = np.linspace(0.00225,0.01,10)
currents = np.linspace(1.0e-3,5.0e-3,5)
dt_runs = np.empty((len(radi),len(currents)),dtype=object)
dt_fmoments = np.empty((len(radi),2,5))
dt_baseruns = np.empty((len(currents)),dtype=object)

if 0:
    for j,sc in enumerate(currents):
        z,y_base,motion,ksol,kquad = run_moments(params, h, z_interval, energy, sc, 0.0, hardedge_flag=1)
        dt_baseruns[j]  = y_base
        for i,pr in enumerate(radi):
            z,y2,motion,ksol,kquad = run_moments(params, h, z_interval, energy, sc, pr, hardedge_flag=1)
            dt_runs[i,j] = y2
            dt_fmoments[i,0,j] = (y2[0,-1]+y2[1,-1])
            dt_fmoments[i,1,j] = (y2[0,-1]-y2[1,-1])

    np.save('Results_1to5mA_225to10mm',np.array([dt_baseruns,dt_runs,dt_fmoments,radi,currents],dtype=object))


dt = np.load('Results_1to5mA_225to10mm.npy',allow_pickle=True)
radi = dt[3]
currents = dt[4]
dt_fmoments=dt[2]
dt_baseruns=dt[0]
dt_runs = dt[1]

plt.figure()
plt.subplot(1,2,1)
lpos = [3.1,3.75,4.45,5.2,6.1]
for i in range(5):
    plt.plot(radi*1e3,dt_fmoments[:,0,i],marker='o',label=str(currents[i]*1e3) +' mA') 
    plt.text(7.1,lpos[i]*1e-6, str(currents[i]*1e3) +' mA')
plt.ylabel('Moments [$m^2$]')
plt.gca().ticklabel_format(axis='y',style='sci',scilimits=[0,0])
plt.title('$\langle x^2 \\rangle$')
plt.xlabel('Pipe radius [mm]')
plt.grid(True)
plt.legend()

plt.subplot(1,2,2)
for i in range(5):
    plt.plot(radi*1e3,dt_fmoments[:,1,i],marker='o',label=str(currents[i]*1e3) +' mA')    
plt.ylabel('Moments [$m^2$]')
plt.gca().ticklabel_format(axis='y',style='sci',scilimits=[0,0])
plt.title('$\langle y^2 \\rangle$')
plt.xlabel('Pipe radius [mm]')   
plt.grid(True) 



###
# plot x^2 y^2 moments
plt.figure()
plt.plot(z,(dt_baseruns[0][0,:]+dt_baseruns[0][1,:]),color='k')
plt.plot(z,(dt_baseruns[0][0,:]-dt_baseruns[0][1,:]),color='k')


plt.plot(z,(dt_runs[0,0][0,:]+dt_runs[0,0][1,:]),color='C0',label='r=2.25 mm')
plt.plot(z,(dt_runs[0,0][0,:]-dt_runs[0,0][1,:]),color='C0')

plt.plot(z,(dt_runs[1,0][0,:]+dt_runs[1,0][1,:]),color='C1',label='r=3.11 mm')
plt.plot(z,(dt_runs[1,0][0,:]-dt_runs[1,0][1,:]),color='C1')

plt.plot(z,(dt_runs[2,0][0,:]+dt_runs[2,0][1,:]),color='C2',label='r=4.00 mm')
plt.plot(z,(dt_runs[2,0][0,:]-dt_runs[2,0][1,:]),color='C2')

plt.plot(z,(dt_runs[3,0][0,:]+dt_runs[3,0][1,:]),color='C3',label='r=4.83 mm')
plt.plot(z,(dt_runs[3,0][0,:]-dt_runs[3,0][1,:]),color='C3')

plt.plot(z,(dt_runs[4,0][0,:]+dt_runs[4,0][1,:]),color='C4',label='r=5.70 mm')
plt.plot(z,(dt_runs[4,0][0,:]-dt_runs[4,0][1,:]),color='C4')

plt.ylabel('Moments [$m^2$]')
plt.gca().ticklabel_format(axis='y',style='sci',scilimits=[0,0])
plt.xlabel('Z position [m]')
plt.legend()
plt.grid(True)
plt.title('Moment runs for 1 mA beam current and different pipe radi')


plt.show()








