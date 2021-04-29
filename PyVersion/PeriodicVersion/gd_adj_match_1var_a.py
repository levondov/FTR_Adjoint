from analysis_util import *
from gd_adj_match_1var_functions import *
               
#############################################################

# a and X parameters for lattice optimization                    
qstr = 0.02661
ql = 0.05425
a0 = [qstr,qstr,qstr]
X0 = Initial_conditions()
    
# compute adjoint equations
z,y,yadj,_,k = gd_adj(a0,X0)
O0,N0,ACT0 = get_ON_and_ACT(z,y,yadj,k,physics_params)

# compute initial values for FoM
[f0,f0p] = gd_FOM(a0,X0)
f00 = f0

# compute initial values for FoM derivatives w/ respect to X and a
dfa0,_ = gd_dFOM(a0,X0,z,y,yadj,O0,N0,ACT0)

# scaling term for gd
gamma_a0 = 1e-3*(f0/np.sum(dfa0**2))

# setup variables for holding history of iterations
gamma_h = np.reshape(gamma_a0,(1,1))
a_h = np.reshape(a0.copy(),(1,len(a0)))
f_h = np.reshape(f0,(1,1))
dfa_h = np.reshape(dfa0.copy(),(1,len(dfa0)))

# adjust starting gamma
# iterate a
a_t = a0 - gamma_a0*dfa0
a_h = np.concatenate((a_h,np.reshape(a_t,(1,len(a_t)))),axis=0)

# compute initial values for FoM
[f_t,_] = gd_FOM(a_h[-1,:],X0)
f_h = np.concatenate((f_h,np.array(f_t).reshape((1,1))),axis=0)
PrintFoM(f_t)

# iterate to find ideal starting gamma
while f_h[-1,0] >= f_h[0,0]:
    gamma_h = np.concatenate((gamma_h,np.array(gamma_h[-1,:]/2.0).reshape((1,1))),axis=0)

    a_t = a0 - gamma_h[-1,0]*dfa0
    a_h = np.concatenate((a_h,np.reshape(a_t,(1,len(a_t)))),axis=0)

    [f_t,_] = gd_FOM(a_h[-1,:],X0)
    f_h = np.concatenate((f_h,np.array(f_t).reshape((1,1))),axis=0)
    PrintFoM(f_t)

# start main optimization loop
while True:
    ii=1
    while f_h[-1,0] < f_h[-2,0]:
        print('Iterating '+str(ii))

        # iterate a
        a_t = a_h[-1,:] - gamma_h[-1,0]*dfa_h[-1,:]
        a_h = np.concatenate((a_h,np.reshape(a_t,(1,len(a_t)))),axis=0)

        # compute FoM
        [f_t,_] = gd_FOM(a_h[-1,:],X0)
        f_h = np.concatenate((f_h,np.array(f_t).reshape((1,1))),axis=0)
        PrintFoM(f_t)

        # if we have iterated too much, start increasing gamma
        ii += 1
        if (ii > 20):
            gamma_h = np.concatenate((gamma_h,np.array(gamma_h[-1,:]*2.0).reshape((1,1))),axis=0)

    # recompute adjoint equations for new search direction
    print("Recomputing adjoint equations")

    # grab last good settings
    a_h = np.concatenate((a_h,np.reshape(a_h[-2,:],(1,len(a_h[-2,:])))),axis=0)
    f_h = np.concatenate((f_h,np.array(f_h[-2,0]).reshape((1,1))),axis=0)

    # calculate adjoint equations
    z,y,yadj,_,k = gd_adj(a_h[-1,:],X0)
    O0,N0,ACT0 = get_ON_and_ACT(z,y,yadj,k,physics_params)

    # calculate new grdient
    dfa_t,_ = gd_dFOM(a_h[-1,:],X0,z,y,yadj,O0,N0,ACT0)
    dfa_h = np.concatenate((dfa_h,np.reshape(dfa_t,(1,len(dfa_t)))),axis=0)

    if (ii == 2): # meaning no improvement from recalculating gradient
        # change gamma
        print('Updating gamma')

        f0n = f_h[-1,0]
        a0n = a_h[-1,:]
        while f_h[-1,0] >= f0n:
            gamma_h = np.concatenate((gamma_h,np.array(gamma_h[-1,:]/2.0).reshape((1,1))),axis=0)

            # iterate a
            a_t = a0n - gamma_h[-1,0]*dfa_h[-1,:]
            a_h = np.concatenate((a_h,np.reshape(a_t,(1,len(a_t)))),axis=0)

            # compute FoM
            [f_t,_] = gd_FOM(a_h[-1,:],X0)
            f_h = np.concatenate((f_h,np.array(f_t).reshape((1,1))),axis=0)
            PrintFoM(f_t)

    if f_h[-1,0] < 1e-16:
        break

    if len(f_h[:,0]) > 5000:
        break

SaveData("tests",[a_h,f_h,dfa_h,gamma_h])






















   

