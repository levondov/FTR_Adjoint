from analysis_util import *
from gd_adj_match_1var_functions import *
               
#############################################################

# a and X parameters for lattice optimization                    
qstr = 0.02661
ql = 0.05425
a0 = [0.00638462, 0.00636695, 0.00626687] #[qstr,qstr,qstr]
X0 = Initial_conditions()
    
# compute adjoint equations
z,y,yadj,_,k = gd_adj(a0,X0)
O0,N0,ACT0 = get_ON_and_ACT(z,y,yadj,k,physics_params)

# compute initial values for FoM
[f0,f0p] = gd_FOM(a0,X0)
f00 = f0

# compute initial values for FoM derivatives w/ respect to X and a
dfa0,dfX0 = gd_dFOM(a0,X0,z,y,yadj,O0,N0,ACT0)

# scaling term for gd
gamma_X0 = (f0/np.sum(dfX0**2))

# setup variables for holding history of iterations
gamma_h = np.reshape(gamma_X0,(1,1))
X_h = np.reshape(X0.copy(),(1,len(X0)))
f_h = np.reshape(f0,(1,1))
dfX_h = np.reshape(dfX0.copy(),(1,len(dfX0)))

# adjust starting gamma
# iterate a
X_t = X0 - gamma_X0*dfX0
X_h = np.concatenate((X_h,np.reshape(X_t,(1,len(X_t)))),axis=0)

# compute initial values for FoM
[f_t,_] = gd_FOM(a0,X_h[-1,:])
f_h = np.concatenate((f_h,np.array(f_t).reshape((1,1))),axis=0)
PrintFoM(f_t)

# iterate to find ideal starting gamma
while f_h[-1,0] >= f_h[0,0]:
    gamma_h = np.concatenate((gamma_h,np.array(gamma_h[-1,:]/2.0).reshape((1,1))),axis=0)

    X_t = X0 - gamma_h[-1,0]*dfX0
    X_h = np.concatenate((X_h,np.reshape(X_t,(1,len(X_t)))),axis=0)

    [f_t,_] = gd_FOM(a0,X_h[-1,:])
    f_h = np.concatenate((f_h,np.array(f_t).reshape((1,1))),axis=0)
    PrintFoM(f_t)

# start main optimization loop
while True:
    ii=1
    while f_h[-1,0] < f_h[-2,0]:
        print('Iterating '+str(ii))

        # iterate a
        X_t = X_h[-1,:] - gamma_h[-1,0]*dfX_h[-1,:]
        X_h = np.concatenate((X_h,np.reshape(X_t,(1,len(X_t)))),axis=0)

        # compute FoM
        [f_t,_] = gd_FOM(a0,X_h[-1,:])
        f_h = np.concatenate((f_h,np.array(f_t).reshape((1,1))),axis=0)
        PrintFoM(f_t)

        # if we have iterated too much, start increasing gamma
        ii += 1
        if (ii > 20):
            gamma_h = np.concatenate((gamma_h,np.array(gamma_h[-1,:]*2.0).reshape((1,1))),axis=0)

    # recompute adjoint equations for new search direction
    print("Recomputing adjoint equations")

    # grab last good settings
    X_h = np.concatenate((X_h,np.reshape(X_h[-2,:],(1,len(X_h[-2,:])))),axis=0)
    f_h = np.concatenate((f_h,np.array(f_h[-2,0]).reshape((1,1))),axis=0)

    # calculate adjoint equations
    z,y,yadj,_,k = gd_adj(a0,X_h[-1,:])
    O0,N0,ACT0 = get_ON_and_ACT(z,y,yadj,k,physics_params)

    # calculate new grdient
    _,dfX_t = gd_dFOM(a0,X_h[-1,:],z,y,yadj,O0,N0,ACT0)
    dfX_h = np.concatenate((dfX_h,np.reshape(dfX_t,(1,len(dfX_t)))),axis=0)

    if (ii == 2): # meaning no improvement from recalculating gradient
        # change gamma
        print('Updating gamma')

        f0n = f_h[-1,0]
        X0n = X_h[-1,:]
        while f_h[-1,0] >= f0n:
            gamma_h = np.concatenate((gamma_h,np.array(gamma_h[-1,:]/2.0).reshape((1,1))),axis=0)

            # iterate a
            X_t = X0n - gamma_h[-1,0]*dfX_h[-1,:]
            X_h = np.concatenate((X_h,np.reshape(X_t,(1,len(X_t)))),axis=0)

            # compute FoM
            [f_t,_] = gd_FOM(a0,X_h[-1,:])
            f_h = np.concatenate((f_h,np.array(f_t).reshape((1,1))),axis=0)
            PrintFoM(f_t)

    if f_h[-1,0] < 5e-18:
        break

    if len(f_h[:,0]) > 5000:
        break

SaveData("testx2",[a0,X_h,f_h,dfX_h,gamma_h])






















   

