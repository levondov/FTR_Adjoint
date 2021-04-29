from analysis_util import *
from gd_adj_match_2var_functions import *
               
#############################################################

# a and X and Y parameters for lattice optimization      
lam = 0.1              
qstr = 0.02661
ql = 0.05425
a0 = [0.00638462, 0.00636695, 0.00626687] #[qstr,qstr,qstr]
X0 = Initial_conditions()
Y0 = Initial_conditions()
    
####### compute adjoint equations for X
z_p,y_p,yadj_p,_,k_p = gd_adj_p(a0,X0)
O0_p,N0_p,ACT0_p = get_ON_and_ACT(z_p,y_p,yadj_p,k_p,physics_params)
# compute initial values for FoM
[f0_p,f00_p] = gd_FOMp(a0,X0)
# compute initial values for FoM derivatives w/ respect to X and a
dfa0_p,dfX0_p = gd_dFOMp(a0,X0,z_p,y_p,yadj_p,O0_p,N0_p,ACT0_p)
# scaling term for gd
gamma_X0 = (f0_p/np.sum(dfX0_p**2))
####### compute adjoint equations info for Y & Fe
z_e,y_e,yadj_e,_,k_e = gd_adj_e(a0,Y0)
O0_e,N0_e,ACT0_e = get_ON_and_ACT(z_e,y_e,yadj_e,k_e,physics_params)
# compute initial values for FoM
f0_e = gd_FOMe(a0,Y0)
# compute initial values for FoM derivatives w/ respect to X and a
dfa0_e,dfY0_e = gd_dFOMe(a0,Y0,z_e,y_e,yadj_e,O0_e,N0_e,ACT0_e)
# since starting values are the same
dfY0_p = dfX0_p

# setup variables for holding history of iterations for X
gammaX_h = np.reshape(gamma_X0,(1,1))
X_h = np.reshape(X0.copy(),(1,len(X0)))
fX_p_h = np.reshape(f0_p,(1,1))
fX_e_h = np.reshape(f0_e,(1,1))
dfX_p_h = np.reshape(dfX0_p.copy(),(1,len(dfX0_p)))
# setup variables for holding history of iterations for Y (initially Y=X)
gammaY_h = np.reshape(gamma_X0,(1,1))
Y_h = np.reshape(X0.copy(),(1,len(X0)))
fY_p_h = np.reshape(f0_p,(1,1))
fY_e_h = np.reshape(f0_e,(1,1))
fY_h = np.reshape(f0_p + f0_e*lam,(1,1))
dfY_p_h = np.reshape(dfY0_p.copy(),(1,len(dfY0_p)))
dfY_e_h = np.reshape(dfY0_e.copy(),(1,len(dfY0_e)))
# setup variablkes for holding history of iterations for a
# setup variables for holding history of iterations
gammaA_h = np.reshape(gamma_X0,(1,1))
a_h = np.reshape(a0.copy(),(1,len(a0)))
dfaX_p_h = np.reshape(dfa0_p.copy(),(1,len(dfa0_p)))
dfaY_p_h = np.reshape(dfa0_p.copy(),(1,len(dfa0_p)))
fa_h = np.reshape(0.0,(1,1))

# adjust starting gamma
# iterate X
X_t = X0 - gamma_X0*dfX0_p
X_h = np.concatenate((X_h,np.reshape(X_t,(1,len(X_t)))),axis=0)
# iterate Y
Y_t = Y0 - gamma_X0*(dfY0_p + lam*dfY0_e)
Y_h = np.concatenate((Y_h,np.reshape(Y_t,(1,len(Y_t)))),axis=0)
# iterate a
a_t = a0 - (1/lam)*gamma_X0*(dfaY_p_h[-1,:] - dfaX_p_h[-1,:])
a_h = np.concatenate((a_h,np.reshape(a_t,(1,len(a_t)))),axis=0)

# compute initial values for FoM, Fp(X) , Fp(Y) , Fe(Y)
[f_t,_] = gd_FOMp(a_h[-1,:],X_h[-1,:])
fX_p_h = np.concatenate((fX_p_h,np.array(f_t).reshape((1,1))),axis=0)
PrintFoM(f_t,'Xp')
[f_t,_] = gd_FOMp(a_h[-1,:],Y_h[-1,:])
fY_p_h = np.concatenate((fY_p_h,np.array(f_t).reshape((1,1))),axis=0)
PrintFoM(f_t,'Yp')
f_t = gd_FOMe(a_h[-1,:],Y_h[-1,:])
fY_e_h = np.concatenate((fY_e_h,np.array(f_t).reshape((1,1))), axis=0)
PrintFoM(f_t,'Ye')
fY_h = np.concatenate((fY_h,np.array(fY_p_h[-1,0]+lam*fY_e_h[-1,0]).reshape((1,1))), axis=0)
fa_h = np.concatenate((fa_h,np.array((1/lam)*(fY_p_h[-1,0]-fX_p_h[-1,0])).reshape((1,1))), axis=0)

# iterate to find ideal starting gamma for X
while fX_p_h[-1,0] >= fX_p_h[0,0]:
    gammaX_h = np.concatenate((gammaX_h,np.array(gammaX_h[-1,:]/2.0).reshape((1,1))),axis=0)

    X_t = X0 - gammaX_h[-1,0]*dfX0_p
    X_h = np.concatenate((X_h,np.reshape(X_t,(1,len(X_t)))),axis=0)

    [f_t,_] = gd_FOMp(a_h[-1,:],X_h[-1,:])
    fX_p_h = np.concatenate((fX_p_h,np.array(f_t).reshape((1,1))),axis=0)
    PrintFoM(f_t,'Xp')

# iterate to find ideal starting gamma for Y
# use same gamma as X for now
gammaY_h = np.concatenate((gammaY_h,np.array(gammaX_h[-1,:]).reshape((1,1))),axis=0)

# iterate to find ideal starting gamma for a
while fa_h[-1,0] > fa_h[0,0]:
    gammaA_h = np.concatenate((gammaA_h,np.array(gammaA_h[-1,:]/2.0).reshape((1,1))),axis=0)

    a_t = a0 - (1/lam)*gammaA_h[-1,0]*(dfaY_p_h[-1,:] - dfaX_p_h[-1,:])
    a_h = np.concatenate((a_h,np.reshape(a_t,(1,len(a_t)))),axis=0)

    [f_t,_] = gd_FOMp(a_h[-1,:],X_h[-1,:])
    fX_p_h = np.concatenate((fX_p_h,np.array(f_t).reshape((1,1))),axis=0)    

    [f_t,_] = gd_FOMp(a_h[-1,:],Y_h[-1,:])
    fY_p_h = np.concatenate((fY_p_h,np.array(f_t).reshape((1,1))),axis=0)  

    fa_h = np.concatenate((fa_h,np.array((1/lam)*(fY_p_h[-1,0]-fX_p_h[-1,0])).reshape((1,1))), axis=0)
    PrintFoM(fa_h[-1,0],'ap')


# start main optimization loop
while True:
    ii=1
    while (fX_p_h[-1,0] < fX_p_h[-2,0]) and (fY_h[-1,0] < fY_h[-2,0]) and (fa_h[-1,0] <= fa_h[-2,0]):
        print('Iterating '+str(ii))

        # iterate X
        X_t = X_h[-1,:] - gammaX_h[-1,0]*dfX_p_h[-1,:]
        X_h = np.concatenate((X_h,np.reshape(X_t,(1,len(X_t)))),axis=0)
        # iterate Y
        Y_t = Y_h[-1,:] - gammaY_h[-1,0]*(dfY_p_h[-1,:] + lam*dfY_e_h[-1,:])
        Y_h = np.concatenate((Y_h,np.reshape(Y_t,(1,len(Y_t)))),axis=0)
        # iterate a
        a_t = a_h[-1,:] - (1/lam)*gammaA_h[-1,0]*(dfaY_p_h[-1,:] - dfaX_p_h[-1,:])
        a_h = np.concatenate((a_h,np.reshape(a_t,(1,len(a_t)))),axis=0)        

        # compute FoM, Fp(X) , Fp(Y) , Fe(Y)
        [f_t,_] = gd_FOMp(a_h[-1,:],X_h[-1,:])
        fX_p_h = np.concatenate((fX_p_h,np.array(f_t).reshape((1,1))),axis=0)
        PrintFoM(f_t,'Xp')
        [f_t,_] = gd_FOMp(a_h[-1,:],Y_h[-1,:])
        fY_p_h = np.concatenate((fY_p_h,np.array(f_t).reshape((1,1))),axis=0)
        PrintFoM(f_t,'Yp')
        f_t = gd_FOMe(a_h[-1,:],Y_h[-1,:])
        fY_e_h = np.concatenate((fY_e_h,np.array(f_t).reshape((1,1))), axis=0)
        PrintFoM(f_t,'Ye')
        fY_h = np.concatenate((fY_h,np.array(fY_p_h[-1,0]+lam*fY_e_h[-1,0]).reshape((1,1))), axis=0)
        fa_h = np.concatenate((fa_h,np.array((1/lam)*(fY_p_h[-1,0]-fX_p_h[-1,0])).reshape((1,1))), axis=0)        

        # if we have iterated too much, start increasing gamma
        ii += 1
        if (ii > 20):
            gammaX_h = np.concatenate((gammaX_h,np.array(gammaX_h[-1,:]*2.0).reshape((1,1))),axis=0)
            gammaY_h = np.concatenate((gammaY_h,np.array(gammaY_h[-1,:]*2.0).reshape((1,1))),axis=0)
            gammaA_h = np.concatenate((gammaA_h,np.array(gammaA_h[-1,:]*2.0).reshape((1,1))),axis=0)

    # recompute adjoint equations for new search direction
    print("Recomputing adjoint equations")

    # grab last good settings
    a_h = np.concatenate((a_h,np.reshape(a_h[-2,:],(1,len(a_h[-2,:])))),axis=0)
    X_h = np.concatenate((X_h,np.reshape(X_h[-2,:],(1,len(X_h[-2,:])))),axis=0)
    Y_h = np.concatenate((Y_h,np.reshape(Y_h[-2,:],(1,len(Y_h[-2,:])))),axis=0)
    fX_p_h = np.concatenate((fX_p_h,np.array(fX_p_h[-2,0]).reshape((1,1))),axis=0)
    fY_h = np.concatenate((fY_h,np.array(fY_h[-2,0]).reshape((1,1))),axis=0)
    fa_h = np.concatenate((fa_h,np.array(fa_h[-2,0]).reshape((1,1))),axis=0)

    # calculate adjoint equations for X
    z_p,y_p,yadj_p,_,k_p = gd_adj_p(a_h[-1,:],X_h[-1,:])
    O0_p,N0_p,ACT0_p = get_ON_and_ACT(z_p,y_p,yadj_p,k_p,physics_params)
    # calculate new grdient for X and a
    dfaX_t,dfX_t = gd_dFOMp(a_h[-1,:],X_h[-1,:],z_p,y_p,yadj_p,O0_p,N0_p,ACT0_p)
    dfX_p_h = np.concatenate((dfX_p_h,np.reshape(dfX_t,(1,len(dfX_t)))),axis=0)
    dfaX_p_h = np.concatenate((dfaX_p_h,np.reshape(dfaX_t,(1,len(dfaX_t)))),axis=0)
    # calculate adjoint equations for Y
    z_p,y_p,yadj_p,_,k_p = gd_adj_p(a_h[-1,:],Y_h[-1,:])
    O0_p,N0_p,ACT0_p = get_ON_and_ACT(z_p,y_p,yadj_p,k_p,physics_params)  
    z_e,y_e,yadj_e,_,k_e = gd_adj_e(a_h[-1,:],Y_h[-1,:])
    O0_e,N0_e,ACT0_e = get_ON_and_ACT(z_e,y_e,yadj_e,k_e,physics_params) 
    # calculate new gradient for Y and a     
    dfaY_t,dfY_t = gd_dFOMp(a_h[-1,:],Y_h[-1,:],z_p,y_p,yadj_p,O0_p,N0_p,ACT0_p)
    dfaYe_t,dfYe_t = gd_dFOMe(a_h[-1,:],Y_h[-1,:],z_e,y_e,yadj_e,O0_e,N0_e,ACT0_e)
    dfY_p_h = np.concatenate((dfY_p_h,np.reshape(dfY_t,(1,len(dfY_t)))),axis=0)
    dfaY_p_h = np.concatenate((dfaY_p_h,np.reshape(dfaY_t,(1,len(dfaY_t)))),axis=0)
    dfY_e_h = np.concatenate((dfY_e_h,np.reshape(dfYe_t,(1,len(dfYe_t)))),axis=0)

    # calculate adjoint equations for Y

    if (ii == 2): # meaning no improvement from recalculating gradient
        # change gamma
        print('Updating gamma')

        f0n = fX_p_h[-1,0]
        X0n = X_h[-1,:]
        while fX_p_h[-1,0] >= f0n:
            gammaX_h = np.concatenate((gammaX_h,np.array(gammaX_h[-1,:]/2.0).reshape((1,1))),axis=0)

            # iterate X
            X_t = X0n - gammaX_h[-1,0]*dfX_p_h[-1,:]
            X_h = np.concatenate((X_h,np.reshape(X_t,(1,len(X_t)))),axis=0)

            # compute FoM
            [f_t,_] = gd_FOMp(a_h[-1,:],X_h[-1,:])
            fX_p_h = np.concatenate((fX_p_h,np.array(f_t).reshape((1,1))),axis=0)
            PrintFoM(f_t)
        gammaY_h =  np.concatenate((gammaY_h,np.array(gammaX_h[-1,:]).reshape((1,1))), axis=0)

        f0n = fa_h[-1,0]
        a0n = a_h[-1,:]
        while fa_h[-1,0] >= f0n:
            gammaA_h = np.concatenate((gammaA_h,np.array(gammaA_h[-1,:]/2.0).reshape((1,1))),axis=0)

            # iterate a
            a_t = a0n - (1/lam)*gammaA_h[-1,0]*(dfaY_p_h[-1,:] - dfaX_p_h[-1,:])
            a_h = np.concatenate((a_h,np.reshape(a_t,(1,len(a_t)))),axis=0)  

            # compute FoM
            [fX_t,_] = gd_FOMp(a_h[-1,:],X_h[-1,:])
            [fY_t,_] = gd_FOMp(a_h[-1,:],Y_h[-1,:])
            fa_h = np.concatenate((fa_h,np.array((1/lam)*(fY_t-fX_t)).reshape((1,1))), axis=0) 
            PrintFoM(fX_t,'Xp')
            PrintFoM(fY_t,'Yp')             

    if fX_p_h[-1,0] < 5e-18:
        break

    if len(fX_p_h[:,0]) > 5000:
        break

SaveData("testx3",[a0,X_h,f_h,dfX_p_h,gammaX_h])






















   

