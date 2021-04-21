from gd_adj_match_1var_functions import RunMoments,getLattice,Initial_conditions
from analysis_util import *


if 0:
    dt = LoadData('test')
    period = 5

    lat,_ = getLattice(dt[0][0,:],period)
    z,y,m,k = RunMoments(lat,Initial_conditions())

    fig,ax=plt.subplots()
    #PlotMoments(z,y,k,fig,ax,linestyles=['--','--'])

    lat,_ = getLattice(dt[0][-1,:],period)
    z,y,m,k = RunMoments(lat,Initial_conditions())

    PlotMoments(z,y,k,fig,ax)

    plt.figure()
    plt.plot(np.log10(dt[1]/dt[1][1,0]),marker='.')
    plt.ylabel('$log( FoM/FoM_0 )$')
    plt.xlabel('Iterations (kind of)')
    plt.title('FoM vs Iterations')
    plt.grid(True)

else:
    dt = LoadData('testx2')
    period = 5

    lat,_ = getLattice(dt[0],period)
    z,y,m,k = RunMoments(lat,dt[1][0,:])

    fig,ax=plt.subplots()
    PlotMoments(z,y,k,fig,ax,linestyles=['--','--'])

    lat,_ = getLattice(dt[0],period)
    z,y,m,k = RunMoments(lat,dt[1][-1,:])

    PlotMoments(z,y,k,fig,ax)

    plt.figure()
    plt.plot(np.log10(dt[2]/dt[2][1,0]),marker='.')
    plt.ylabel('$log( FoM/FoM_0 )$')
    plt.xlabel('Iterations (kind of)')
    plt.title('FoM vs Iterations')
    plt.grid(True)


plt.show()


plt.show()