import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from tabulate import tabulate
from moment_equations_util import *
from moment_equations import *


def PrintTable(ary):
    table = tabulate(ary)
    print(table)

def SaveData(fname,ary):
    np.save("Data/"+fname,ary)

def LoadData(fname):
    return np.load("Data/"+fname+".npy",allow_pickle=True)

def PlotMoments(z,y,k,fig=None,ax=None,colors=None,linestyles=None):
    if not fig or not ax:
        fig, ax = plt.subplots()
    if not colors:
        colors = ['C0','C1']
    if not linestyles:
        linestyles = ['-','-']

    # plot moments
    xrms=(y[0,:]+y[1,:])
    yrms=(y[0,:]-y[1,:])
    ax.plot(z,xrms,label='$\langle x^2 \\rangle$',color=colors[0],linestyle=linestyles[0])
    ax.plot(z,yrms,label='$\langle y^2 \\rangle$',color=colors[1],linestyle=linestyles[1])

    ax.legend()
    ax.grid(True)
    ax.set_ylabel('moments $[m^2]$')
    ax.set_xlabel('Z position [m]')    

    # plot k values
    ax2 = ax.twinx()
    ax2.plot(z,k,color='k',linestyle=linestyles[0])
    ax2.set_ylim([-0.05,0.05])
    ax2.set_yticklabels([])