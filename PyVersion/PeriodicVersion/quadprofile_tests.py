import matplotlib.pyplot as plt
from moment_equations_util import *

# number of quads
numQuads = 6

amplitude = 1.0 # quadrupole amplitude
qlength = [0.2] # quadrupole length
dlength = 0.2 # drift space length
polarity = [1,-1] # quadrupole polarity

stepsize=0.01


z,y = CreateQuadProfile(amplitude,qlength,dlength,numQuads,polarity,stepsize)


plt.figure()
plt.plot(z,y)
plt.ylabel('$K_q$ (a.u.)')
plt.xlabel('Z position [m]')
plt.show()
