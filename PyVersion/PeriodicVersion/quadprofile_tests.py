import matplotlib.pyplot as plt
from moment_equations_util import *

# number of quads
repeat = 5

amplitude = [1.0,1.0,1.0] # quadrupole amplitude
qlength = [0.05,0.1,0.05] # quadrupole length
dlength = [0.0,0.2,0.2] # drift space length
polarity = [1,-1,1] # quadrupole polarity


lattice = CreateLatticeProfile(amplitude,qlength,dlength,polarity,repeat)
print(lattice)

PlotLatticeProfile(lattice)
