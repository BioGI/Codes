import matplotlib.pyplot as plt
import numpy as np

sl = np.loadtxt('singleLattice/output')
dl = np.loadtxt('dualLattice/interpolation-00001.dat')

# for i in range(15):
#     fig = plt.figure()
#     plt.plot(dl[i::30,-1]-dl[i+15::30,-1])


fig = plt.figure()
plt.plot(sl[:,1]-dl[:,1])

fig = plt.figure()
plt.plot(sl[:,0]-dl[:,0])

plt.show()




