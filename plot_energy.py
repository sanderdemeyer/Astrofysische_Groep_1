import numpy as np
import matplotlib.pyplot as plt

E = np.loadtxt("two_body_motion_energy_RK2.txt", unpack=True)
times = range(len(E))

E_rel_error = (E-E[0])/E[0]
plt.plot(times, E_rel_error)
plt.xlabel('timestep')
plt.ylabel('Energy error')
plt.grid()
plt.show()