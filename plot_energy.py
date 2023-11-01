import numpy as np
import matplotlib.pyplot as plt

# times, E = np.loadtxt("Practice Session 3/two_body_motion_energy.txt", unpack=True)
E = np.loadtxt("Astrofysische_Groep_1/two_body_motion_energy_RK2.txt", unpack=True)
times = range(len(E))

E_rel_error = (E-E[0])/E[0]
plt.plot(times, E_rel_error)
plt.xlabel('timestep')
plt.ylabel('Energy error')
plt.show()