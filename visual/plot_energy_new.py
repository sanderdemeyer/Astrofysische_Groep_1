import numpy as np
import matplotlib.pyplot as plt
import os

Integrators = ["Forward Euler", "RK2", "Heun", "Heun3", "Ralston", "Ralston3", "RK3", "RK4", "Forest Ruth", "PEFRL", "Velocity Verlet", "Position Verlet", "Leapfrog", "Yoshida 4"]
H = [0.1, 0.01, 0.001]

system = "Solar-System-New"

if not os.path.exists("plot_energy/{}".format(system)):
    os.mkdir("plot_energy/{}".format(system))

for integ in Integrators:
    plt.figure(figsize=(8,6))
    plt.title("{} with {}".format(system, integ))
    for h in H:
        t, E =  np.loadtxt('energy/{}/{}_{}_50.000000_{:.6f}.txt'.format(system, system, integ, h), unpack=True)
        E_rel_error = abs((E - E[0])/E[0])
        plt.plot(t, E_rel_error, label=str(h))
    plt.yscale('log')
    plt.xlim(0, 50)
    plt.xlabel('t')
    plt.ylabel('$\Delta E_{rel}$')
    plt.legend(loc="lower right")
    plt.savefig("plot_energy/{}/{}_{}.png".format(system, system, integ))
