import numpy as np
import matplotlib.pyplot as plt

Systems = ["4-body-star", "binary", "rings", "two-body-plane"]

"""
plt.figure(figsize=(8,6))
plt.title("RK45 with $h_{start} = 0.001$")
#t, E =  np.loadtxt('energy/4-body-star/4-body-star_RK45_100.000000_0.010000.txt', unpack=True)
#E_rel_error = abs((E - E[0])/E[0])
#plt.plot(t, E_rel_error, label="4-body-star, $h=0.01$")
for system in Systems:
    t, E =  np.loadtxt('energy/{}/{}_RK45_100.000000_0.001000.txt'.format(system, system), unpack=True)
    E_rel_error = abs((E - E[0])/E[0])
    plt.plot(t, E_rel_error, label=system)
plt.yscale('log')
plt.xlim(0, 50)
plt.xlabel('t')
plt.ylabel('$\Delta E_{rel}$')
plt.legend(loc="upper right")
plt.savefig("plot_energy/RK45.png")
"""

Integrators = ["Forward Euler", "RK2", "Heun", "Heun3", "Ralston", "Ralston3", "RK3", "RK4", "Forest Ruth", "PEFRL", "Velocity Verlet", "Position Verlet", "Leapfrog", "Yoshida 4"]

Systems_sel = ["binary", "rings"]

H = [0.01, 0.001, 0.0001]

for integ in Integrators:
    plt.figure(figsize=(8,6))
    plt.title("{}".format(integ))
    for system in Systems_sel:
        for h in H:
            t, E =  np.loadtxt('energy/{}/{}_{}_50.000000_{:.6f}.txt'.format(system, system, integ, h), unpack=True)
            E_rel_error = abs((E - E[0])/E[0])
            plt.plot(t, E_rel_error, label="{}, h={}".format(system, str(h)))
    plt.yscale('log')
    plt.xlim(0, 50)
    plt.xlabel('t')
    plt.ylabel('$\Delta E_{rel}$')
    plt.legend(loc="upper right")
    plt.savefig("plot_energy/{}.png".format(integ))