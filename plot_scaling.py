import numpy as np
import pathlib
import matplotlib.pyplot as plt
import scipy as sp

files_friend = [f for f in pathlib.Path('scaling_friend').iterdir()]
files_general = [f for f in pathlib.Path('scaling_general').iterdir()]

num = np.arange(2, 100)

def fun(x,a,b):
    return a*(x**2) + b

t = np.arange(0, 100, 0.1)

symplectic = ["Leapfrog", "Velocity Verlet", "Position Verlet", "Yoshida 4", "Forest Ruth", "PEFRL"]
primitive = ["Heun", "Heun 3", "Ralston", "Ralston3", "Forward Euler", "RK2", "RK3", "RK4"]

line_symp = []
line_prim = []

plt.figure(figsize=(8,6))
for friend in files_friend:
    # to plot the expected behavior
    times = np.loadtxt(friend, unpack=True)
    l = str(friend).removesuffix('.txt').lstrip('scaling_friend\\')
    par, cov = sp.optimize.curve_fit(fun, num, times)
    p = plt.scatter(num, times, label=l, s=5)
    if l in symplectic:
        line_symp.append(p)
    elif l in primitive:
        line_prim.append(p)
    plt.plot(t,fun(t, par[0], par[1]), '--')
plt.xlabel('N')
plt.ylabel('t per steps [$\mu s$]')
#plt.legend(loc='upper left')
l1 = plt.legend(line_symp, symplectic, loc=(0.23, 0.733))
l2 = plt.legend(line_prim, primitive, loc=(0.01, 0.693))
plt.gca().add_artist(l1)
plt.tight_layout()
plt.savefig('performance/figs/scaling_friend.png')
plt.show()


plt.figure(figsize=(8,6))
for friend in files_general:
    # to plot the expected behavior
    times = np.loadtxt(friend, unpack=True)
    l = str(friend).removesuffix('.txt').lstrip('scaling_general\\')
    par, cov = sp.optimize.curve_fit(fun, num, times)
    label = l + ": $  \propto x^{}$".format(int(par[0]))
    plt.scatter(num, times, label=l, s=5)
    plt.plot(t,fun(t, par[0], par[1]), '--')
plt.xlabel('N')
plt.ylabel('t per steps [$\mu s$]')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig('performance/figs/scaling_general.png')
plt.show()
