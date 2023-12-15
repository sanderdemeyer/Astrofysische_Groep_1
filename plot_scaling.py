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

plt.figure(figsize=(8,6))
for friend in files_friend:
    # to plot the expected behavior
    times = np.loadtxt(friend, unpack=True)
    l = str(friend).removesuffix('.txt').lstrip('scaling_friend\\')
    par, cov = sp.optimize.curve_fit(fun, num, times)
    plt.scatter(num, times, label=l, s=5)
    plt.plot(t,fun(t, par[0], par[1]), '--')
plt.xlabel('N')
plt.ylabel('t per steps [$\mu s$]')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
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
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.savefig('performance/figs/scaling_general.png')
plt.show()
