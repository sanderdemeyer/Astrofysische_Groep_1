import numpy as np
import pathlib
import matplotlib.pyplot as plt
import scipy as sp

files_friend = [f for f in pathlib.Path('scaling_friend').iterdir()]
files_general = [f for f in pathlib.Path('scaling_general').iterdir()]

num = np.arange(2, 100)

def fun(x,n,b):
    return x**n + b

t = np.arange(0, 100, 0.1)

plt.figure(figsize=(8,6))
for friend in files_friend:
    # to plot the expected behavior
    times = np.loadtxt(friend, unpack=True)
    l = str(friend).removesuffix('.txt').lstrip('scaling_friend\\')
    par, cov = sp.optimize.curve_fit(fun, num, times)
    label = l + ": $x^{} + {:.0f}$".format(int(par[0]), par[1])
    plt.scatter(num, times, label=label, s=5)
    plt.plot(t,fun(t, par[0], par[1]), '--')
plt.xlabel('N')
plt.ylabel('t per 10000 steps [ms]')
plt.legend()
plt.savefig('performance/scaling_friend.png')
plt.show()

'''
plt.figure(figsize=(8,6))
for friend in files_general:
    # to plot the expected behavior
    times = np.loadtxt(friend, unpack=True)
    l = str(friend).removesuffix('.txt').lstrip('scaling_general\\')
    par, cov = sp.optimize.curve_fit(fun, num, times)
    label = l + ": $x^{} + {:.0f}$".format(int(par[0]), par[1])
    plt.scatter(num, times, label=label, s=5)
    plt.plot(t,fun(t, par[0], par[1]), '--')
plt.xlabel('N')
plt.ylabel('t per 10000 steps [ms]')
plt.legend()
plt.savefig('performance/scaling_general.png')
plt.show()
'''