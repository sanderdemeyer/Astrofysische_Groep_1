import numpy as np
import matplotlib.pyplot as plt
import os
directory = os.getcwd()

print(directory)
print('This program plots the total energy or the realtive energy error of the N-body system.')
integrator = input('Please provide the file:\n')
t, E = np.loadtxt('{}/energy/{}'.format(directory, integrator), unpack=True)
type_plot = input('Would you like to plot the total energy or the relative error:\n')
E_rel_error = (E - E[0])/E[0]
name = integrator.rstrip('.txt')
print(integrator)

if type_plot == 'total energy':
    plt.plot(t, E)
    plt.xlabel('t')
    plt.ylabel('E')
    plt.grid()
    plt.savefig('{}/plot_energy/{}_total.png'.format(directory,name), dpi=300)
    plt.show()
if type_plot == 'relative error':
    plt.plot(t, E_rel_error)
    plt.xlabel('t')
    plt.ylabel('$\Delta E_{rel}$')
    plt.grid()
    plt.savefig('{}/plot_energy/{}_relerror.png'.format(directory,name), dpi=300)
    plt.show()
if type_plot == 'both':
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(t, E)
    ax1.set(xlabel='t', ylabel='E')
    ax2.grid()
    
    ax2.plot(t, E_rel_error)
    ax2.set(xlabel='t', ylabel='$\Delta E_{rel}$')
    ax2.grid()

    fig.tight_layout()
    plt.savefig('{}/plot_energy/{}.png'.format(directory,name), dpi=300)
    plt.show()