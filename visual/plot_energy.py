import numpy as np
import matplotlib.pyplot as plt
import os
directory = os.getcwd()

print(directory)
print('This program plots the total energy or the realtive energy error of the N-body system.')
#integrator = input('Please provide the file:\n')

integrator = "RK4"

#t, E = np.loadtxt('{}/energy/{}'.format(directory, integrator), unpack=True)
#t, E = np.loadtxt('energy/perturbed-criss-cross_RK4_100.000000_0.100000.txt', unpack=True)
E = np.loadtxt('energy_reg/Reg_2D_Burrau_scaled_RK4_10000000_0.001000.txt', unpack=True)
E = np.loadtxt('energy_reg/Reg_2D_collision_RK4_1000000_0.001000_0.500000.txt', unpack=True)



t = range(len(E))
#type_plot = input('Would you like to plot the total energy or the relative error:\n')
type_plot = "total energy"
E_rel_error = (E - E[0])/E[0]

#name = integrator.rstrip('.txt')
#args = name.split('_')
#filename = args[0] + '_' + args[1]

if type_plot == 'total energy':
    #plt.plot(t, E)
    plt.plot(t, [np.log10(abs(e)) for e in E])
    plt.xlabel('t')
    plt.ylabel('E')
    plt.grid()
    #plt.savefig('{}/plot_energy/{}_total.png'.format(directory, filename), dpi=300)
    plt.show()
if type_plot == 'relative error':
    plt.plot(t, E_rel_error)
    plt.xlabel('t')
    plt.ylabel('$\Delta E_{rel}$')
    plt.grid()
    #plt.savefig('{}/plot_energy/{}_relerror.png'.format(directory,filename), dpi=300)
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
    plt.savefig('{}/plot_energy/{}.png'.format(directory,filename), dpi=300)
    plt.show()