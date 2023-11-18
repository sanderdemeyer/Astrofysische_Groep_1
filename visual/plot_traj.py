import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import random
import mpl_toolkits.mplot3d.axes3d as p3
from scipy import interpolate
import os

directory = os.getcwd()

def trajectories(file, dim):
    """
    This function takes a txt file where the lines contains the timesteps and the positions of the N partilces at that timestep and gives the timesteps and a numpy array with the coordinates of each particle at a timestep.
    
    Parameters
    -----------------------------
    file
        The txt file containing the times teps and positions
        The input needs to have the form: t  x_1  y_1  z_1 ... x_N  y_N  z_N
    dim
        Dimension of the Nbody simulation
    
    Output
    -----------------------------
    t
        An array containing the timesteps
    ind_traj
        An array defining the traject of the N particles.
        The output ind_traj has the form [pos_0, ..., pos_t] with pos_i=[[x_0i, y_0i, z_0i], ..., [x_Ni, y_Ni, z_Ni]]
    """
    traj = np.loadtxt('{}/traj/{}'.format(directory, file))
    # get the number of particles and timesteps
    N = int((len(traj[0]) - 1)/dim)
    n_t = len(traj)
    # getting the timesteps
    t=np.zeros(n_t)
    # an arry with the trajectories of the N particles
    ind_traj = np.zeros((n_t, N, dim))
    for i in range(n_t):
        t[i] = traj[i][0]
        part_traj = np.split(traj[i][1:], N)
        for j in range(N):
            ind_traj[i][j]= part_traj[j]
    return t, ind_traj

def animate(file, dim, name='animation', lim= None, dpi=300):
    """
    This function takes a txt file where the lines contains the timesteps and the positions of the N partilces at that timestep and animates the trajectories of the particles.
    
    Parameters
    -----------------------------
    file
        The txt file containing the times teps and positions
        The input needs to have the form: t  x_1  y_1  z_1 ... x_N  y_N  z_N
    dim
        Dimension of the Nbody simulation, has to be 2 or 3
    name
        Name of the animation file
    lim
        An array with the axis limits
    dpi
        DPI of the animation
    
    Output
    -----------------------------
    Creates an mkv file a animating the trajectories of the N particles.
    """
    t, traj = trajectories(file, dim)
    n_t = len(traj)
    N = len(traj[0])
    def update_2(i):
        for j in range(N):
            scatters[j].set_offsets([traj[i][j][0], traj[i][j][1]])
        return scatters

    def update_3(i):
        for j in range(N):
            scatters[j]._offsets3d=(traj[i][j,0:1], traj[i][j,1:2], traj[i][j,2:])
        return scatters

    fig = plt.figure()
    # TODO: add a time slider
    if dim == 2:
        ax = fig.add_subplot(111)
        ax.set_title('N={}'.format(N))
        if lim:
            ax.set_xlim(lim[0][0], lim[0][1])
            ax.set_ylim(lim[1][0], lim[1][1])
        ax.grid(False)
        scatters = [ax.scatter(traj[0][i][0], traj[0][i][1], s=2) for i in range(N)]
        ani = animation.FuncAnimation(fig, update_2, frames=n_t)
        ani.save('{}/ani_traj/{}_2D.mkv'.format(directory,name), fps=30, dpi=dpi)

    if dim == 3:
        ax = fig.add_subplot(projection='3d')
        ax.set_title('N={}'.format(N))
        if lim:
            ax.set_xlim3d(lim[0][0], lim[0][1])
            ax.set_ylim3d(lim[1][0], lim[1][1])
            ax.set_zlim3d(lim[2][0], lim[2][1])
        ax.grid(False)
        scatters = [ax.scatter(traj[0][i][0], traj[0][i][1], traj[0][i][2], s=2) for i in range(N)]
        ani = animation.FuncAnimation(fig, update_3, frames=n_t)
        ani.save('{}/ani_traj/{}_3D.mkv'.format(directory,name), fps=30, dpi=dpi)
    plt.show()

def plot(file, dim, name='trajectories', lim =None, dpi=300):
    """
    This function takes a txt file where the lines contains the timesteps and the positions of the N partilces at that timestep and plots the trajectories of the particles in 2D or 3D.
    
    Parameters
    -----------------------------
    file
        The txt file containing the times teps and positions
        The input needs to have the form: t  x_1  y_1  z_1 ... x_N  y_N  z_N
    dim
        Dimension of the Nbody simulation, has to be 2 or 3
    name
        Name of the plot
    lim
        An array with the axis limits
    dpi
        DPI of the image
    
    Output
    -----------------------------
    Creates a png file plotting the trajectories of the N particles.
    """
    t, traj = trajectories(file, dim)
    n_t = len(traj)
    N = len(traj[0])
    fig = plt.figure()
    if dim == 2:
        x = np.zeros((N, n_t))
        y = np.zeros((N, n_t))
        for i in range(N):
            for j in range(n_t):
                x[i][j] = traj[j][i][0]
                y[i][j] = traj[j][i][1]
        ax = fig.add_subplot(111)
        ax.set_title('N={}'.format(N))
        if lim:
            ax.set_xlim(lim[0][0], lim[0][1])
            ax.set_ylim(lim[1][0], lim[1][1])
        ax.grid(False)
        for i in range(N):
            ax.plot(x[i],y[i])
        plt.savefig('{}/plot_traj/{}_2D.png'.format(directory,name), dpi=dpi)

    if dim == 3:
        x = np.zeros((N, n_t))
        y = np.zeros((N, n_t))
        z = np.zeros((N, n_t))
        for i in range(N):
            for j in range(n_t):
                x[i][j] = traj[j][i][0]
                y[i][j] = traj[j][i][1]
                z[i][j] = traj[j][i][2]
        ax = fig.add_subplot(projection='3d')
        ax.set_title('N={}'.format(N))
        if lim:
            ax.set_xlim3d(lim[0][0], lim[0][1])
            ax.set_ylim3d(lim[1][0], lim[1][1])
            ax.set_zlim3d(lim[2][0], lim[2][1])
        for i in range(N):
            ax.plot(x[i],y[i], z[i])
        plt.savefig('{}/plot_traj/{}_3D.png'.format(directory,name), dpi=dpi)
    plt.show()

print('This program animtes or plots the trajectories of the particles in an N-body simulation.')
integrator = input('Please provide the file:\n')
type_plot = input('Would you like to animate or plot the trajectories:\n')
dim = int(input('What is the dimension of the simulation:\n'))
lim_nolim = input('Would you like to limit the axis of the simulation?\n')
if lim_nolim == 'y':
    lim = []
    x_in = input('xlim: ')
    xlim = x_in.split(",")
    lim.append([int(xlim[0]), int(xlim[1])])
    y_in = input('ylim: ')
    ylim = y_in.split(",")
    lim.append([int(ylim[0]), int(ylim[1])])
    if dim == 3:
        z_in = input('zlim: ')
        zlim = z_in.split(",")
        lim.append([int(zlim[0]), int(zlim[1])])
else:
    lim = None

name = integrator.rstrip('.txt')
if type_plot == 'animate':
    animate(integrator, dim, name= name, lim=lim)
elif type_plot == 'plot':
    plot(integrator, dim, name= name, lim=lim)
elif type_plot == 'both':
    animate(integrator, dim, name= name, lim=lim)
    plot(integrator, dim, name= name, lim=lim)