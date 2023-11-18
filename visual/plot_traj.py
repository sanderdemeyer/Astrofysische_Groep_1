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

def animate(file, dim, tstep, name='animation', lim= None, dpi=300):
    """
    This function takes a txt file where the lines contains the timesteps and the positions of the N partilces at that timestep and animates the trajectories of the particles.
    
    Parameters
    -----------------------------
    file
        The txt file containing the times teps and positions
        The input needs to have the form: t  x_1  y_1  z_1 ... x_N  y_N  z_N
    dim
        Dimension of the Nbody simulation, has to be 2 or 3
    tstep
        The timescale of the simulation, or in other words how much should one second of video be equal to the time used in simulation
    name
        Name of the animation file
    dpi
        DPI of the animation
    
    Output
    -----------------------------
    Creates an mkv file a animating the trajectories of the N particles.
    """
    t, traj_org = trajectories(file, dim)
    t_max = t[-1]
    tstep_per_frame = tstep/30
    time = np.arange(0, t_max, step=tstep_per_frame)
    f = interpolate.interp1d(t, traj_org, axis=0)
    traj = f(time)
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
        x_all = np.array([traj[i][:,0] for i in range(n_t)]).flatten()
        y_all = np.array([traj[i][:,1] for i in range(n_t)]).flatten()
        ax.set_xlim(x_all.min(), x_all.max())
        ax.set_ylim(y_all.min(), y_all.max())
        ax.grid(False)
        scatters = [ax.scatter(traj[0][i][0], traj[0][i][1], s=2) for i in range(N)]
        ani = animation.FuncAnimation(fig, update_2, frames=n_t)
        ani.save('{}/ani_traj/{}_2D.mkv'.format(directory,name), fps=30, dpi=dpi)

    if dim == 3:
        ax = fig.add_subplot(projection='3d')
        ax.set_title('N={}'.format(N))
        x_all = np.array([traj[i][:,0] for i in range(n_t)]).flatten()
        y_all = np.array([traj[i][:,1] for i in range(n_t)]).flatten()
        z_all = np.array([traj[i][:,2] for i in range(n_t)]).flatten()
        ax.set_xlim3d(x_all.min(), x_all.max())
        ax.set_ylim3d(y_all.min(), y_all.max())
        ax.set_zlim3d(z_all.min(), z_all.max())
        ax.grid(False)
        scatters = [ax.scatter(traj[0][i][0], traj[0][i][1], traj[0][i][2], s=2) for i in range(N)]
        ani = animation.FuncAnimation(fig, update_3, frames=n_t)
        ani.save('{}/ani_traj/{}_3D.mkv'.format(directory,name), fps=30, dpi=dpi)
    plt.show()

def plot(file, dim, name='trajectories', dpi=300):
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
        x_all = x.flatten()
        y_all = y.flatten()
        ax.set_xlim(x_all.min(), x_all.max())
        ax.set_ylim(y_all.min(), y_all.max())
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
        x_all = x.flatten()
        y_all = y.flatten()
        z_all = z.flatten()
        ax.set_xlim3d(x_all.min(), x_all.max())
        ax.set_ylim3d(y_all.min(), y_all.max())
        ax.set_zlim3d(z_all.min(), z_all.max())
        for i in range(N):
            ax.plot(x[i],y[i], z[i])
        plt.savefig('{}/plot_traj/{}_3D.png'.format(directory,name), dpi=dpi)
    plt.show()

print('This program animtes or plots the trajectories of the particles in an N-body simulation.')
integrator = input('Please provide the file:\n')
type_plot = input('Would you like to animate or plot the trajectories:\n')
dim = int(input('What is the dimension of the simulation:\n'))
tstep = int(input('How much should one second of video be equal to the time used in simulation: \n'))

name = integrator.rstrip('.txt')
if type_plot == 'animate':
    animate(integrator, dim, tstep, name= name)
elif type_plot == 'plot':
    plot(integrator, dim,name= name)
elif type_plot == 'both':
    animate(integrator, dim, tstep, name= name)
    plot(integrator, dim, name= name)