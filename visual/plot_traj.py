import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
import random
import mpl_toolkits.mplot3d.axes3d as p3
from scipy import interpolate
import os
import warnings
# ignore this specific warning
warnings.filterwarnings("ignore", message="Attempting to set identical low and high zlims makes transformation singular; automatically expanding.")

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

def animate(file, dim, tstep, title, lim=None, line=False, project=False, label='auto', name='animation', dpi=300):
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
    title
        Title of the plot
    lim
        Limit the axis
    line
        Whether or not to plot the lines
    project
        Whether or not to plot a 2D projection
    label
        Labels of the bodies. Defaults to automatic numbering
    name
        Name of the animation file
    dpi
        DPI of the animation
    
    Output
    -----------------------------
    Creates a gif animating the trajectories of the N particles.
    """
    t, traj_org = trajectories(file, dim)
    n_t_org = len(traj_org)
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

    plt.style.use('dark_background')
    fig = plt.figure()
    # TODO: add a time slider
    if dim == 2:
        ax = fig.add_subplot(111)
        ax.set_title('N={}, '.format(N) + title)
        x_all = np.array([traj[i][:,0] for i in range(n_t)]).flatten()
        y_all = np.array([traj[i][:,1] for i in range(n_t)]).flatten()
        if lim:
            ax.set_xlim(lim[0])
            ax.set_ylim(lim[1])
        else:
            ax.set_xlim(x_all.min(), x_all.max())
            ax.set_ylim(y_all.min(), y_all.max())
        ax.set_xlabel("x [AU]")
        ax.set_ylabel("y [AU]")
        scatters = [ax.scatter(traj[0][i][0], traj[0][i][1], s=2) for i in range(N)]
        if label=='auto':
            label = ['Body {}'.format(i) for i in range(1, N+1)]
            ax.legend(label, loc='center left', bbox_to_anchor=(1.2, 0.5))
        elif label == None:
            ax.legend_=None
        else:
            ax.legend(label, loc='center left', bbox_to_anchor=(1.2, 0.5))
        fig.set_tight_layout(True)
        if line:
            x = np.zeros((N, n_t_org))
            y = np.zeros((N, n_t_org))
            for i in range(N):
                for j in range(n_t_org):
                    x[i][j] = traj_org[j][i][0]
                    y[i][j] = traj_org[j][i][1]
            for i in range(N):
                ax.plot(x[i],y[i], linestyle='dotted', linewidth=0.3)
        ax.grid(False)
        ani = animation.FuncAnimation(fig, update_2, frames=n_t, repeat=False)
        ani.save('{}/ani_traj/{}_2D.gif'.format(directory,name), fps=30, dpi=dpi)

    if dim == 3:
        ax = fig.add_subplot(projection='3d')
        ax.set_title('N={}, '.format(N) + title)
        x_all = np.array([traj[i][:,0] for i in range(n_t)]).flatten()
        y_all = np.array([traj[i][:,1] for i in range(n_t)]).flatten()
        z_all = np.array([traj[i][:,2] for i in range(n_t)]).flatten()
        if lim:
            ax.set_xlim3d(lim[0])
            ax.set_ylim3d(lim[1])
            ax.set_zlim3d(lim[1])
        else:
            ax.set_xlim3d(x_all.min(), x_all.max())
            ax.set_ylim3d(y_all.min(), y_all.max())
            ax.set_zlim3d(z_all.min(), z_all.max())
        ax.set_xlabel("x [AU]")
        ax.set_ylabel("y [AU]")
        ax.set_zlabel("z [AU]")
        scatters = [ax.scatter(traj[0][i][0], traj[0][i][1], traj[0][i][2], s=2) for i in range(N)]
        if label=='auto':
            label = ['Body {}'.format(i) for i in range(1, N+1)]
            ax.legend(label, loc='center left', bbox_to_anchor=(1.2, 0.5))
        elif label == None:
            ax.legend_=None
        else:
            ax.legend(label, loc='center left', bbox_to_anchor=(1.2, 0.5))
        fig.set_tight_layout(True)
        if line:
            x = np.zeros((N, n_t_org))
            y = np.zeros((N, n_t_org))
            z = np.zeros((N, n_t_org))
            for i in range(N):
                for j in range(n_t_org):
                    x[i][j] = traj_org[j][i][0]
                    y[i][j] = traj_org[j][i][1]
                    z[i][j] = traj_org[j][i][2]
                ax.plot(x[i],y[i], z[i], linestyle='dotted', linewidth=0.3)
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.grid(False)
        ani = animation.FuncAnimation(fig, update_3, frames=n_t, repeat=False)
        ani.save('{}/ani_traj/{}_3D.gif'.format(directory,name), fps=30, dpi=dpi)
        if project:
            ax.clear()
            ax.remove()
            ax = fig.add_subplot(111)
            ax.set_title('N={}, '.format(N) + title)
            x_all = np.array([traj[i][:,0] for i in range(n_t)]).flatten()
            y_all = np.array([traj[i][:,1] for i in range(n_t)]).flatten()
            if lim:
                ax.set_xlim(lim[0])
                ax.set_ylim(lim[1])
            else:
                ax.set_xlim(x_all.min(), x_all.max())
                ax.set_ylim(y_all.min(), y_all.max())
            ax.set_xlabel("x [AU]")
            ax.set_ylabel("y [AU]")
            scatters = [ax.scatter(traj[0][i][0], traj[0][i][1], s=2) for i in range(N)]
            if label=='auto':
                label = ['Body {}'.format(i) for i in range(1, N+1)]
                ax.legend(label, loc='center left', bbox_to_anchor=(1.2, 0.5))
            elif label == None:
                ax.legend_=None
            else:
                ax.legend(label, loc='center left', bbox_to_anchor=(1.2, 0.5))
            if line:
                x = np.zeros((N, n_t_org))
                y = np.zeros((N, n_t_org))
                for i in range(N):
                    for j in range(n_t_org):
                        x[i][j] = traj_org[j][i][0]
                        y[i][j] = traj_org[j][i][1]
                for i in range(N):
                    ax.plot(x[i],y[i], linestyle='dotted', linewidth=0.3)
            ax.grid(False)
            fig.set_tight_layout(True)
            ani = animation.FuncAnimation(fig, update_2, frames=n_t, repeat=False)
            ani.save('{}/ani_traj/{}_projection.gif'.format(directory,name), fps=30, dpi=dpi)
    plt.show()

def plot(file, dim, title, lim=None,project=False, name='trajectories', label= 'auto',  dpi=300):
    """
    This function takes a txt file where the lines contains the timesteps and the positions of the N partilces at that timestep and plots the trajectories of the particles in 2D or 3D.
    
    Parameters
    -----------------------------
    file
        The txt file containing the times teps and positions
        The input needs to have the form: t  x_1  y_1  z_1 ... x_N  y_N  z_N
    dim
        Dimension of the Nbody simulation, has to be 2 or 3
    title
        Title of the plot
    lim
        Limit the axis
    project
        Whether to also plot a 2D projection
    name
        Name of the plot
    label
        Labels of trajectories. Defaults to automatic numbering
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
    plt.style.use('dark_background')
    if dim == 2:
        x = np.zeros((N, n_t))
        y = np.zeros((N, n_t))
        for i in range(N):
            for j in range(n_t):
                x[i][j] = traj[j][i][0]
                y[i][j] = traj[j][i][1]
        ax = fig.add_subplot(111)
        ax.set_title('N={}, '.format(N) + title)
        x_all = x.flatten()
        y_all = y.flatten()
        if lim:
            ax.set_xlim(lim[0])
            ax.set_ylim(lim[1])
        else:
            ax.set_xlim(x_all.min(), x_all.max())
            ax.set_ylim(y_all.min(), y_all.max())
        ax.set_xlabel("x [AU]")
        ax.set_ylabel("y [AU]")
        for i in range(N):
            ax.plot(x[i],y[i], linewidth=0.5)
        if label=='auto':
            label = ['Body {}'.format(i) for i in range(1, N+1)]
            ax.legend(label, loc='center left', bbox_to_anchor=(1.2, 0.5))
        elif label == None:
            ax.legend_=None
        else:
            ax.legend(label, loc='center left', bbox_to_anchor=(1.2, 0.5))
        fig.set_tight_layout(True)
        ax.grid(False)
        plt.savefig('{}/plot_traj/{}_2D.png'.format(directory,name), dpi=dpi, bbox_inches='tight')

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
        ax.set_title('N={}, '.format(N) + title)
        x_all = x.flatten()
        y_all = y.flatten()
        z_all = z.flatten()
        if lim:
            ax.set_xlim3d(lim[0])
            ax.set_ylim3d(lim[1])
            ax.set_zlim3d(lim[2])
        else:
            ax.set_xlim3d(x_all.min(), x_all.max())
            ax.set_ylim3d(y_all.min(), y_all.max())
            ax.set_zlim3d(z_all.min(), z_all.max())
        ax.set_xlabel("x [AU]")
        ax.set_ylabel("y [AU]")
        ax.set_zlabel("z [AU]")
        for i in range(N):
            ax.plot(x[i],y[i], z[i], linewidth=0.5)
        ax.grid(False)
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        if label=='auto':
            label = ['Body {}'.format(i) for i in range(1, N+1)]
            ax.legend(label, loc='center left', bbox_to_anchor=(1.2, 0.5))
        elif label == None:
            ax.legend_=None
        else:
            ax.legend(label, loc='center left', bbox_to_anchor=(1.2, 0.5))
        fig.set_tight_layout(True)
        plt.savefig('{}/plot_traj/{}_3D.png'.format(directory,name), dpi=dpi, bbox_inches='tight')
        
        if project:
            x = np.zeros((N, n_t))
            y = np.zeros((N, n_t))
            for i in range(N):
                for j in range(n_t):
                    x[i][j] = traj[j][i][0]
                    y[i][j] = traj[j][i][1]
            ax.clear()
            ax.remove()
            ax = fig.add_subplot(111)
            ax.set_title('N={}, '.format(N) + title)
            x_all = x.flatten()
            y_all = y.flatten()
            if lim:
                ax.set_xlim(lim[0])
                ax.set_ylim(lim[1])
            else:
                ax.set_xlim(x_all.min(), x_all.max())
                ax.set_ylim(y_all.min(), y_all.max())
            ax.set_xlabel("x [AU]")
            ax.set_ylabel("y [AU]")
            for i in range(N):
                ax.plot(x[i],y[i], linewidth=0.5)
            ax.grid(False)
            if label=='auto':
                label = ['Body {}'.format(i) for i in range(1, N+1)]
                ax.legend(label, loc='center left', bbox_to_anchor=(1.2, 0.5))
            elif label == None:
                ax.legend_=None
            else:
                ax.legend(label, loc='center left', bbox_to_anchor=(1.2, 0.5))
            fig.set_tight_layout(True)
            plt.savefig('{}/plot_traj/{}_projection.png'.format(directory,name), dpi=dpi, bbox_inches='tight')
    plt.show()

# Here the parameters can be changed
# ----------------------------------------------------------------------------
## file to animate/plot
trajectory = 'rings/rings_PEFRL_50.000000_0.001000.txt'
type_plot = 'plot'
## dimension of the file
dim = 3
## whether or not to plot 2D projection
project = True
## labels of the bodies. Options are: 'auto' which does automatic numbering, None which does not add labels, and a list ['body1', 'body2', ...] with the labels
label = 'auto'
## limit of the axis, default is None or thus automatic
lim = None
#lim = [[-10,10], [-10,10], [-10,10]]

# for animation
## plot trajectorires in animation
line = True
## timescale of the animation: how much should one second of animation be in simulation time
tstep = 5


# DO NOT CHANGE ANYTHING UNDER THIS
# ----------------------------------------------------------------------------
name = trajectory.rstrip('.txt')
args = name.split('_')
system = args[0].split('/')[0]
title = "{} integrated using {} \n tmax= {}, h= {}".format(system, args[1], args[2], args[3])
if len(args) == 5:
    title += "\n with adaptive timestep"

filename = system + '_' + args[1] 

if type_plot == 'animate':
    #tstep = int(input('How much should one second of video be equal to the time used in simulation: \n'))
    if line:
        filename += '_with_traj'
    print('Animating the trajectories...')
    animate(trajectory, dim, tstep, title=title, lim=lim, line= line, project= project, label=label, name= filename)
    
elif type_plot == 'plot':
    print('Plotting the trajectories...')
    plot(trajectory, dim, title=title,lim=lim , project=project, name= filename, label=label)
    
elif type_plot == 'both':
    #tstep = int(input('How much should one second of video be equal to the time used in simulation: \n'))
    print('Plotting the trajectories...')
    plot(trajectory, dim, title=title, lim=lim, project=project, name= filename, label=label)
    if line:
        filename += '_with_traj'
    print('Animating the trajectories...')
    animate(trajectory, dim, tstep, title=title, lim=lim, line=line, project= project, label=label, name= filename)