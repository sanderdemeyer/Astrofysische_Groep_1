# Astrophysical Simulations Group 1
Code for the N-body simulation.

The trajectories and energies for different integrators are stored in **`traj`** en **`energy`** folders respectively.

## **`main.cpp`**
Final README for main comes here.

## Contents of the different folders
* **`ani_traj`**: The animated trajectories for some chosen initial conditions. 
    * Naming convention used: *SystemName_integrator_(with_traj_)(3D/projection).mkv*

* **`energy`**: Text files with the total energy at each step for given initial conditions and integrator.
    * Naming convention used: *SystemName_integrator_tmax_h(_adaptive).txt*

* **`energy_reg`**: Text files with total energy at each step for simulation with regularization.

* **`initial_conditions`**: Text files with some interesting initial conditions. The filename should not contain any underscores as they are used as separators. The files have the given format:
$$
\begin{gather*}
m_1 \quad x_1 \quad y_1 \quad z_1 \quad v_1 \quad v_1 \quad v_1 \\
m_2 \quad x_2 \quad y_2 \quad z_2 \quad v_2 \quad v_2 \quad v_2\\
\vdots\\
m_n \quad x_n \quad y_n \quad z_n \quad v_n \quad v_n \quad v_n
\end{gather*}
$$

* **`legacy`**: Old discarded code.

* **`performance`**: Plots used to compare the different integrators.

* **`plot_energy`**: Plots of the relative energy error in function of the simulated time for some chosen initial conditions and integrators.

* **`plot_traj`**: Plotted trajectories of some chosen for initial conditions with a given integrator. The timestep used and the simulated time can be found in the titles of the plots.
    * Naming convention used: *SystemName_integrator_(3D/projection).png*

* **`random_gauss`**: 2- to 99-body initial conditions according to a Gaussian distribution, generated using **`random_gauss.py`**. These are used to study the scaling behavior of the different integrators used.

* **`scaling_friend`**: Text files with the execution time in milliseconds for 1000 steps for 2- to 99-body systems for all the ```void friend``` class integrators (see **`main.cpp`**) with a timestep of 0.01 y.

* **`scaling_general`**: Text files with the execution time in milliseconds for 1000 steps for 2- to 99-body systems for all the ```General_integrator``` class integrators (see **`main.cpp`**) with a timestep of 0.01 y.

* **`temperatures`**: IDK

* **`traj`**: Text files containing the calculated trajectories for a given system with a given integrator, timestep and simulation time. 
  * Naming convention used: *SystemName_integrator_tmax_h(_adaptive).txt*
  * The files have the given format:
$$
\begin{gather*}
t_0 \quad x_{1, t_0} \quad y_{1, t_0} \quad z_{1, t_0} \quad x_{2, t_0} \quad y_{2, t_0} \quad z_{2, t_0} \quad \cdots \quad x_{n, t_0} \quad y_{n, t_0} \quad z_{n, t_0}\\
t_1 \quad x_{1, t_1} \quad y_{1, t_1} \quad z_{1, t_1} \quad x_{2, t_1} \quad y_{2, t_1} \quad z_{2, t_1} \quad \cdots \quad x_{n, t_1} \quad y_{n, t_1} \quad z_{n, t_1}\\
\vdots \\
t_{m} \quad x_{1, t_{m}} \quad y_{1, t_{m}} \quad z_{1, t_{m}} \quad x_{2, t_{m}} \quad y_{2, t_{m}} \quad z_{2, t_{m}} \quad \cdots \quad x_{n, t_{m}} \quad y_{n, t_{m}} \quad z_{n, t_{m}}\\
\end{gather*}
$$

* **`traj_reg`**: Text files containing the calculated trajectories for a given system with a given integrator, timestep and simulation time with regularization.

* **`visual`**: Python code used for the visualization of the N-body trajectories.

## **`main.cpp`**
### Input
```Provide the initial conditions:```  
The initial conditions are stored in the **`Initial_conditions`** folder. Only the filename should be provided. The txt file should have the given format:

$$
\begin{gather*}
m_1 \quad x_1 \quad y_1 \quad z_1 \quad v_1 \quad v_1 \quad v_1 \\
m_2 \quad x_2 \quad y_2 \quad z_2 \quad v_2 \quad v_2 \quad v_2\\
\vdots\\
m_n \quad x_n \quad y_n \quad z_n \quad v_n \quad v_n \quad v_n
\end{gather*}
$$


```Which integrator would you like to use:```  
Answer with one of the following:
- Forward Euler
- Heun
- Heun3
- Ralston
- Ralston3
- RK3
- RK4
- Forest Ruth
- PEFRL
- Velocity Verlet
- Position Verlet
- Leapfrog
- Yoshida_4

```Provide the initial timestep: ```  
Give the initial timestep.

```Provide the number of iterations: ```  
Give the number of iterations.

### Output
The trajectories are stored in the **`traj`** folder and the energies are stored in the **`energy`** folder. The naming convention is:
*SystemName_integrator_tmax_h(_adaptive).txt*

## **`scaling_friend.cpp`**
Calculates the execution time in milliseconds for 1000 steps for 2- to 99-body systems for all the ```void friend``` class integrators with a timestep of 0.01 y. The integrators studied are:
- Forest Ruth
- Forward Euler
- Heun
- Heun3
- Leapfrog
- Position Verlet
- PEFRL
- Ralston
- Ralston3
- RK2
- RK3
- RK4
- Velocity Verlet
- Yoshida 4

## **`scaling_general.cpp`**
Calculates the execution time in milliseconds for 1000 steps for 2- to 99-body systems for all the ```General_integrator``` class integrators with a timestep of 0.01 y. The integrators studied are:
- 3 over 8
- Ralston4
- RK4
- RK5
- RK6
- RK8
- SSPRK3
- Wray3