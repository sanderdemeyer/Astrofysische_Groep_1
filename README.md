# Astrophysical Simulations Group 1
Code for the N-body simulation.

## **`N_body_sim.cpp`**
This is the final code for the N-body simulation. The units used in the simulations are:

* length: [AU]
* time: [yr]
* mass: $[M_{\odot}]$
* G = 39.473107 $[AU]^3 M_{\odot}^{-1} [yr]^{-2}$

The initial conditions are read from the **`Initial_conditions`** folder. The filename should not contain any underscores as they are used as separators. The initial condition files have the given format:

$$
\begin{gather*}
m_1 \quad x_1 \quad y_1 \quad z_1 \quad v_1 \quad v_1 \quad v_1 \\
m_2 \quad x_2 \quad y_2 \quad z_2 \quad v_2 \quad v_2 \quad v_2\\
\vdots\\
m_n \quad x_n \quad y_n \quad z_n \quad v_n \quad v_n \quad v_n
\end{gather*}
$$

The trajectories and energies for are stored in folders with the names of the corresponding systems in the **`traj`** en **`energy`** folders respectively. The naming convention used in each subfolder is:
``SystemName_integrator_tmax_h(_adaptive).txt``. The energy files contain two columns with the time and energy. The format of the trajectory files is:

$$
\begin{gather*}
t_0 \quad x_{1, t_0} \quad y_{1, t_0} \quad z_{1, t_0} \quad x_{2, t_0} \quad y_{2, t_0} \quad z_{2, t_0} \quad \cdots \quad x_{n, t_0} \quad y_{n, t_0} \quad z_{n, t_0}\\
t_1 \quad x_{1, t_1} \quad y_{1, t_1} \quad z_{1, t_1} \quad x_{2, t_1} \quad y_{2, t_1} \quad z_{2, t_1} \quad \cdots \quad x_{n, t_1} \quad y_{n, t_1} \quad z_{n, t_1}\\
\vdots \\
t_{m} \quad x_{1, t_{m}} \quad y_{1, t_{m}} \quad z_{1, t_{m}} \quad x_{2, t_{m}} \quad y_{2, t_{m}} \quad z_{2, t_{m}} \quad \cdots \quad x_{n, t_{m}} \quad y_{n, t_{m}} \quad z_{n, t_{m}}\\
\end{gather*}
$$

Two types of integrators have been defined
1. ```friend```: These integrators are defined using recursive formula's and are defined as a ``friend void`` class for the ```NSystem``` class in **`classes.cpp`**. This is the default integrator type and is used in the ```integrate``` function (see further).
2. ```general```: These are integrators defined as a new class: ```General_integrator```. These integrators are defined using Butcher tableau's and used in the ```integrate_general``` function (see further).

Some functions are defined to run the simulation:
* ```integrate```: This function takes some given N-body initial conditions and calculates their trajectories and total energy for a given maximum time using a given integrator and initial timestep. All integrators can be used with an adaptive timestep and if ADAPTIVE_RK45 is true then RK45 will be used. The available integrators are:
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
    - RK45: called by a boolean, any of the above given integrators can be given for the ```integrator``` parameter as the function won't work without this!
* ```integrate_general```: This function takes some given N-body initial conditions and calculates their trajectories and total energy for a given maximum time using a given integrator and initial timestep. All integrators can be used with an adaptive timestep. This is essentially the same function as `integrate` with the only difference being that here the integrators are defined using Butcher tableau's. The available integrators are:
    - 3 over 8
    - Ralston4
    - RK4
    - RK5
    - RK6
    - RK8
    - SSPRK3
    - Wray3
* ```loop_h```: This function takes some given N-body initial conditions and runs the `integrate` function for fixed timestep using different timesteps in the given range.
* ```loop_h_general```: This function takes some given N-body initial conditions and runs the `integrate_general` function for a number of fixed timesteps in the given range.


To use the code, the initial conditions, integrator type, integrator, (initial) timestep and the maximum time to simulate should be given. Along with whether to use adaptive timestep and whether to use RK45 integrator which has embedded adaptive timestep. All integrators can be used with adaptive timestep. To use RK45 the boolean ```ADAPTIVE_RK45``` should be set to true. To use ```ADAPTIVE_RK45``` the integrator type needs to be 'friend'. An example to integrate the Burrau initial conditions for 70 years using an initial timestep of 0.001 year using the RK4 integrator of the type ```General_integrator``` with adaptive timestep is:

```cpp
// File with the initial conditions to be read from the `Initial_conditions` folder
std::string in_cond = "Burrau.txt";
 // The type of integrator to be used. By default the `friend void` type integrators are used.
std::string type_integ = "general";
// he type of integrator to be used. For each type, available integrators are listed in the README file.
std::string integrator = "RK4";
// The (initial) timestep.
double h = 0.001;
 // Total time considered in the simulation
double tmax = 70;
// Whether or not to use RK45 embedded adaptive timestep. Only implemented in the `integrate` function.
bool ADAPTIVE_RK45 = false;
// Whether or not to use an adaptive timestep
bool ADAPTIVE_TIME_STEP = true;
// The following should only be changed to loop over run the simulation for different timesteps.
bool h_loop = false;
// the maximum timestep
double hmax = 0.1;
// the factor by which to loop from `h` to `hmax`. Should be greater than 1.
int step = 10;
```
Only these need to be changed in the code to run other simulations. To loop over different timesteps ```h_loop``` should be set to true and ```hmax``` and ```step``` should be given.

## Regularization

The files **`main_2D.cpp`** and **`classes_regularization.cpp`** can be used to simulate 2D systems with regularization. By default, z = 0.

The parameters to be used are the same as for **`main.cpp`**, with the exception of: 

    - type_integ: this is always "function"
    - ADAPTIVE_RK45: this is always false.
    - h_loop: This is always false.
    - transform_distance: this is an extra parameter that denotes below which distance 2 bodies are regularized.

Disclaimers:
    - The code is always 2-dimensional, with z = 0 and v_z = 0.
    - The system loses energy upon each regularization. This means that the code is not fully correct.
    - when transform_distance = 0.0, this is equivalent to  **`N_body_sim.cpp`** and works in 3D too.
    - The trajectories and energies are saved in the folders **`traj_reg`** and **`energy_reg`** respectively.
    - Since the number of bodies can change, the trajectory file does not work with **`plot_traj`**, which assumes a fixed number of bodies. Thus, alternatively, the files **`plot_motion.py`** and **`plot_energy.py`** should be used.


## Contents of the different folders
* **`ani_traj`**: The animated trajectories for some chosen initial conditions. 
    * Naming convention used: ``SystemName_integrator_(with_traj_)(3D/projection).gif``

* **`dist_mercury`**: Distance and angle of all the bodies in the system with respect to the second body. This is used to study the Solar System and isn't valid for other initial conditions. 

* **`energy`**: Text files with the total energy at each step for given initial conditions and integrator.
    * Naming convention used: ``SystemName/SystemName_integrator_tmax_h(_adaptive).txt``

* **`energy_reg`**: Text files with total energy at each step for simulation with regularization.

* **`Initial_conditions`**: Text files with some interesting initial conditions. The filename should not contain any underscores as they are used as separators. The files have the given format:

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
    * Naming convention used: ``SystemName_integrator_(3D/projection).png``

* **`random_gauss`**: 2- to 99-body initial conditions according to a Gaussian distribution, generated using **`random_gauss.py`**. These are used to study the scaling behavior of the different integrators used.

* **`scaling_friend`**: Text files with the execution time in milliseconds for 1000 steps for 2- to 99-body systems for all the ```void friend``` class integrators (see **`main.cpp`**) with a timestep of 0.01 y.

* **`scaling_general`**: Text files with the execution time in milliseconds for 1000 steps for 2- to 99-body systems for all the ```General_integrator``` class integrators (see **`main.cpp`**) with a timestep of 0.01 y.

* **`temperatures`**: This is only relevant for systems with three stars and a planet. This folder contains text files with the temperature of the planet.

* **`traj`**: Text files containing the calculated trajectories for a given system with a given integrator, timestep and simulation time. 
  * Naming convention used: ``SystemName/SystemName_integrator_tmax_h(_adaptive).txt``
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