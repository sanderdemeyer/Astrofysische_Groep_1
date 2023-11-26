# Astrofysische_Groep_1
Code voor het project van Astrofysische Simulaties - Groep 1

De oplossingen van de verschillende integrators worden opgeslagen in **`traj`** en **`energy`**.

De folder **`visual`** bevat de python bestanden gebruikt voor visualisatie. De animaties en plots worden opgeslagen in **`ani_traj`**, **`plot_traj`** en **`plot_energy`**.

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
*SystemName_integrator_iterations_h(_adaptive).txt*
