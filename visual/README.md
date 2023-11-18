## **` plot_traj.py `**
This program plots the calculated trajectories of the N-body simulation. The file with trajectories is read from the **`traj`** folder in the parent directory.

### Input
```Please provide the file:```  
Only the filename should be provided given the trajectories are stored in the **`traj`** folder.


```Would you like to animate or plot the trajectories:```  
Answer with ```animate``` or ```plot```.


```What is the dimension of the simulation:```  
Answer with ```2``` or ```3```.


```Would you like to limit the axis of the simulation?```  
To set the limits answer with ```y``` and then provide the limits in the form ```lower_limit, upper_limit```. Otherwise, press any key.

### Output
The animations are stored in the **`ani_traj`** folder and the plots are stored in the **`plot_traj`** folder.

## **` plot_energy.py `**
This program plots the total energy and the relative energy error of the N-body simulation. The file with energy is read from the **`energy`** folder in the parent directory.

### Input
```Please provide the file:```  
Only the filename should be provided given the trajectories are stored in the **`energy`** folder.


```Would you like to plot the total energy or the relative error:```  
Answer with ```total energy```, ```relative error``` or ```both```.

### Output
The plots are stored in the **`plot_energy`** folder.