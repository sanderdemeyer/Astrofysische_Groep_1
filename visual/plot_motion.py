### LEGACY ###

import numpy as np
import matplotlib.pyplot as plt
#times, x1, y1, x2, y2 = np.loadtxt("two_body_motion_RK2.txt", unpack=True)
#times, x1, y1, z1, x2, y2, z2, x3, y3, z3 = np.loadtxt("three_body_motion_RK4.txt", unpack=True)
#data = np.loadtxt("two_body_motion_RK4.txt", unpack=True)



#data = np.loadtxt("traj/Solar-System_RK4_2500000_0.001000_adaptive.txt", unpack=True)
#data = np.loadtxt("traj/3_body_problem_RK4_500000_0.001000_adaptive.txt", unpack=True)
#data = np.loadtxt("traj/3_body_RK4.txt", unpack=True)
#data = np.loadtxt("two_body_motion_RK4.txt", unpack = True)
data = np.loadtxt("traj_reg/Reg_2D_two-body_RK4_2500000_0.001000.txt", unpack = True)
data = np.loadtxt("traj_reg/Reg_2D_two-body_RK4_2500000_0.001000.txt", unpack = True)
data = np.loadtxt("traj_reg/Reg_2D_Burrau_RK4_2500000_0.001000.txt", unpack = True)
data = np.loadtxt("traj/lemniscate_RK4_1000000_0.001000_adaptive.txt", unpack=True)
#data = np.loadtxt("traj/two-body_RK4_2500000_0.001000_adaptive.txt", unpack = True)
data = np.array(data)

print(np.shape(data))

times = data[:,0]
columns = np.shape(data)[0]
number_of_bodies = (columns-1)//3

"""

(r, c) = np.shape(data)
print(r)
print(c)
sampling = 10000
data_new = np.zeros((r, c//sampling))

for i in range(r):
    for j in range(c//sampling):
        data_new[i,j] = data[i,sampling*j]

print(np.shape(data_new))
"""

for body in range(number_of_bodies):
    plt.scatter(data[1+3*body], data[2+3*body], s = 0.01, c = range(np.shape(data)[1]), label = f'Body {body+1}')
    #plt.scatter(data_new[1+3*body], data_new[2+3*body], s = 0.01, c = range(np.shape(data_new)[1]), label = f'Body {body+1}')
#plt.scatter(x1, y1, s = 0.01, label = 'Body 1')
#plt.scatter(x2, y2, s = 0.01, label = 'Body 2')
#plt.scatter(x3, y3, s = 0.01, label = 'Body 3')
plt.grid()
plt.legend()
plt.show()