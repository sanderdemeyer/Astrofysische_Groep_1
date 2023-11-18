import numpy as np
import matplotlib.pyplot as plt

times, x1, y1, x2, y2 = np.loadtxt("two_body_motion_RK2.txt", unpack=True)
#times, x1, y1, z1, x2, y2, z2, x3, y3, z3 = np.loadtxt("three_body_motion_RK4.txt", unpack=True)
data = np.loadtxt("three_body_motion_RK4.txt", unpack=True)
times = data[:,0]
columns = np.shape(data)[0]

number_of_bodies = (columns-1)//3

for body in range(number_of_bodies):
    plt.scatter(data[1+3*body], data[2+3*body], s = 0.01, label = f'Body {body+1}')
#plt.scatter(x1, y1, s = 0.01, label = 'Body 1')
#plt.scatter(x2, y2, s = 0.01, label = 'Body 2')
#plt.scatter(x3, y3, s = 0.01, label = 'Body 3')
plt.grid()
plt.legend()
plt.show()