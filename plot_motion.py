import numpy as np
import matplotlib.pyplot as plt

times, x1, y1, x2, y2 = np.loadtxt("two_body_motion_RK2.txt", unpack=True)

plt.scatter(x1, y1, s = 0.01, label = 'Body 1')
plt.scatter(x2, y2, s = 0.01, label = 'Body 2')
plt.grid()
plt.legend()
plt.show()