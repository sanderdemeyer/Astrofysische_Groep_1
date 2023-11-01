import numpy as np
import matplotlib.pyplot as plt

# times, x, y = np.loadtxt("Practice Session 3/two_body_motion.txt", unpack=True)
times, x1, y1, x2, y2 = np.loadtxt("Astrofysische_Groep_1/two_body_motion_RK2.txt", unpack=True)


# plt.scatter(x1, y1, s = 0.01, c = np.array(range(len(x1)))/len(x1), label = 'Body 1')
# plt.scatter(x2, y2, s = 0.01, c = np.array(range(len(x2)))/len(x2), label = 'Body 2')
plt.scatter(x1, y1, s = 0.01, label = 'Body 1')
plt.scatter(x2, y2, s = 0.01, label = 'Body 2')
plt.legend()
plt.show()