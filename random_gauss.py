import numpy as np

for i in range(2, 100):
    pos = np.random.randn(i, 3)
    vel = np.random.randn(i, 3)

    with open('random_gauss/{:02d}.txt'.format(i), 'w') as d:
            for t in range(i):
                d.write('1 {} {} {} {} {} {}\n'.format(pos[t][0], pos[t][1], pos[t][2], vel[t][0], vel[t][1], vel[t][2]))