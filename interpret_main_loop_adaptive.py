import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm
import seaborn as sns

"""

data = pd.read_csv("accuracy_vs_cost_adaptive_general_integrator.txt", sep=" ", header=None)
data = np.array(data)

data_points = (np.shape(data)[1]-2)//3
rows = np.shape(data)[0]


colors = plt.cm.tab20.colors

colors = sns.color_palette(n_colors=rows)


colors = sns.color_palette(n_colors=rows)

for row in range(rows):
    cost = [((data[row,2+3*j])) for j in range(data_points)]    
    accuracy = [abs((data[row,3+3*j])) for j in range(data_points)]
    plt.loglog(cost, accuracy, label = data[row,0], c = colors[row])

"""
data = pd.read_csv("accuracy_vs_cost_adaptive.txt", sep=" ", header=None)
data = np.array(data)

data_points = (np.shape(data)[1]-2)//3
rows = np.shape(data)[0]

colors = sns.color_palette('husl', n_colors=rows)
colors = sns.color_palette('Set3', n_colors=rows-2)
colors += ['cyan', 'magenta']
colors = plt.cm.tab20.colors
colors = sns.color_palette(n_colors=rows)

data_points = (np.shape(data)[1]-2)//3
rows = np.shape(data)[0]

#colors = sns.color_palette('husl', n_colors=rows)
#colors += ['cyan', 'magenta']

factor = 1
for row in [r for r in range(rows) if r not in [0, 11, 12, 1, 2, 4, 10]]:
    if row == 1:
        print(data[row],0)
        factor = -2
    elif row == 6:
        print(data[row],0)
        factor = -3
    elif row == 9:
        print(data[row],0)
        factor = -4
    else:
        factor = 1
    cost = [((data[row,2+3*j]))*factor for j in range(data_points)]    
    accuracy = [abs((data[row,3+3*j])) for j in range(data_points)]
    plt.loglog(cost, accuracy, label = data[row,0])#, c = colors[row])
plt.xlabel('number of driver evaluations', fontsize=15)
plt.ylabel('accuracy', fontsize=15)
plt.legend()
plt.title('Lemniscate - adaptive timestep - high orders', fontsize=15)
plt.show()



"""
cost = [(abs(data[:,2+3*j])) for j in range(data_points)]    
accuracy = [(abs(data[:,3+3*j])) for j in range(data_points)]
plt.loglog(cost, accuracy, label = data[:,0])
plt.xlabel('number of driver evaluations')
plt.ylabel('accuracy')
plt.legend()
plt.show()
"""