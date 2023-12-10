import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.cm as cm
import seaborn as sns

data = pd.read_csv("accuracy_vs_cost_non_adaptive_more_data.txt", sep=" ", header=None)
data = pd.read_csv("accuracy_vs_cost_adaptive_more_data.txt", sep=" ", header=None)
data = np.array(data)

data_points = (np.shape(data)[1]-2)//3
rows = np.shape(data)[0]

colors = sns.color_palette('husl', n_colors=rows)

colors = sns.color_palette('Set3', n_colors=rows-2)

colors += ['cyan', 'magenta']

colors = plt.cm.tab20.colors

colors = sns.color_palette(n_colors=rows)

for row in range(rows):
    cost = [(abs(data[row,2+3*j])) for j in range(data_points)]    
    accuracy = [(abs(data[row,3+3*j])) for j in range(data_points)]
    plt.loglog(cost, accuracy, label = data[row,0], c = colors[row])

"""
data = pd.read_csv("accuracy_vs_cost_non_adaptive_more_data_general_integrator.txt", sep=" ", header=None)
data = np.array(data)

data_points = (np.shape(data)[1]-2)//3
rows = np.shape(data)[0]

for row in range(rows):
    cost = [(abs(data[row,2+3*j])) for j in range(data_points)]    
    accuracy = [(abs(data[row,3+3*j])) for j in range(data_points)]
    plt.loglog(cost, accuracy, label = data[row,0], c = colors[row])
"""
    
plt.xlabel('number of driver evaluations')
plt.ylabel('accuracy')
plt.legend()
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