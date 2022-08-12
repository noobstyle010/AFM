import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('./data/xyz.csv',names=['x', 'y', 'z'],sep="\s")
fig = plt.figure()


x = data['x'].values
y = data['y'].values
z = data['z'].values

ax = fig.add_subplot(projection='3d')
ax.scatter(x, y, z)
#plt.show()
plt.savefig("img/point_cloud.jpg")