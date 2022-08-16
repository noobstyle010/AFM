import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('./data/xyz.csv',names=['x', 'y', 'z'],sep="\s")
fig = plt.figure(figsize=plt.figaspect(1))
ax = fig.add_subplot(projection='3d')

X = data['x'].values
Y = data['y'].values
Z = data['z'].values

max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() * 0.5

mid_x = (X.max()+X.min()) * 0.5
mid_y = (Y.max()+Y.min()) * 0.5
mid_z = (Z.max()+Z.min()) * 0.5
ax.set_xlim(mid_x - max_range, mid_x + max_range)
ax.set_ylim(mid_y - max_range, mid_y + max_range)
ax.set_zlim(mid_z - max_range, mid_z + max_range)


ax.scatter(X, Y, Z,c=Z,s=1,cmap="jet")
plt.savefig("img/point_cloud.jpg")