import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

data = pd.read_csv('./data/scanXY.csv',header=None,sep="\s")
fig = plt.figure()
print(data.values)
plt.imshow(data.values)
plt.colorbar()
plt.show()