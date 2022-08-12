import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

#複数の原子の塊を管理
class Atoms():
  def __init__(self):
    self.xyzs = np.empty((0,3))
    self.D = 4
  def add_atom(self, xyz):
    xyz = np.array([xyz])
    self.xyzs = np.append(self.xyzs,xyz,axis=0)
  def add_atoms(self, Atoms):
    self.xyzs = np.append(self.xyzs, Atoms.xyzs,axis=0)
  # 描画
  def show(self):
    x = self.xyzs[:,0]
    y = self.xyzs[:,1]
    z = self.xyzs[:,2]
    fig = plt.figure(figsize = (8, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, color = (1,1,0,0.7),s=0.1)
    max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() * 0.5
    mid_x = (x.max()+x.min()) * 0.5
    mid_y = (y.max()+y.min()) * 0.5
    mid_z = (z.max()+z.min()) * 0.5
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    ax.set_xlabel("Am")
    ax.set_ylabel("Am")
    ax.set_zlabel("Am")
    plt.show()

  def plane(self, xs, ys, zs):
    self.xyzs = self.make_plane(xs, ys, zs)
  
  def slide(self, x, y ,z):
    self.xyzs[:,0] += x
    self.xyzs[:,1] += y 
    self.xyzs[:,2] += z
  
  def rotateX(self, x):
    x = np.pi * x / 180
    M = np.array(
      [[1,0,0],
      [0,np.cos(x),-np.sin(x)],
      [0,np.sin(x),np.cos(x)]])
    self.xyzs = (M@self.xyzs.T).T

  def rotateY(self,y):
    y = np.pi * y / 180
    M = np.array(
      [[np.cos(y),0,np.sin(y)],
      [0,1,0],
      [-1*np.sin(y),0,np.cos(y)]])
    self.xyzs = (M@self.xyzs.T).T

  def rotateZ(self,z):
    z = np.pi * z / 180
    M = np.array(
      [[np.cos(z),-np.sin(z),0],
      [np.sin(z),np.cos(z),0],
      [0,0,1]])
    self.xyzs = (M@self.xyzs.T).T
  # 平面状の配列の生成(面心立方格子)

  def make_plane(self, xs, ys, zs):
    xyzs = np.zeros((4*xs*ys*zs,3))
    counter = 0
    for iz in range(zs):
      for iy in range(ys):
        for ix in range(xs):
          xyzs[counter] = [self.D * ix, self.D *iy , self.D * iz]
          xyzs[counter+1] = [self.D * (0.5 + ix), self.D * (0.5 + iy), self.D * iz]
          xyzs[counter+2] = [self.D * ix, self.D * (0.5 + iy), self.D * (0.5 + iz)]
          xyzs[counter+3] = [self.D * (0.5 + ix), self.D * iy, self.D * (0.5 + iz)]
          counter += 4
    xyzs[:,0] -= self.D*(xs-0.5)/2
    xyzs[:,1] -= self.D*(ys-0.5)/2
    xyzs[:,2] -= self.D*(zs-0.5)/2
    return xyzs
  #指定した条件を満たさない点を削除
  def trim_atoms(self,func):
    newxyzs = np.empty((0,3))
    for xyz in self.xyzs:
      x,y,z = xyz
      if(func(x,y,z)):
        newxyzs = np.append(newxyzs,[xyz],axis=0)
    self.xyzs = newxyzs




