from sys import stdin
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
from scipy.interpolate import griddata

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
def zscore(x, axis = None):
  xmean = x.mean(axis=axis, keepdims=True)
  xstd  = np.std(x, axis=axis, keepdims=True)
  zscore = (x-xmean)/xstd
  return zscore

# ./a.out >> python3 run.py >> test.txt こんな感じでやれば良さげ
# gifイメージの作成用
input = stdin.readline

###  具体的な処理　###
sn,pn = map(int, input().split())
nx,ny,nz = map(int, input().split())
dx,dy,dz = map(float, input().split())


fps = 10
frn = nz

F = np.ones((nx,ny,nz))
ims = []

# 表面
for i in range(sn):
  x_,y_,z_ = map(float, input().split())
# プローブ
for i in range(pn):
  x_,y_,z_ = map(float, input().split())

# 力場の変化のanimation
num = 0
while num < nz:
  # Forceのplot
  for iy in range(ny):
    for ix in range(nx):
      x_,y_,z_,f_ =  map(float, input().split())
      F[ix][iy][num] = f_
  num+=1

fig = plt.figure()
ax = fig.add_subplot(111)
div = make_axes_locatable(ax)
cax = div.append_axes('right', '5%', '5%')

def plot(data):
  cax.cla()
  ax.set_xlabel('X[\N{ANGSTROM SIGN}m]')
  ax.set_ylabel('Y[\N{ANGSTROM SIGN}m]')
  f = F[:,:,data]
  norm = colors.DivergingNorm(vmin=-1*abs(max(f.min(),f.max(),key=abs)),vcenter=0., vmax=abs(max(f.min(),f.max(),key=abs)))
  ax.set_title('Atomic Force at {:.1f} [\N{ANGSTROM SIGN}m]'.format(2+data*dz))
  im = ax.imshow(f,extent=[0,nx*dx,0,ny*dy],origin='lower',norm=norm,cmap='bwr')
  fig.colorbar(im,cax=cax)


ani = animation.FuncAnimation(fig, plot, frames= nz,interval=500, repeat=False)
fn = 'plot_surface_animation_funcanimation'
ani.save("anim.mp4", writer="ffmpeg")


