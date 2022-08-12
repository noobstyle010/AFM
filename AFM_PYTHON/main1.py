import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize
import matplotlib.colors as mcolors
from AtomicForce import Atoms
import sys
import os
'''
分子修飾の影響をシミュレーションする
'''
args = sys.argv
needle = int(args[1])
probe_r = int(args[2])
DD = float(args[3])
needle_to_needle = float(args[4])
needle_l = int(args[5])
rotates = float(args[6])

xy_path = './imgs/' + args[1] + '-' + args[2] + '-' + args[3] + '-' + args[4] + '-' + args[5] + '-' + args[6] + '.jpg'
setup_path = './setups/' + args[1] + '-' + args[2] + '-' + args[3] + '-' + args[4] + '-' + args[5] + '-' + args[6] + '.jpg'

#単位は無次元化する
eV = 1.609e-19
k = 1.38062e-23
m0 = 1e-26
d = 1e-10

# z方向のみ
def LennardJonesForce(R,epsilon=0.2294,sigma = 2.951):
  r = np.linalg.norm(R)
  F = -48*epsilon*((sigma/r)**12 - 0.5 * (sigma/r)**6)*R[2]/(r**2)
  F = F * d / eV
  return F

# AFMのシミュレータクラス
class AFM():
  probe_atoms = 0
  surface_atoms = 0
  def __init__(self, probe, surface, DD):
    self.probe_atoms = probe.xyzs
    self.surface_atoms = surface.xyzs
    self.DD = DD
    probe_min_z = np.min(probe.xyzs[:,2])
    surface_max_z = np.max(surface.xyzs[:,2])  
    self.surface_atoms[:,2] -= surface_max_z
    self.probe_atoms[:,2] += DD-probe_min_z
    # print("probeのZで小さいもの")
    # print(np.sort(probe.xyzs[:,2])[:3])

  def show_set_up(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x1=  self.probe_atoms[:,0]
        y1 = self.probe_atoms[:,1]
        z1 = self.probe_atoms[:,2]
        ax.scatter(x1,y1,z1,color=(1,0,0,0.7),s=3)
        x2 =  self.surface_atoms[:,0]
        y2 = self.surface_atoms[:,1]
        z2 = self.surface_atoms[:,2]
        ax.scatter(x2,y2,z2,color=(0,0,1,0.7),s=0.3)
        x = np.concatenate([x1, x2])
        y = np.concatenate([y1, y2])
        z = np.concatenate([z1, z2])
        max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()]).max() * 0.5
        mid_x = (x.max()+x.min()) * 0.5
        mid_y = (y.max()+y.min()) * 0.5
        mid_z = (z.max()+z.min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)
        plt.savefig(setup_path)
  
  def scanZ(self):
    #probe_atomsを動かす
    step = 0.1
    z = []
    f = []
    for i in range(100):
      z.append(step*i)
      f.append(self.calculate_attractive_force())
      self.probe_atoms[:,2] += step
    plt.plot(z, f)
    plt.show()

  def scanXY(self,nx, ny, step):
    plt.clf()
    x = np.arange(nx) * step
    y = np.arange(ny) * step
    f = np.zeros([len(y),len(x)])
    for i in range(nx):
      self.probe_atoms[:,0] += step
      for j in range(ny):
        self.probe_atoms[:,1] += step
        f[j][i] = self.calculate_attractive_force()
      self.probe_atoms[:,1] -= step * ny
    norm = mcolors.TwoSlopeNorm(vmin=-1*abs(max(f.min(),f.max(),key=abs)),vcenter=0., vmax=abs(max(f.min(),f.max(),key=abs)))
    plt.imshow(f,extent=[0,nx*step,0,ny*step],origin='lower',norm=norm,cmap='bwr')
    
    plt.xlabel('X[\N{ANGSTROM SIGN}m]')
    plt.ylabel('Y[\N{ANGSTROM SIGN}m]')
    plt.title('Atomic Force at {} [\N{ANGSTROM SIGN}m]'.format(self.DD))
    plt.colorbar()
    plt.savefig(xy_path)

  def calculate_attractive_force(self):
    #行列計算による高速化？してるか　要計測
    lines = self.surface_atoms.shape[0]
    linep = self.probe_atoms.shape[0]
    new_probe_atoms = np.tile(np.ravel(self.probe_atoms),(lines,1))
    new_surface_atoms = np.tile(self.surface_atoms,(1,linep))
    Rs = new_probe_atoms - new_surface_atoms
    Rs = np.reshape(Rs,[-1,3])
    #距離のベクトルを出したのでファンデルワールス力を求めていく
    attractive_force = 0
    for R in Rs:
      attractive_force += LennardJonesForce(R)
    return attractive_force


base_surface = Atoms()
base_probe = Atoms()

### プローブの準備
g_p = Atoms()
g_p.plane(int(2*probe_r/2.7)+6,int(2*probe_r/2.7)+6,int(2*probe_r/2.7)+10)
g_p.trim_atoms(lambda x,y,z: x**2 + y**2 + z**2 <= probe_r*probe_r and z<=-1*probe_r+15)
base_probe.add_atoms(g_p)

if needle:
  ### 針をつける
  #１層目
  for i in range(needle_l):
    base_probe.add_atom([0,0,-(probe_r+1+i)])
  #二層目
  sita = needle_to_needle / probe_r 
  tmp = needle_to_needle / (probe_r*np.sin(sita))
  fain = int(2*np.pi/tmp)
  for fai in range(fain):
    #åprint(fai)
    fai = fai * tmp
    for i in range(needle_l):
      base_probe.add_atom([(probe_r+1+i)*np.sin(sita)*np.cos(fai),(probe_r+1+i)*np.sin(sita)*np.sin(fai),-(probe_r+1+i)*np.cos(sita)])
base_probe.rotateX(rotates)


base_surface.plane(3,3,1)
base_surface.trim_atoms(lambda x,y,z: x>=0 and x<8 and y>=0 and y<=8)

#距離の調整
AFM = AFM(base_probe, base_surface, DD)
AFM.show_set_up()
AFM.scanXY(100,100,0.1)
