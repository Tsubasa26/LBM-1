import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import shift
# from mpi4py import MPI


# def streaming(f):
#   # e:0,0 
#   fsB[:,:,0] = f[:,:,0]
#   # e:1,0 
#   fsB[:,:,1] = shift(f[:,:,1], [1,  0], cval=0.0)
#   fsB[:,:,2] = shift(f[:,:,2], [0,  1], cval=0.0)
#   fsB[:,:,3] = shift(f[:,:,3], [-1, 0], cval=0.0)
#   fsB[:,:,4] = shift(f[:,:,4], [0, -1], cval=0.0)
#   # e:1,1
#   fsB[:,:,5] = shift(f[:,:,5], [1,  1], cval=0.0)
#   fsB[:,:,6] = shift(f[:,:,6], [-1, 1], cval=0.0)
#   fsB[:,:,7] = shift(f[:,:,7], [-1,-1], cval=0.0)
#   fsB[:,:,8] = shift(f[:,:,8], [1, -1], cval=0.0)
#   return fsB
def streaming(nx,ny,f,e):
  # input  : f,e,nx,ny
  # output : fsB
  for i in range(nx):
    for j in range(ny):
      for l in range(9):
        if((i+e[l,0]) < 0 or (j+e[l,1]) < 0 or (i+e[l,0]) > nx-1 or (j+e[l,1]) > ny-1):
            pass
        else:
            fsB[i+e[l,0],j+e[l,1],l] = f[i,j,l] 
  return fsB


def setBounceBack(f, rho, u, v):
  # # 6, 2, 5
  # # 3, 0, 1
  # # 7, 4, 8
  # Left  6,3,7
  f[0,:,1] = f[0,:,3]
  f[0,:,5] = f[0,:,7]
  f[0,:,8] = f[0,:,6]
  # Right 5,1,8
  f[nx-1,:,3] = f[nx-1,:,1]
  f[nx-1,:,6] = f[nx-1,:,8]
  f[nx-1,:,7] = f[nx-1,:,5]
  # Bottom 7,4,8
  f[:,0,2] = f[:,0,4]
  f[:,0,5] = f[:,0,7]
  f[:,0,6] = f[:,0,8]
  # Top (Zoe-He)
  f[:,ny-1, 4] = f[:, ny-1, 2]
  f[:,ny-1, 7] = f[:, ny-1, 5] \
                 + 0.50 * (f[:, ny-1, 1] - f[:, ny-1, 3]) \
                 - 0.50 * rho[:, ny-1] * u[:, ny-1]
  f[:,ny-1, 8] = f[:, ny-1, 6] \
                 - 0.50 * (f[:, ny-1, 1] - f[:, ny-1, 3]) \
                 + 0.50 * rho[:, ny-1] * u[:, ny-1]
  # reset
  # Left  6,3,7
  f[0,:,3] = 0.0
  f[0,:,6] = 0.0
  f[0,:,7] = 0.0
  # Right 5,1,8
  f[nx-1,:,1] = 0.0
  f[nx-1,:,5] = 0.0
  f[nx-1,:,8] = 0.0
  # Bottom 7,4,8
  f[:,0,4] = 0.0
  f[:,0,7] = 0.0
  f[:,0,8] = 0.0
  # Top 6,2,5
  f[:,ny-1,2] = 0.0
  f[:,ny-1,5] = 0.0
  f[:,ny-1,6] = 0.0
  return f


def calc_rho(f):
  rho = f.sum(axis=2)
  return rho

def calc_vel(f, e, rho):
  u = (f * e[:,0]).sum(axis=2) / rho
  v = (f * e[:,1]).sum(axis=2) / rho
  return u, v

def calc_feq(nx, ny, w, rho, u, v):
  feq = np.zeros([nx,ny,9])
  for i in range(nx):
    for j in range(ny):
      for l in range(9):
        eu = e[l,0] * u[i,j] + e[l,1] * v[i,j]
        uu = u[i,j]**2 + v[i,j]**2
        s = w[l] * (3.0 * eu + 4.5 * eu**2 - 1.5 * uu)
        feq[i,j,l] = w[l] * rho[i,j] + rho[i,j] * s
  return feq

def setboundary_vel(nx,ny,u,v,rho,fsB):
  # bottom
  u[:,0] = 0.0
  v[:,0] = 0.0
  rho[:,0] = rho[:,1]

  # Left
  u[0,:] = 0.0
  v[0,:] = 0.0
  rho[0,:] = rho[1,:]

  # Right
  u[nx-1,:] = 0.0
  v[nx-1,:] = 0.0
  rho[nx-1,:] = rho[nx-2,:]

  # Top (Zoe-He)
  u[:,ny-1] = 1.0
  v[:,ny-1] = 0.0
  # 3,0,1 + 2.0 * 6,2,5
  rho[:,ny-1] = (1 / (1 + v[:,ny-1])) * (fsB[:, ny-1, 0] + fsB[:, ny-1, 1] + fsB[:, ny-1, 3]) \
                                + 2.0 * (fsB[:, ny-1, 2] + fsB[:, ny-1, 5] + fsB[:, ny-1, 6])

  return u,v,rho


# LBM for Cavity
# Boltzmann equation with BGK model
# Sigle-relaxation-time
# g: Boltzmann-Maxwellian distribution function
# gに低マッハ近似を施すことで2次精度のfeqを得る（高次精度化可能）

nx = 10
ny = 10
Lx = 10.0
Ly = 10.0
ntime = 5
Re = 10
nu = 1 / Re
tau= (6.0 * nu + 1) / 2.0
print("tau is "+str(tau))

# D2Q9 model
e = np.array([[0,0],
              [1,0],
              [0,1],
              [-1,0],
              [0,-1],
              [1,1],
              [-1,1],
              [-1,-1],
              [1,-1]])
w = np.array([4/9, 1/9, 1/9, 1/9, 1/9, 1/36,1/36, 1/36, 1/36])

# Step 1
# Init rho, u, fi, feqi
x = np.linspace(0, nx-1)
y = np.linspace(0, ny-1)
# rho, u, v
rho = np.zeros((nx,ny))
u = np.zeros((nx,ny))
v = np.zeros((nx,ny))
# Distribution
f   = np.zeros((nx,ny,9))
fs  = np.zeros((nx,ny,9))
fc  = np.zeros((nx,ny,9))
fsB = np.zeros((nx,ny,9))
feq = np.zeros((nx,ny,9))

# Initial Condition is weight
for i in range(9):
  f[:,:,i]   = w[i]
rho = calc_rho(f)
u, v = calc_vel(f,e,rho)
u,v,rho = setboundary_vel(nx,ny,u,v,rho,f)
# print("rho  is "+str(rho.mean()))
# print("Uvel is "+str(u.mean()))
# print("Vvel is "+str(v.mean()))

for num in range(ntime):
  if(np.mod(num, 10)==0):
    print(num+1)
    print("rho is "+str(rho.min()))
  # Collision
  feq = calc_feq(nx,ny,w,rho,u,v)
  fc = f - (1/tau) * (f - feq)
  # Streaming
  fs = streaming(nx,ny,fc,e)
  # BounceBack
  f = setBounceBack(fs, rho, u, v)
  # macroscopic quantities
  rho = calc_rho(f)
  u,v = calc_vel(f,e,rho)
  # Boundary for macroscopic quantities
  u,v,rho = setboundary_vel(nx,ny,u,v,rho,f)


print(np.shape(u))
p = 1/3 * rho
plt.subplots(ncols=3, figsize=(10,4))
plt.subplot(131)
# rho
plt.imshow(np.rot90(p), cmap="jet")
plt.colorbar()
# U
plt.subplot(132)
plt.imshow(np.rot90(u), cmap="jet")
plt.colorbar()
# V
plt.subplot(133)
plt.imshow(np.rot90(v), cmap="jet")
plt.colorbar()
plt.show()
