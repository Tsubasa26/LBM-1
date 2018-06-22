
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import shift
# from mpi4py import MPI


def streaming(f):
  # e:0,0 
  fsB[:,:,0] = f[:,:,0]
  # e:1,0 
  fsB[:,:,1] = shift(f[:,:,1], [1,  0], cval=0.0)
  fsB[:,:,2] = shift(f[:,:,2], [0,  1], cval=0.0)
  fsB[:,:,3] = shift(f[:,:,3], [-1, 0], cval=0.0)
  fsB[:,:,4] = shift(f[:,:,4], [0, -1], cval=0.0)
  # e:1,1
  fsB[:,:,5] = shift(f[:,:,5], [1,  1], cval=0.0)
  fsB[:,:,6] = shift(f[:,:,6], [-1, 1], cval=0.0)
  fsB[:,:,7] = shift(f[:,:,7], [-1,-1], cval=0.0)
  fsB[:,:,8] = shift(f[:,:,8], [1, -1], cval=0.0)
  return fsB

def setBounceBack(fsB, rho, u, v):
  # Right 5,1,8
  fsB[nx-1,:,3] += fsB[nx-1,:,1]
  fsB[nx-1,:,7] += fsB[nx-1,:,5]
  fsB[nx-1,:,6] += fsB[nx-1,:,8]
  # Left  6,3,7
  fsB[0,:,1] += fsB[0,:,3]
  fsB[0,:,5] += fsB[0,:,7]
  fsB[0,:,8] += fsB[0,:,6]
  # Bottom 7,4,8
  fsB[:,0,6] += fsB[:,0,8]
  fsB[:,0,2] += fsB[:,0,4]
  fsB[:,0,5] += fsB[:,0,7]

  # Top (Zoe-He)
  # 6, 2, 5
  # 3, 0, 1
  # 7, 4, 8
  fsB[:,ny-1, 4] += fsB[:, ny-1, 2] + 2/3 * rho[:, ny-1] * v[:, ny-1]
  fsB[:,ny-1, 7] += fsB[:, ny-1, 5] + 0.50 * (fsB[:, ny-1, 1] - fsB[:, ny-1, 3]) + (1/6)*rho[:, ny-1]*u[:, ny-1] + 0.50*rho[:, ny-1]*v[:, ny-1]
  fsB[:,ny-1, 8] += fsB[:, ny-1, 6] + 0.50 * (fsB[:, ny-1, 1] - fsB[:, ny-1, 3]) + (1/6)*rho[:, ny-1]*u[:, ny-1] - 0.50*rho[:, ny-1]*v[:, ny-1]

  # reset
  # Right 5,1,8
  fsB[nx-1,:,1] = 0.0
  fsB[nx-1,:,5] = 0.0
  fsB[nx-1,:,8] = 0.0
  # Left  6,3,7
  fsB[0,:,6] = 0.0
  fsB[0,:,3] = 0.0
  fsB[0,:,7] = 0.0
  # Bottom 7,4,8
  fsB[:,0,7] = 0.0
  fsB[:,0,4] = 0.0
  fsB[:,0,8] = 0.0

  # Top 6,2,5
  fsB[:,ny-1,6] = 0.0
  fsB[:,ny-1,2] = 0.0
  fsB[:,ny-1,5] = 0.0


  return fsB


def calc_rho(f):
  rho = f.sum(axis=2)
  return rho

def calc_vel(f, e, rho):
  u = (f * e[:,0]).sum(axis=2) / rho
  v = (f * e[:,1]).sum(axis=2) / rho
  return u, v

def calc_feq(nx,ny, f, w, rho, u, v):
  feq = np.zeros([nx,ny,9])
  for i in range(nx):
    for j in range(ny):
      for l in range(9):
        s1 = w[l] * (3.0*e[l,0]*u[i,j] + 4.5 * (e[l,0]*u[i,j])**2 - 1.5 * u[i,j]**2)
        s2 = w[l] * (3.0*e[l,1]*v[i,j] + 4.5 * (e[l,1]*v[i,j])**2 - 1.5 * v[i,j]**2)
        feq[i,j,l] = w[l] * rho[i,j] + rho[i,j] * np.sqrt(s1**2 + s2**2)
  return feq

def setboundary_vel(nx,ny,u,v):
  # top
  u[:,ny-1] = 1.0
  v[:,ny-1] = 0.0
  # bottom
  u[:,0] = 0.0
  v[:,0] = 0.0
  # Left
  u[0,:] = 0.0
  v[0,:] = 0.0
  # Right
  u[nx-1,:] = 0.0
  v[nx-1,:] = 0.0

  return u,v




nx = 10
ny = 10
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

# rho, u, v
x = np.linspace(0, nx-1)
y = np.linspace(0, ny-1)
rho = np.zeros((nx,ny))
u = np.zeros((nx,ny))
v = np.zeros((nx,ny))
# Distribution
f   = np.zeros((nx,ny,9))
fs  = np.zeros((nx,ny,9))
fsB = np.zeros((nx,ny,9))
feq = np.zeros((nx,ny,9))

# Initial Condition is weight
for i in range(9):
  f[:,:,i]=w[i]

nu = 1 / 18
tau= (6.0 * nu + 1) / 2.0
for num in range(1000):
  if(np.mod(num,10)==0):
    print(num)
    print(f)
  # step2
  fsB = streaming(f)
  # BounceBack?
  fs = setBounceBack(fsB,rho, u, v)
  # step3
  rho = calc_rho(fs)
  u,v = calc_vel(fs,e,rho)
  # Boundary for velocity ?
  u,v = setboundary_vel(nx,ny,u,v)
  # step4
  feq = calc_feq(nx,ny,fs,w,rho,u,v)
  # step5
  f = fs - (1/tau) * (fs - feq)


print(np.shape(u))
plt.subplots(ncols=2, figsize=(10,4))
plt.imshow(np.rot90(rho), cmap="jet")
# U
plt.imshow(np.rot90(u), cmap="jet")
# V
plt.imshow(np.rot90(v), cmap="jet")
plt.colorbar()
plt.show()
