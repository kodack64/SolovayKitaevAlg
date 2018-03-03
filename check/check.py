import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
import glob

flist = glob.glob("result_*.txt")
ids = [int(fn.replace(".txt","").replace("result_","")) for fn in flist]
ids.sort()

fig = plt.figure(figsize=(8,6))
ax = Axes3D(fig)
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_zlim(-1,1)
for ind in ids:
    data = np.loadtxt("result_{}.txt".format(ind))
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]
    ax.scatter(x,y,z,s=1)
    plt.savefig("result_{}.png".format(ind))
plt.clf()

fig = plt.figure(figsize=(8,6))
ax = Axes3D(fig)
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_zlim(-1,1)
for ind in ids:
    data = np.loadtxt("result_{}.txt".format(ind))
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]
    norm = np.sqrt(x**2+y**2+z**2)
    ax.scatter(x/norm,y/norm,z/norm,s=1)
    plt.savefig("result_{}_norm.png".format(ind))
plt.clf()
