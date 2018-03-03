
import numpy as np

I = [1,0,0,0]
X = [0,1,0,0]
sX = [1,1,0,0]/np.sqrt(2)
sXd = [1,-1,0,0]/np.sqrt(2)
sY = [1,0,1,0]/np.sqrt(2)
sYd = [1,0,-1,0]/np.sqrt(2)
Z = [0,0,0,1]
sZ = [1,0,0,1]/np.sqrt(2)
sZd = [1,0,0,-1]/np.sqrt(2)

def prod(p,q):
    return [
        p[0]*q[0] - p[1]*q[1] - p[2]*q[2] - p[3]*q[3],
        p[0]*q[1] + p[1]*q[0] - p[2]*q[3] + p[3]*q[2],
        p[0]*q[2] + p[1]*q[3] + p[2]*q[0] - p[3]*q[1],
        p[0]*q[3] - p[1]*q[2] + p[2]*q[1] + p[3]*q[0]
    ]
def dag(p):
    return [p[0],-p[1],-p[2],-p[3]]
c1 = [I,X,sX,sXd,sY,sYd]
c2 = [I,Z,sZ,sZd]
clif = []
for p in c2:
    for q in c1:
        clif.append(prod(p,q))

def dist(p,q):
    a = np.sum((np.array(p)-np.array(q))**2)
    b = np.sum((np.array(p)+np.array(q))**2)
    return min(a,b)

for ind in range(len(clif)):
    p = dag(clif[ind])
    ind2 = np.argmin([dist(p,q) for q in clif])
    print(ind,ind2)
