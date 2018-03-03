
import numpy as np
import math
from sympy import *

t = symbols("t",real=True,positive=True)
c,s = cos(t),sin(t)
I = np.array([[1,0],[0,1]])
X = np.array([[0,1],[1,0]])
Z = np.array([[1,0],[0,-1]])
Y = np.array([[0,-1j],[1j,0]])

u = np.dot(np.dot(np.dot((c*I - 1j*s*X),(c*I - 1j*s*Y)),(c*I + 1j*s*X)),(c*I + 1j*s*Y))
#u = np.dot(np.dot((c*I - 1j*s*Y),(c*I + 1j*s*X)),(c*I + 1j*s*Y))
ti = np.trace(u)/2
tx = np.trace(np.dot(1j*X,u))/2
ty = np.trace(np.dot(1j*Y,u))/2
tz = np.trace(np.dot(1j*Z,u))/2

print(simplify(ti))
print(simplify(tx))
print(simplify(ty))
print(simplify(tz))

tnx = tx/math.sqrt(1-simplify(ti)**2)
tny = ty/math.sqrt(1-simplify(ti)**2)
tnz = tz/math.sqrt(1-simplify(ti)**2)
print(simplify(tnx))
print(simplify(tny))
print(simplify(tnz))
