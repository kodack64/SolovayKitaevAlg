
import numpy as np
import math
from sympy import *

#a,b,c,d,q,r,s = symbols("a b c d q r s",real=True)
t1,t2,t3,t4,t5 = symbols("t1 t2 t3 t4 t5",real=True)
a = cos(t1)
b = sin(t1)*cos(t2)*cos(t3)
c = sin(t1)*cos(t2)*sin(t3)
d = sin(t1)*sin(t2)
q = sin(t1)*cos(t4)*cos(t5)
r = sin(t1)*cos(t4)*sin(t5)
s = sin(t1)*sin(t4)

I = np.array([[1,0],[0,1]])
X = np.array([[0,1],[1,0]])
Z = np.array([[1,0],[0,-1]])
Y = np.array([[0,-1j],[1j,0]])

u = a*I - 1j*b*X - 1j*c*Y - 1j*d*Z
v = a*I - 1j*q*X - 1j*r*Y - 1j*s*Z
un = sqrt(1-a**2)
nx = b/un
ny = c/un
nz = d/un
mx = q/un
my = r/un
mz = s/un
wx = (nx+mx)/2
wy = (ny+my)/2
wz = (nz+mz)/2
wn = sqrt(wx**2+wy**2+wz**2)
w = -1j*wx/wn*X -1j*wy/wn*Y -1j*wz/wn*Z
uu = np.dot(np.dot(w,u),w)

ti = np.trace(uu)/2
tx = np.trace(np.dot(1j*X,uu))/2
ty = np.trace(np.dot(1j*Y,uu))/2
tz = np.trace(np.dot(1j*Z,uu))/2
#print(simplify(ti))
#print(simplify(tx))
#print(simplify(ty))
#print(simplify(tz))
vi = np.trace(v)/2
vx = np.trace(np.dot(1j*X,v))/2
vy = np.trace(np.dot(1j*Y,v))/2
vz = np.trace(np.dot(1j*Z,v))/2
#print(simplify(vi))
#print(simplify(vx))
#print(simplify(vy))
#print(simplify(vz))
ui = np.trace(u)/2
ux = np.trace(np.dot(1j*X,u))/2
uy = np.trace(np.dot(1j*Y,u))/2
uz = np.trace(np.dot(1j*Z,u))/2


udd = [0.987489,0.121766,-0.0813614,0.058474]
vwvw = [0.987489,0.0426913,-0.0426913,-0.145673]
vt1 = -acos(udd[0])
vt2 = asin(udd[3]/sin(vt1))
vt3 = acos(udd[1]/sin(vt1)/cos(vt2))
vt4 = asin(vwvw[3]/sin(vt1))
vt5 = acos(vwvw[1]/sin(vt1)/cos(vt4))
sets = [(t1,vt1),(t2,vt2),(t3,vt3),(t4,vt4),(t5,vt5)]

print(simplify(ti.subs(sets)))
print(simplify(tx.subs(sets)))
print(simplify(ty.subs(sets)))
print(simplify(tz.subs(sets)))
print()
print(simplify(ui.subs(sets)))
print(simplify(ux.subs(sets)))
print(simplify(uy.subs(sets)))
print(simplify(uz.subs(sets)))
print()
print(simplify(vi.subs(sets)))
print(simplify(vx.subs(sets)))
print(simplify(vy.subs(sets)))
print(simplify(vz.subs(sets)))
print()
print(nx.subs(sets))
print(ny.subs(sets))
print(nz.subs(sets))
print()
print(mx.subs(sets))
print(my.subs(sets))
print(mz.subs(sets))
print()
print((wx/wn).subs(sets))
print((wy/wn).subs(sets))
print((wz/wn).subs(sets))
