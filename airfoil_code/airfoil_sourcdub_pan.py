import numpy as np
import matplotlib.pyplot as plt
import math
from math import log
from airfoil_util import *


def sourcePot(p1, p2, p, sig=1):
    # rotate to panel frame
    theta = -math.atan2((p2[1]-p1[1]),(p2[0]-p1[0]))
    p_new = coordRot2D(p, theta, p1)

    x1 = 0
    x2 = pt_dist(p1, p2)

    x = p_new[0]
    z = p_new[1]

    r12 = x**2+z**2
    r22 = (x-x2)**2+z**2
    th1 = math.atan2(z,x)
    th2 = math.atan2(z,x-x2)
    f = x*log(r12) - (x-x2)*log(r22) + 2*z*(th2-th1)
    return sig/(math.pi*4) * f

def dubPot(p1, p2, p, mu=1):
    if p2 is None: # means the panel is wake, with p2 at infinity:
        theta1 = math.atan2(p[1]-p1[1],p[0]-p1[0])
        theta2 = math.atan2(p[1]-p1[1],p[0]-10000000*p1[0])
        return mu/(2*math.pi) * (theta1-theta2)
    else:
        # rotate to panel frame
        theta = -math.atan2((p2[1]-p1[1]),(p2[0]-p1[0]))
        p_new = coordRot2D(p, theta, p1)

        x1 = 0
        x2 = pt_dist(p1, p2)

        x = p_new[0]
        z = p_new[1]

        eta1 = math.atan2(z,x)
        eta2 = math.atan2(z,x-x2)

        return -mu/(2*math.pi)*(eta2-eta1)

U = 1
chord = 1
alfa = 2
alfa = math.radians(alfa)
Uinf = U * np.array([math.cos(alfa), math.sin(alfa)])
roh = 1.225
#pts = read_csv('airfoil.csv')
pts = read_csv('4521.csv')
#pts = repanel(pts,100, chord_dist = 'cosine', cos_wgt=0.8, show=True)
co_pts, norms, tans, lens, thetas = proc_panels(pts, debug=True)
num_pan = pts.shape[0]-1


A = np.zeros((num_pan,num_pan))
for i in range(num_pan):
    for j in range(num_pan):
        if i == j:
            c = 1/2
        else:
            c = dubPot(pts[j], pts[j+1], co_pts[i], mu=1)
        if j == 0:
            ciw = dubPot(pts[0], None, co_pts[i], mu=1)
            c = c - ciw
        elif j == num_pan-1:
            ciw = dubPot(pts[0], None, co_pts[i], mu=1)
            c = c + ciw

        A[i,j] = c

B = np.zeros((num_pan,num_pan))
for j in range(num_pan):
    for i in range(num_pan):
        B[i,j] = sourcePot(pts[j], pts[j+1], co_pts[i], sig=1)

source_strenghts = np.zeros(num_pan)
for i in range(num_pan):
    source_strenghts[i] = np.dot(-norms[i], Uinf)

RHS = np.matmul(B, source_strenghts)

sol = np.linalg.solve(A,  RHS)

cps = np.zeros(num_pan-1)

phi = np.zeros(num_pan)
for i in range(num_pan):
    phi[i] = co_pts[i,0]*U*math.cos(alfa) + co_pts[i,1]*U*math.sin(alfa)+sol[i]

cl = 0
for i in range(num_pan-1):
    dl = pt_dist(co_pts[i], co_pts[i+1])
    qti = (phi[i]-phi[i+1])/dl
    cps[i] = 1 - qti**2/U**2
    ps = cps[i]*1/2*roh*U**2
    f = -ps*dl*norms[i,1]
    cl = cl + f/(1/2*roh*U**2*chord)

print('CL= {}'.format(cl))
plt.plot(pts[1:-1,0], cps)
plt.gca().invert_yaxis()
plt.show()