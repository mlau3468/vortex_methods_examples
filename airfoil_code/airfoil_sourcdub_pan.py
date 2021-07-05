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

    r1 = math.sqrt(x**2+z**2)
    r2 = math.sqrt((x-x2)**2+z**2)
    th1 = math.atan2(z,x)
    th2 = math.atan2(z,x-x2)

    return sig/(2*math.pi)*(x*log(r1)* -(x-x2)*log(r2) +z*(th2-th1))

def dubPot(p1, p2, p, mu=1):
    if p2 is None: # means the panel is wake, with p2 at infinity:
        p_new = p-p1
        x1 = 0
        x = p_new[0]
        z = p_new[1]
        return -mu/(2*math.pi) * (-math.atan(z/(x-x1)))
    else:
        # rotate to panel frame
        theta = -math.atan2((p2[1]-p1[1]),(p2[0]-p1[0]))
        p_new = coordRot2D(p, theta, p1)

        x1 = 0
        x2 = pt_dist(p1, p2)

        x = p_new[0]
        z = p_new[1]

        return -mu/(2*math.pi)*(math.atan2(z,(x-x2))-math.atan2(z,(x-x1)))

def proc_panels(pts, debug=False):
    # Processes panels, gives collocation points, orientations, tangents, and normals
    # Points must be numpy array
    # Computes collocation points for panels defined by list of points
    co_pts = np.array([np.mean(pts[i:i+2,:], axis=0) for i in range(pts.shape[0]-1)])
    # computes panel orientation angles a for panels defined by list of points
    # positive going clockwise, starting at -x axis
    dz = np.diff(pts[:,1])
    dx = np.diff(pts[:,0])
    theta = -np.arctan2(dz, dx) #panel orientation. 
    norms = np.transpose(np.array([np.sin(theta), np.cos(theta)]))
    tans = np.transpose(np.array([np.cos(theta), -np.sin(theta)]))

    lens = np.sqrt(dx**2 + dz**2)
    if debug:
        # DEBUG: show normals and tangents
        plt.plot(pts[:,0], pts[:,1], 'k', marker='o')
        plt.plot(co_pts[:,0], co_pts[:,1], marker='x')
        for i in range(len(pts)-1):
            m = 0.1
            plt.plot([co_pts[i][0],co_pts[i][0]+norms[i][0]*m], [co_pts[i][1],co_pts[i][1]+norms[i][1]*m], 'r')
            plt.plot([co_pts[i][0],co_pts[i][0]+tans[i][0]*m], [co_pts[i][1],co_pts[i][1]+tans[i][1]*m], 'b')
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    return co_pts, norms, tans, lens, theta

U = 1
chord = 1
alfa = 5
alfa = math.radians(alfa)
u_vec = U * np.array([math.cos(alfa), math.sin(alfa)])
roh = 1.225
#pts = read_dat('airfoil.dat')
pts = read_csv('airfoil.csv')
#pts = repanel(pts,100, chord_dist = 'cosine', cos_wgt=0.8, show=True)
co_pts, norms, tans, lens, thetas = proc_panels(pts, debug=True)
num_pan = pts.shape[0]-1

U = 1
alpha = 5
alpha = math.radians(alpha)
Uinf = U*np.array([math.cos(alpha), math.sin(alpha)])


A = np.zeros((num_pan+1,num_pan+1))
for i in range(num_pan):
    for j in range(num_pan):
        if i == j:
            c = 1/2
        else:
            c = dubPot(pts[j], pts[j+1], co_pts[i], mu=1)

        A[i,j] = c

A[num_pan,0] = -1
A[num_pan,num_pan] = 1
A[num_pan,num_pan-1] = -1 


B = np.zeros((num_pan,num_pan))
for j in range(num_pan):
    for i in range(num_pan):
        B[i,j] = sourcePot(pts[j], pts[j+1], co_pts[i], sig=1)

source_strenghts = np.zeros(num_pan)
for i in range(num_pan):
    source_strenghts[i] = np.dot(-norms[i], Uinf)

RHS = np.matmul(B, source_strenghts)

RHS = np.append(RHS, 0)

for i in range(num_pan):
    A[i,num_pan] = dubPot(pts[0], None, co_pts[i], mu=1)

sol = np.linalg.solve(A,  RHS)

cps = np.zeros(num_pan-1)
cl_s = np.zeros(num_pan-1)

phi = np.zeros(num_pan)
for i in range(num_pan):
    phi[i] = co_pts[i,0]*math.cos(alfa) + co_pts[i,1]*math.sin(alfa)+sol[i]

for i in range(num_pan-1):
    dl = pt_dist(co_pts[i], co_pts[i+1])
    vel = (phi[i]-phi[i+1])/dl
    cps[i] = 1 - vel**2
    cl_s[i] = -cps[i]*dl*math.cos(thetas[i])/chord

print('Cl: {}'.format(np.sum(cl_s)))
print('Cl: {}'.format(U*sol[-1]))
plt.plot(pts[1:-1,0], cps)
plt.gca().invert_yaxis()
plt.show()