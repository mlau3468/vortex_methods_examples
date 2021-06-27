import math
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
from numpy import cross, dot, mean

def vrtxline(p1, p2, p, gam=1):
    r1 = p-p1
    r2 = p2-p1
    r0 = r1 - r2
    # check for singular conditions
    e = 0.0001
    if norm(r1) < e or norm(r2) < e or norm(cross(r1,r2))**2 < e:
        vel = np.array([0,0,0])
    else:
        K = gam/(4*math.pi*cross(r1,r2)**2)*(dot(r0,r1)/norm(r1)-dot(r0,r2)/norm(r2))
        vel = K*cross(r1,r2)
    return vel

def vrtxring(p1, p2, p3, p4, p, gam=1):
    # clockwise panels
    vel1 = vrtxline(p1, p2, p)
    vel2 = vrtxline(p2, p3, p)
    vel3 = vrtxline(p3, p4, p)
    vel4 = vrtxline(p4, p1, p)
    return vel1 + vel2 + vel3 + vel4

class Panel:
    def __init__(self, pts, vel):
        # defined using right hand rotation, upstream, left point is first
        self.gam = 0
        self.dgdt = 0
        self.pts = pts
        self.vel = vel
        A = pts[2] - pts[0]
        B = pts[3] - pts[1]
        self.norm = cross(A,B)/norm(cross(A,B))
        # place collocation point at 3 quarter chord
        self.cpt = 0.75*(pts[2]+pts[3])/2 + 0.25 * (pts[0]+pts[1])/2 # collocation point
        # place vortex ring at quarter chord, 2D kutta condition satisfid along the chord
        mcv = ((pts[3]-pts[0]) + (pts[2]-pts[1]))/2 # mean chord vector
        self.rpts = [p + 0.25*mcv for p in self.pts]# vortex ring points

nspan = 10
nchord = 5

chord = 1
span = 10

# create panels for rectangular wing
panels = []
for i in range(nchord):
    for j in range(nspan):
        p1 = i*chord/nchord, j*span/nspan
        p2 = i*chord/nchord, (j+1)*span/nspan
        p3 = (i+1)*chord/nchord, (j+1)*span/nspan
        p4 = (i+1)*chord/nchord, j*span/nspan
        p1 = np.array(p1)
        p2 = np.array(p2)
        p3 = np.array(p3)
        p4 = np.array(p4)
        vel = np.array([-1, 0, 0])
        new_pan = Panel([p1, p2, p3, p4], vel)
        panels.append(new_pan)

for p in panels:
    # show panels
    plt.plot([p[0] for p in p.pts], [p[1] for p in p.pts], color='black', linestyle='-')
    plt.plot([p.pts[3][0],p.pts[0][0]], [p.pts[3][1],p.pts[0][1]], color='black', linestyle='-')
for p in panels:
    # show vortex rings
    plt.plot([p[0] for p in p.rpts], [p[1] for p in p.rpts], color='red', linestyle='--')
    plt.plot([p.rpts[3][0],p.rpts[0][0]], [p.rpts[3][1],p.rpts[0][1]], color='red', linestyle='--')
for p in panels:
    # show collocations points
    plt.plot(p.cpt[0], p.cpt[1], color='blue', markersize=5, marker='o')

plt.show()

# -------------------------------------------------------

