import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm
from numpy import cross, dot, mean

def vrtxline(p1, p2, p, gam=1):
    r1 = p-p1
    r2 = p-p2
    r0 = r1 - r2
    # check for singular conditions
    e = 0.001
    if norm(r1) < e or norm(r2) < e or norm(cross(r1,r2))**2 < e:
        vel = np.array([0,0,0])
    else:
        K = gam/(4*math.pi*norm(cross(r1,r2))**2)*(dot(r0,r1)/norm(r1)-dot(r0,r2)/norm(r2))
        vel = K*cross(r1,r2)
    return vel

def vrtxring(p1, p2, p3, p4, p, gam=1):
    # clockwise panels
    vel1 = vrtxline(p1, p2, p, gam=gam)
    vel2 = vrtxline(p2, p3, p, gam=gam)
    vel3 = vrtxline(p3, p4, p, gam=gam)
    vel4 = vrtxline(p4, p1, p, gam=gam)
    return vel1 + vel2 + vel3 + vel4

def vrtxringstreamwise(p1, p2, p3, p4, p, gam=1):
    # clockwise panels
    vel2 = vrtxline(p2, p3, p, gam=gam)
    vel4 = vrtxline(p4, p1, p, gam=gam)
    return vel2 + vel4

def view(panels, wake_rings):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for p in panels:
        # show panels
        ax.plot3D([p[0] for p in p.pts], [p[1] for p in p.pts], [p[2] for p in p.pts],color='black', linestyle='-')
        ax.plot3D([p.pts[3][0],p.pts[0][0]], [p.pts[3][1],p.pts[0][1]], [p.pts[3][2],p.pts[0][2]],color='black', linestyle='-')
    for p in panels:
        # show vortex rings
        ax.plot3D([p[0] for p in p.rpts], [p[1] for p in p.rpts], [p[2] for p in p.rpts],color='red', linestyle='--')
        ax.plot3D([p.rpts[3][0],p.rpts[0][0]], [p.rpts[3][1],p.rpts[0][1]],  [p.rpts[3][2],p.rpts[0][2]],color='red', linestyle='--')
    for p in panels:
        # show collocations points
        ax.scatter3D(p.cpt[0], p.cpt[1], p.cpt[2], color='blue')
    for p in wake_rings:
        # show wake rings
        ax.plot3D([p[0] for p in p.pts], [p[1] for p in p.pts], [p[2] for p in p.pts],color='green', linestyle='-')
        ax.plot3D([p.pts[3][0],p.pts[0][0]], [p.pts[3][1],p.pts[0][1]],  [p.pts[3][2],p.pts[0][2]],color='green', linestyle='-')
    plt.show()


class Panel:
    def __init__(self, pts, vel):
        # defined using right hand rotation, upstream, left point is first
        self.gam = 0
        self.last_gam = 0
        self.dgdt = 0
        self.pts = pts
        self.vel = vel
        A = pts[2] - pts[0]
        B = pts[3] - pts[1]
        self.normal = cross(A,B)/norm(cross(A,B))
        # place collocation point at 3 quarter chord
        self.cpt = 0.75*(pts[2]+pts[3])/2 + 0.25 * (pts[0]+pts[1])/2 # collocation point
        # place vortex ring at quarter chord, 2D kutta condition satisfid along the chord
        mcv = ((pts[3]-pts[0]) + (pts[2]-pts[1]))/2 # mean chord vector
        self.rpts = [p + 0.25*mcv for p in self.pts]# vortex ring points
    def new_gam(self, gamma):
        self.last_gam = self.gam
        self.gam = gamma
        self.dgdt = (self.gam-self.last_gam)/dt

class WakePanel():
    def __init__(self, pts, gam):
        self.pts = pts
        self.gam = gam


nspan = 13
nchord = 4

chord = 1
span = 8

U = 1
alpha = 5

dt = 1/16*chord/U
tsteps = 100



U_inf = U*np.array([math.cos(math.radians(alpha)), 0, math.sin(math.radians(alpha))])

# create panels for rectangular wing
panels = []
te_idx = []
last_te_vort = []
for i in range(nchord):
    for j in range(nspan):
        p1 = i*chord/nchord, j*span/nspan, 0
        p2 = i*chord/nchord, (j+1)*span/nspan, 0
        p3 = (i+1)*chord/nchord, (j+1)*span/nspan, 0
        p4 = (i+1)*chord/nchord, j*span/nspan, 0
        p1 = np.array(p1)
        p2 = np.array(p2)
        p3 = np.array(p3)
        p4 = np.array(p4)
        vel = np.array([0, 0, 0])
        new_pan = Panel([p1, p2, p3, p4], vel)
        panels.append(new_pan)
        if i+1 == nchord:
            te_idx.append(i*nspan+j)

# -------------------------------------------------------
# Initizalize matrices
A = np.zeros([len(panels), len(panels)])
# B matrix used for induced drag calculations
B = np.zeros([len(panels), len(panels)])
wake_rings = []

# influence coefficients
for i in range(len(panels)):
    for j in range(len(panels)):
        # influence of jth panel on ith collcation point
        vel =  vrtxring(*panels[j].rpts, panels[i].cpt)
        A[i,j] = dot(vel, panels[i].normal)

        vel2 = vrtxringstreamwise(*panels[j].rpts, panels[i].cpt)
        B[i,j] = dot(vel2, panels[i].normal)

view(panels, wake_rings)

for t in range(tsteps):
    # move wake
    for i in range(len(wake_rings)):
        for j in [0,1,2,3]:
            pt = wake_rings[i].pts[j]
            vel = np.array([0,0,0])
            for k in range(len(panels)):
                vel = vel + vrtxring(*panels[k].rpts, pt, panels[k].gam)
            for k in range(len(wake_rings)):
                vel = vel + vrtxring(*wake_rings[k].pts, pt, wake_rings[k].gam)
            wake_rings[i].pts[j] = wake_rings[i].pts[j] + vel*dt

    # build RHS vector
    RHS = np.zeros(len(panels))
    for i in range(len(panels)):
        RHS[i] = -dot(U_inf,panels[i].normal)
        for j in range(len(wake_rings)):
            RHS[i] = RHS[i] - dot(vrtxring(*wake_rings[j].pts, panels[i].cpt, wake_rings[j].gam), panels[i].normal)
    
    # solve matrix for panel gamma
    sol = np.linalg.solve(A,RHS)
    for i in range(len(panels)):
        panels[i].new_gam(sol[i])

    # shed wake
    for i, idx in enumerate(te_idx):
        p1 = panels[idx].rpts[3]
        p2 = panels[idx].rpts[2]
        p3 = p2 + 0.3*U_inf
        p4 = p1 + 0.3*U_inf
        '''
        if t == 0:
            gam = panels[idx].gam
            last_te_vort.append(gam)
        else:
            gam = last_te_vort[i]
            last_te_vort[i]
        '''
        gam = panels[idx].last_gam
        new_wake_pan = WakePanel([p1,p2,p3,p4], gam)
        wake_rings.append(new_wake_pan)

    view(panels, wake_rings)
