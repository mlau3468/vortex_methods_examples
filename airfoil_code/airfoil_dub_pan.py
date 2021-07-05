import numpy as np
import matplotlib.pyplot as plt
import math
from airfoil_util import *


def proc_panels(pts, debug=False):
    # Processes panels, gives collocation points, orientations, tangents, and normals
    # Points must be numpy array
    # Computes collocation points for panels defined by list of points
    co_pts = np.array([np.mean(pts[i:i+2,:], axis=0) for i in range(pts.shape[0]-1)])
    # computes panel orientation angles a for panels defined by list of points
    # positive going clockwise, starting at -x axis
    dz = np.diff(pts[:,1])
    dx = np.diff(pts[:,0])
    theta = np.arctan2(dz, dx)
    s_theta = -np.sin(theta)
    c_theta = np.cos(theta)
    norms = np.transpose(np.array([s_theta, c_theta]))
    tans = np.transpose(np.array([c_theta, -s_theta]))

    #tangent vector but aligned with free stream
    orients = [math.atan2( (pts[i,1]-pts[i+1,1])  ,  (pts[i+1, 0]-pts[i, 0]))   for i in range(pts.shape[0]-1)]
    tan_fs = [np.array([math.cos(a), -math.sin(a)]) for a in orients]
    tan_fs = np.array(tan_fs)
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

    return co_pts, norms, tans, tan_fs, lens

def dub2D(mu, p1, p2, p):
    # p1 is set as the coordinate for doublet origin
    # p2 then rotated to land on the doublet x axis
    if p2 is None:
        # assume that theta is 0 as this is a wake
        theta = 0
        # x2 -> Inf, z -> 0
        x1 = 0
        x = p[0] - p1[0]
        z = p[1] - p1[1]
        r_2 = (x-x1)**2+z**2
        u = mu/2/math.pi/r_2 * z
        w = -mu/2/math.pi / r_2 * (x)
        uw = np.array([u,w])
    else:
        # coordinate transform
        theta = math.atan2(p2[1]-p1[1] , p2[0]-p1[0])
        theta = -theta
        t1 = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
        # take p1 as coordinate of the origin of the dublet
        t2 = p - p1
        p_new = np.matmul(t1,t2)
        
        x1 = 0
        x2 = pt_dist(p1, p2)
        x = p_new[0]
        z = p_new[1]
        if abs(z) < 0.000001 and (x1 < x and x < x2):  # Influence on itself 
            u = 0
            w = -mu/2/math.pi * (1/(x-x1) - 1/(x-x2))
        elif abs(z) < 0.000001 and (x1 == x or x == x2):
            raise Exception('Singularity')
        else:
            u = mu/2/math.pi * (z/((x-x1)**2+z**2) - z/((x-x2)**2+z**2))
            w = -mu/2/math.pi * ((x-x1)/((x-x1)**2+z**2) - (x-x2)/((x-x2)**2+z**2))
            
        uw = np.array([u,w])
        # coordinate transform back into global frame
        t1 = np.array([[math.cos(theta), math.sin(theta)], [-math.sin(theta), math.cos(theta)]])
        uw = np.matmul(t1,uw)

    return uw


U = 1
chord = 1
alfa = 5
alfa = math.radians(alfa)
u_vec = U * np.array([math.cos(alfa), math.sin(alfa)])
roh = 1.225
#pts = read_dat('airfoil.dat')
pts = read_csv('airfoil.csv')
#pts = repanel(pts, 50, chord_dist = 'cosine', cos_wgt=1.0, show=True)
co_pts, norms, tans,  tan_fs, lens = proc_panels(pts, debug=False)

# Initialize matrices
num_pan = len(pts) - 1
a = np.zeros([num_pan+1, num_pan+1])
b = np.copy(a)
rhs = np.zeros(num_pan+1)

# Set kutta condition
a[-1,0] = 1
a[-1, -2] = -1
a[-1, -1] = 1
rhs[-1] = 0
for i in range(num_pan):
    rhs[i] = - np.matmul(u_vec, norms[i])
    for j in range(num_pan + 1):
        if j == num_pan: # wake
            uw = dub2D(1, pts[num_pan], None, co_pts[i])
        else:
           uw = dub2D(1, pts[j], pts[j+1], co_pts[i])
        a[i,j] = np.matmul ( np.transpose(uw), norms[i])
        b[i,j] = np.matmul(np.transpose(uw), tans[i])
# Remove one panel

rem_idx = int(round(num_pan/2)) -1
rhs = np.delete(rhs, rem_idx)
a = np.delete(a, rem_idx, 0)
a = np.delete(a, rem_idx, 1)
b = np.delete(b, rem_idx, 0)
b = np.delete(b, rem_idx, 1)
norms = np.delete(norms, rem_idx, 0)
tans = np.delete(tans, rem_idx, 0)
tan_fs = np.delete(tan_fs, rem_idx, 0)
pts = np.delete(pts, rem_idx, 0)
co_pts = np.delete(co_pts, rem_idx, 0)
lens = np.delete(lens, rem_idx, 0)

num_pan = num_pan - 1

# solve
mu = np.linalg.solve(a, rhs)
# Induced tangential velocities at each collocation point
q_ti = np.matmul(b, mu)
q_ti = q_ti[:-1]  # dump q_ti of wake
# At each local tangent velocity, add freestream contribution and panel contribution on itself
#11.38
for i in range(num_pan):
    
    if i == 0:  # First panel
        v_pan = 0.5*(mu[i+1]-mu[i])/pt_dist(co_pts[i], co_pts[i+1])
    elif i == num_pan - 1:  # last panel
        v_pan = 0.5*(mu[i]-mu[i-1])/pt_dist(co_pts[i-1], co_pts[i])
    else:
         v_pan = 0.5*(mu[i+1]-mu[i-1])/pt_dist(co_pts[i-1], co_pts[i+1])
    
    q_ti[i] = q_ti[i] + v_pan + np.dot(u_vec, tan_fs[i])
    
    #q_ti[i] = q_ti[i] + np.dot(u_vec, tan_fs[i])
# calculate Cps
cl = 0
cp = 1 - q_ti**2/U**2
ps = cp*0.5*roh*U**2
fs = 0
for i in range(len(ps)):
    f = ps[i]*lens[i]*-1*norms[i,1]
    fs = fs + f
cl = np.sum(fs)/(1/2*roh*chord*U**2)
print('Cl: {:.6f}'.format(cl))

plt.plot(co_pts[:,0], cp)
plt.grid(True)
plt.gca().invert_yaxis()
plt.ylabel('Cp')
plt.xlabel('x')
plt.show()