import numpy as np
import matplotlib.pyplot as plt
import math

def ptDist(p1, p2):
    return math.sqrt((p2[1]-p1[1])**2 + (p2[0]-p1[0])**2)

def coordRot2D(p, angle, origin):
    t = [[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]]
    t = np.array(t)
    return np.matmul(t, p-origin)


def sourcePot(p1, p2, p, sig=1):
    # rotate to panel frame
    theta = -np.arctan((p2[1]-p1[1])/(p2[0]-p1[0]))
    p_new = coordRot2D(p, theta, p1)

    x1 = 0
    x2 = ptDist(p1, p2)

    x = p_new[0]
    z = p_new[1]
    if z == 0:
        a = (x-x1)*np.log((x-x1)**2)
        b = (x-x2)*np.log((x-x2)**2)
        return sig/(4*math.pi)*(a-b)
    else:

        a = (x-x1)*np.log((x-x1)**2+z**2)
        b = (x-x2)*np.log((x-x2)**2+z**2)
        c = np.arctan(z/(x-x2))-np.arctan(z/(x-x1))

        return sig/(4*math.pi)*(a-b+2*z*c)

def dubPot(p1, p2, p, mu=1):
    if p2 is None: # means the panel is wake, with p2 at infinity:
        p_new = p-p1
        x1 = 0
        x = p_new[0]
        z = p_new[1]
        return -mu/(2*math.pi) * (-np.arctan(z/(x-x1)))
    else:
        # rotate to panel frame
        theta = -np.arctan((p2[1]-p1[1])/(p2[0]-p1[0]))
        p_new = coordRot2D(p, theta, p1)

        x1 = 0
        x2 = ptDist(p1, p2)

        x = p_new[0]
        z = p_new[1]

        return -mu/(2*math.pi)*(np.arctan(z/(x-x2))-np.arctan(z/(x-x1)))

def read_dat(filename):
    with open(filename, 'r') as fh:
        points = []
        lines = fh.readlines()
        for l in lines:
            l = l.split()
            l = [float(x) for x in l]
            points.append(np.array(l))
    return np.array(points)

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
pts = read_dat('airfoil.dat')
co_pts, norms, tans, lens, thetas = proc_panels(pts, debug=True)
num_pan = pts.shape[0]-1

U = 1
alpha = 5
alpha = math.radians(alpha)
Uinf = U*np.array([math.cos(alpha), math.sin(alpha)])
c = 0

A = np.zeros((num_pan,num_pan))
for i in range(num_pan):
    for j in range(num_pan):
        if i == j:
            c = 1/2
        else:
            c = dubPot(pts[j], pts[j+1], co_pts[i], mu=1)

        if j == 0:
            ciw = dubPot(pts[0], None, co_pts[i], mu=1)
            A[i,j] = c - ciw
        elif j == num_pan-1:
            ciw = dubPot(pts[0], None, co_pts[i], mu=1)
            A[i,j] = c + ciw
        else:
            A[i,j] = c

B = np.zeros((num_pan,num_pan))
for i in range(num_pan):
    for j in range(num_pan):
        B[i,j] = sourcePot(pts[j], pts[j+1], co_pts[i], sig=1)

source_strenghts = np.zeros(num_pan)
for i in range(num_pan):
    source_strenghts[i] = np.dot(-norms[i], Uinf)

RHS = -np.matmul(B, source_strenghts)

sol = np.linalg.solve(A,  RHS)


cps = np.zeros(num_pan-1)
cl_s = np.zeros(num_pan-1)

qtj = np.zeros(num_pan-1)
for i in range(num_pan-1):
    dl = ptDist(co_pts[i], co_pts[i+1])
    qtj[i] = (sol[i]-sol[i+1])/dl + np.dot(Uinf, tans[i])
    cps[i] = 1 - qtj[i]**2/U**2
    cl_s[i] = -cps[i]*dl*math.cos(thetas[i])/c

print(cps)
print('Cl: {}'.format(np.sum(cl_s)))
plt.plot(co_pts[1:,0], cps)
plt.show()