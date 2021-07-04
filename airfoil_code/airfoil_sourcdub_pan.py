import numpy as np
import matplotlib.pyplot as plt
import math
from math import log
from shapely.geometry import LineString, Point


class Plotter:
    def __init__(self):
        pass

    def plot(self, item, marker=False):
        if type(item) is list or type(item) is tuple:
            for it in item:
                self.plot(it)
        else:
            if item.geom_type == 'MultiLineString':
                for string in item:
                    coord = list(string.coords)
                    x = [coord[i][0] for i in range(len(coord))]
                    y = [coord[i][1] for i in range(len(coord))]
                    if marker:
                        plt.plot(x, y, marker='o')
                    else:
                        plt.plot(x, y)
            elif item.geom_type == 'Polygon':
                plt.plot(*item.exterior.xy)
            elif item.geom_type == 'LineString':
                coord = list(item.coords)
                x = [coord[i][0] for i in range(len(coord))]
                y = [coord[i][1] for i in range(len(coord))]
                if marker:
                    plt.plot(x, y, marker='o')
                else:
                   plt.plot(x, y)
            elif item.geom_type == 'Point':
                plt.plot(item.x, item.y, 'ro', markersize=1)
            elif item.geom_type == 'MultiPolygon':
                for poly in item:
                    plt.plot(*poly.exterior.xy)
            elif item.geom_type == 'MultiPoint':
                for i in item:
                    plt.plot(i.x, i.y, 'ro', markersize=1)

    def show_plot(self, title=None):
        plt.gca().set_aspect('equal', adjustable='box')
        if title is not None:
            plt.title(title)
        plt.show()


def repanel(points, n_chord, chord_dist='cosine', cos_wgt=1, show=False):

    # check if we have sharp trailing edge
    sharp_Te = points[0][0] == points[-1][0] and  points[0][1] == points[-1][1]
    # find index of trailing edge:
    anc_idx = np.argmax(points[:,0])
    # check that anc_idx is either the first or last point
    if anc_idx != 0 and anc_idx != len(points) - 1:
        raise Exception('Error: Provdie a dat file with start or end point at the trailing edge')

    if sharp_Te: # If trailing edge ends at a singular point
        # Find which way to wrap
        if points[1,1] > points[-2,1]:
            points = points[::-1]
        else:
            pass
    else: # If the trailing edge is open
        if points[0,1] > points[-1,1]:
            points = points[::-1]
        else:
            pass

    # Note: Point Data ALWAYS goes downwards (loops in positive cartesian direction)
    # Get LE Point
    LE_x = min(points[:,0])
    LE_ys = [points[i,1] for i in range(len(points)) if points[i,0] == LE_x]
    LE_y = np.mean(LE_ys)
    LE_coords = [LE_x, LE_y]

    # Resample
    line = LineString(points)
    LE_point = Point(LE_x, LE_y)
    d1 = line.project(LE_point) # Bottom sec full distance
    d2 = line.length - d1 # Top sec full distance

    if chord_dist == 'uniform':
        dists_bot = list(np.linspace(0, d1, num=n_chord + 1))
        dists_top = list(np.linspace(d1, line.length, num=n_chord + 1))
        dists = dists_bot + dists_top[1:-1] + [dists_bot[0]]  # Enforce end where we start
    elif chord_dist == 'cosineLE':
        x = np.linspace(0, math.pi/2*cos_wgt, num=n_chord + 1)
        s = np.sin(x)
        s = (s - min(s)) 
        s = s / max(s)
        dists_bot = list(np.multiply(s, d1))

        x = np.linspace(-math.pi/2*cos_wgt, 0, num=n_chord + 1)
        s = np.sin(x)
        s = (s - min(s)) 
        s = s / max(s)
        dists_top = np.multiply(s, d2)
        dists_top = list(np.add(dists_top, d1))
        dists = dists_bot + dists_top[1:-1] + [dists_bot[0]]
    elif chord_dist == 'cosineTE':
        x = np.linspace(-math.pi/2*cos_wgt, 0, num=n_chord + 1)
        s = np.sin(x)
        s = (s - min(s)) 
        s = s / max(s)
        dists_bot = list(np.multiply(s, d1))

        x = np.linspace(0, math.pi/2*cos_wgt, num=n_chord + 1)
        s = np.sin(x)
        s = (s - min(s)) 
        s = s / max(s)
        dists_top = np.multiply(s, d2)
        dists_top = list(np.add(dists_top, d1))
        dists = dists_bot + dists_top[1:-1] + [dists_bot[0]]
    elif chord_dist == 'cosine':
        x = np.linspace(-math.pi/2*cos_wgt, math.pi/2*cos_wgt, num=n_chord + 1)
        s = np.sin(x)
        s = (s - min(s)) 
        s = s / max(s)
        dists_bot = list(np.multiply(s, d1))
        dists_top = np.multiply(s, d2)
        dists_top = list(np.add(dists_top, d1))
        dists = dists_bot + dists_top[1:-1] + [dists_bot[0]]

    points_rsmpl = []
    for j in range(len(dists)):
        new_pt = line.interpolate(dists[j])
        new_pt = [new_pt.x, new_pt.y]
        points_rsmpl.append(new_pt)

    line_rsmpl = LineString(points_rsmpl)
    if show:
        p = Plotter()
        p.plot(line)
        p.plot(line_rsmpl, marker=True)
        p.plot(LE_point)
        p.show_plot()

    return np.array(points_rsmpl)

def ptDist(p1, p2):
    return math.sqrt((p2[1]-p1[1])**2 + (p2[0]-p1[0])**2)

def coordRot2D(p, angle, origin):
    t = [[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]]
    t = np.array(t)
    return np.matmul(t, p-origin)


def sourcePot(p1, p2, p, sig=1):
    # rotate to panel frame
    theta = -math.atan2((p2[1]-p1[1]),(p2[0]-p1[0]))
    p_new = coordRot2D(p, theta, p1)

    x1 = 0
    x2 = ptDist(p1, p2)

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
        x2 = ptDist(p1, p2)

        x = p_new[0]
        z = p_new[1]

        return -mu/(2*math.pi)*(math.atan2(z,(x-x2))-math.atan2(z,(x-x1)))

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
pts = repanel(pts,100, chord_dist = 'cosine', cos_wgt=0.8, show=True)
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
    dl = ptDist(co_pts[i], co_pts[i+1])
    vel = (phi[i]-phi[i+1])/dl
    cps[i] = 1 - vel**2
    cl_s[i] = -cps[i]*dl*math.cos(thetas[i])/chord

print('Cl: {}'.format(np.sum(cl_s)))
print('Cl: {}'.format(U*sol[-1]))
plt.plot(pts[1:-1,0], cps)
plt.gca().invert_yaxis()
plt.show()