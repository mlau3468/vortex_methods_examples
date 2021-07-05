from shapely.geometry import LineString, Point
import numpy as np
import matplotlib.pyplot as plt
import math

def read_csv(filename):
    return np.genfromtxt(filename, delimiter=',')


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

def read_dat(filename):
    with open(filename, 'r') as fh:
        points = []
        lines = fh.readlines()
        for l in lines:
            l = l.split()
            l = [float(x) for x in l]
            points.append(np.array(l))
    return np.array(points)

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

def pt_dist(pt1, pt2):
    return math.sqrt( (pt2[0] - pt1[0])**2 + (pt2[1] - pt1[1])**2)

def coordRot2D(p, angle, origin):
    t = [[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]]
    t = np.array(t)
    return np.matmul(t, p-origin)