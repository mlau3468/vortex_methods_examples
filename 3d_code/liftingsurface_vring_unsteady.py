import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.linalg import norm
from numpy import cross, dot, mean
import vtk
from vtk.util import numpy_support

def vrtxline(p1, p2, p, gam=1):
    r1 = p-p1
    r2 = p-p2
    r0 = p2-p1
    # check for singular conditions
    e = 1e-8
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

def writevtk(panels, wake_rings, fname):
    n = 0
    points = vtk.vtkPoints()
    quads = vtk.vtkCellArray()
    pressures = np.zeros(len(panels))
    for i,p in enumerate(panels):      
        pressures[i] = p.dp
        for pt in p.rpts:
            points.InsertNextPoint(pt)
        quad = vtk.vtkQuad()
        quad.GetPointIds().SetId(0, n+0)
        quad.GetPointIds().SetId(1, n+1)
        quad.GetPointIds().SetId(2, n+2)
        quad.GetPointIds().SetId(3, n+3)
        quads.InsertNextCell(quad)
        n = n + 4

    polydata = vtk.vtkPolyData()
    data = numpy_support.numpy_to_vtk(pressures)
    data.SetName('Pressure')
    polydata.SetPolys(quads)
    polydata.SetPoints(points)
    polydata.GetCellData().SetScalars(data)
    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(polydata)
    writer.SetFileName(fname + '_panels.vtk')
    writer.Write()

    n = 0
    points = vtk.vtkPoints()
    quads = vtk.vtkCellArray()
    for p in wake_rings:      
        for pt in p.pts:
            points.InsertNextPoint(pt)
        quad = vtk.vtkQuad()
        quad.GetPointIds().SetId(0, n+0)
        quad.GetPointIds().SetId(1, n+1)
        quad.GetPointIds().SetId(2, n+2)
        quad.GetPointIds().SetId(3, n+3)
        quads.InsertNextCell(quad)
        n = n + 4
    polydata = vtk.vtkPolyData()
    polydata.SetPolys(quads)
    polydata.SetPoints(points)
    writer = vtk.vtkPolyDataWriter()
    writer.SetInputData(polydata)
    writer.SetFileName(fname + '_wake.vtk')
    writer.Write()


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
        B = pts[1] - pts[3]
        self.area = norm(cross(pts[1]-pts[0],pts[3]-pts[0]))
        self.tanj_vec = pts[1]-pts[0]
        self.tanj_len = norm(self.tanj_vec)
        self.tanj_uvec = self.tanj_vec/self.tanj_len
        self.tani_vec = pts[3]-pts[0]
        self.tani_len = norm(self.tani_vec)
        self.tani_uvec = self.tani_vec/self.tani_len
        self.normal = cross(A,B)/norm(cross(A,B))
        self.wake_vel = np.zeros(3) # wake induced velocity
        self.dp = 0 # pressure differential
        self.df = np.zeros(3) # force
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
        self.ptvels = np.zeros((4,3))
        self.atTE = True
    def step(self):
        if not self.atTE:
            for i in [0,1,2,3]:
                self.pts[i] = self.pts[i] + self.ptvels[i,:]*dt
        else:
            for i in [0,1]:
                self.pts[i] = self.pts[i] + te_scale*self.ptvels[i,:]*dt
            for i in [2,3]:
                self.pts[i] = self.pts[i] + self.ptvels[i,:]*dt
            self.atTE = False



# -------------------------------------------------------
nspan = 13
nchord = 4
n_wake = 100

chord = 1
span = 8

S = span*chord

U = 50
alpha = 5

rho = 1.225

dt = 1/16*chord/U
tsteps = 50
te_scale = 0.3

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

for t in range(tsteps):
    
    # build RHS vector
    RHS = np.zeros(len(panels))
    for i in range(len(panels)):
        RHS[i] = -dot(U_inf,panels[i].normal)
        RHS[i] = RHS[i] -dot(panels[i].wake_vel,panels[i].normal)

    # solve matrix for panel gamma
    sol = np.linalg.solve(A,RHS)
    for i in range(len(panels)):
        panels[i].new_gam(sol[i])
    '''
    # local velocity check
    P0 = 1/2*rho*U**2
    L = 0
    for i in range(len(panels)):
        pan_vel = panels[i].wake_vel + U_inf
        for j in range(len(panels)):
            pan_vel = pan_vel + vrtxring(*panels[j].rpts, panels[i].cpt, panels[j].gam)
        print(pan_vel)
    '''

    # shed wake
    new_wake = []
    if len(wake_rings) < n_wake*len(te_idx):
        for i, idx in enumerate(te_idx):
            p1 = panels[idx].rpts[3]
            p2 = panels[idx].rpts[2]

            vel1 = np.copy(U_inf)
            for p in panels:
                vel1 = vel1 + vrtxring(*p.rpts, p1, gam=p.gam)
            for p in wake_rings:
                vel1 = vel1 + vrtxring(*p.pts, p1, gam=p.gam)

            vel2 = np.copy(U_inf)
            for p in panels:
                vel2 = vel2 + vrtxring(*p.rpts, p2, gam=p.gam)
            for p in wake_rings:
                vel2 = vel2 + vrtxring(*p.pts, p2, gam=p.gam)
            p3 = p2 + te_scale*vel2*dt
            p4 = p1 + te_scale*vel1*dt
            gam = panels[idx].last_gam
            new_wake_pan = WakePanel([p1,p2,p3,p4], gam)
            new_wake.append(new_wake_pan)

    # calculate induced velocities at wake points
    for i in range(len(wake_rings)):
        for j in [0,1,2,3]:
            pt = wake_rings[i].pts[j]
            vel = np.copy(U_inf)
            for k in range(len(panels)):
                vel = vel + vrtxring(*panels[k].rpts, pt, panels[k].gam)
            for k in range(len(wake_rings)):
                vel = vel + vrtxring(*wake_rings[k].pts, pt, wake_rings[k].gam)
            wake_rings[i].ptvels[j,:] = vel

    # move the wake
    for i in range(len(wake_rings)):
        wake_rings[i].step()

    wake_rings = wake_rings + new_wake

    # calculate wake induced velocity at each panel
    for i in range(len(panels)):
        wake_vel = np.zeros(3)
        for j in range(len(wake_rings)):
            wake_vel = wake_vel + vrtxring(*wake_rings[j].pts, panels[i].cpt, wake_rings[j].gam)
        panels[i].wake_vel = wake_vel
    
    # pressure calculation
    for i in range(nchord):
        for j in range(nspan):
            idx = i*nspan + j
            val = 0
            if i > 0:
                gam2 = panels[(i-1)*nspan+j].gam
            else:
                gam2 = 0
            val = val + dot(U_inf+panels[idx].wake_vel, panels[idx].tani_uvec)* (panels[idx].gam-gam2)/panels[idx].tani_len
            
            if j > 0:
                gam2 = panels[i*nspan+(j-1)].gam
            else: 
                gam2 = 0
            val = val + dot(U_inf+panels[idx].wake_vel, panels[idx].tanj_uvec)* (panels[idx].gam-gam2)/panels[idx].tanj_len
            val = val + panels[idx].dgdt
            panels[idx].dp = rho*val
            panels[idx].df = -panels[idx].dp*panels[idx].area*panels[idx].normal
        
    # total forces
    total_force = np.zeros(3)
    for i in range(len(panels)):
        total_force = total_force + panels[i].df
    # lift coefficient
    cl = math.cos(math.radians(alpha))*total_force[2] - math.sin(math.radians(alpha))*total_force[0]
    cl = cl/(1/2*rho*U**2*S)
    print("Timestep: {}, CL={}, F={}".format(t+1, cl, total_force))

    writevtk(panels, wake_rings, './3d_code/viz/out_{}'.format(t))