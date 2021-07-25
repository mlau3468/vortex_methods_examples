import numpy as np
import matplotlib.pyplot as plt
import math
from math import log
from airfoil_util import *
from scipy.optimize import minimize

def sourc2D(p1, p2, p, sig=1):
    theta = -math.atan2(p2[1]-p1[1] , p2[0]-p1[0])
    t1 = np.array([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
    # take p1 as coordinate of the origin of the dublet
    t2 = p - p1
    p_new = np.matmul(t1,t2)
    
    x1 = 0
    x2 = pt_dist(p1, p2)
    x = p_new[0]
    z = p_new[1]

    top = (x-x1)**2+z**2
    bot = (x-x2)**2+z**2
    u = sig/(4*math.pi)*log(top/bot)
    th2 = np.arctan2(z, x-x2)
    th1 = np.arctan2(z,x-x1)
    w = sig/(2*math.pi)*(th2-th1)
    uw = np.array([u,w])
    # coordinate transform back into global frame
    t1 = np.array([[math.cos(theta), math.sin(theta)], [-math.sin(theta), math.cos(theta)]])
    uw = np.matmul(t1,uw)

    return uw


def dub2D(p1, p2, p, mu=1):
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
        theta = -math.atan2(p2[1]-p1[1] , p2[0]-p1[0])
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

def f1(H):
    if H <= 4:
        return 1.515 + 0.076*(H-4)**2/H
    else:
        return 1.515 + 0.04*(H-4)**2/H

def f2(H):
    if H <= 7.4:
        return -0.067 + 0.01977*(7.4-H)**2/(H-1)
    else:
        return -0.067 + 0.022*(1-1.4/(H-6))**2

def f3(H):
    if H <= 4:
        return 0.207 + 0.00205*(4-H)**5.5
    else:
        return 0.207 - 0.003*(H-4)**2

def f1Wake(H):
    if H <= 3.5:
        return 1.50 + 0.025*(3.5-H)**3+0.001*(3.5-H)**5
    else:
        return 1.50 + 0.015*(H-3.5)**2/H

def f2Wake(H):
    return 0

def f3Wake(H):
    return 1.52*(H-1)**2/(3+H**3)

U = 1
chord = 1
alfa = 5
alfa = math.radians(alfa)
Uinf = U * np.array([math.cos(alfa), math.sin(alfa)])
roh = 1.225
pts = read_csv('airfoil.csv')
#pts = read_csv('4521.csv')
#pts = repanel(pts,100, chord_dist = 'cosine', cos_wgt=0.8, show=True)
co_pts, norms, tans, lens, thetas = proc_panels(pts, debug=False)
num_pan = pts.shape[0]-1

Cpot = np.zeros((num_pan+1,num_pan+1)) # doublet influence coefficients
for i in range(num_pan):
    for j in range(num_pan+1):
        if j < num_pan:
            if i == j:
                c = 1/2
            else:
                c = dubPot(pts[j], pts[j+1], co_pts[i], mu=1)
        else:
            c = dubPot(pts[0], None, co_pts[i], mu=1)
        Cpot[i,j] = c

# kutta condition
Cpot[num_pan,0] = 1
Cpot[num_pan, num_pan] = 1
Cpot[num_pan, num_pan-1] = -1

Bpot = np.zeros((num_pan+1,num_pan+1))
for j in range(num_pan):
    for i in range(num_pan):
        Bpot[i,j] = sourcePot(pts[j], pts[j+1], co_pts[i], sig=1)

source_strenghts = np.zeros(num_pan+1)
for i in range(num_pan):
    source_strenghts[i] = np.dot(-norms[i], Uinf)

RHS = np.matmul(Bpot, source_strenghts)

invSol = np.linalg.solve(Cpot,  RHS) # inviscid solution

cps = np.zeros(num_pan-1)

phi = np.zeros(num_pan)
for i in range(num_pan):
    phi[i] = co_pts[i,0]*U*math.cos(alfa) + co_pts[i,1]*U*math.sin(alfa)+invSol[i]

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

# --------------------------------------------------------

nu = 1.46e-5

# matrices for tangent velocity influence
Ctan = np.zeros((num_pan+1, num_pan+1)) # doublet panel tangent influence
for i in range(num_pan):
    for j in range(num_pan+1):
        if j < num_pan:
            vel = dub2D(pts[j], pts[j+1], co_pts[i],mu=1)
            Ctan[i,j] = np.dot(vel, tans[i])
        else: # wake
            vel = dub2D(pts[0], [1e12,0], co_pts[i], mu=1)
            Ctan[i,j] = np.dot(vel, tans[i])

Atan = np.zeros((num_pan+1, num_pan+1)) # source panel tangent influence
for i in range(num_pan): 
    for j in range(num_pan + 1):
        if j < num_pan:
            vel = sourc2D(pts[j], pts[j+1], co_pts[i],sig=1)
            Atan[i,j] = np.dot(vel, tans[i])
        else:
            vel = sourc2D(pts[0], [1e12,0], co_pts[i], sig=1)
            Atan[i,j] = np.dot(vel, tans[i])

Apot2 = np.zeros((num_pan+1, num_pan+1)) # Modified source potential influence
for i in range(num_pan):
    for j in range(num_pan+1):
        if j > 0:
            Apot2[i,j] = (Bpot[i,j]-Bpot[i,j-1])/pt_dist(pts[j], pts[j-1])
        else:
            Apot2[i,j] = (Bpot[i,j])/pt_dist(pts[0], co_pts[j])

mu = invSol # doublet strength
th = 0.01*np.ones(num_pan+1) # momentum thickness
m = 0.01*np.ones(num_pan+1) # mass defect

x0 = np.concatenate([mu, th, m])

lens = np.append(lens, 1e12)
pts = np.vstack([pts, [1e12,0]])
tans = np.vstack([tans, [1,0]])

def calcResidual(in_vec):
    mu = in_vec[:num_pan+1]
    th = in_vec[num_pan+1:2*(num_pan+1)]
    m = in_vec[2*(num_pan+1):3*(num_pan+1)]
    Uelast = None
    Hlast = None
    residual = np.zeros((num_pan+1,3))
    for i in range(num_pan+1):
            Uei = np.dot(Uinf,tans[i]) + np.matmul(Atan[i,:],np.divide(m, lens)) + np.matmul(Ctan[i,:], mu)
            deli = m[i]/Uei
            Hi = deli/th[i]

            if i == num_pan: #wake panel
                Hi_star = f1Wake(Hi)
                cf2 = nu/(Uei*th[i])*f2Wake(Hi)
                cf2H = nu/(Uei*th[i])*f3Wake(Hi)
            else:
                Hi_star = f1(Hi)
                cf2 = nu/(Uei*th[i])*f2(Hi)
                cf2H = nu/(Uei*th[i])*f3(Hi)

            if i > 0:
                del_th = th[i] - th[i-1]
                th_avg = 1/2*(th[i] + th[i-1])
                del_Ue = Uei-Uelast
                del_x = pt_dist(pts[i+1], pts[i])
                # del_x = pt_dist(co_pts[i+1], co_pts[i])
                del_h = Hi-Hlast
            else:
                del_th = th[i]
                th_avg = th[i]/2
                del_Ue = Uei
                del_x = pt_dist(pts[i+1], pts[i])
                # del_X = pt_dist(pts[0], co_pts[i])
                del_h = Hi
                
            Uelast = Uei
            Hlast = Hi
            # residuals
            R1 = del_th/th_avg + (Hi+2)*del_Ue/Uei - cf2*del_x/th_avg
            R2 = del_h/Hi_star + (1-Hi)*del_Ue/Uei + (cf2-cf2H)*del_x/th_avg
            R3 = np.matmul(Apot2[i,:],m) + np.matmul(Cpot[i,:], mu) - RHS[i]
            residual[i,:] = [R1, R2, R3]
        
    res_sum = np.sum(np.sum(np.abs(residual.flatten())))
    print(res_sum)
    return res_sum


options_dict = {'maxiter': 2000, 'disp': True}
result = minimize(calcResidual, x0, tol=1e-5, method='Nelder-Mead', options=options_dict)

