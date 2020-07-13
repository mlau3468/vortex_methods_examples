import numpy as np
import math
import matplotlib.pyplot as plt
c = 1
b = 5
sweep = math.radians(45)
n_span = 4
x_end = 10
u_ref = 10
a = math.radians(5)
fi = 0
roh_ref = 1.225

# bound_vortex points
bnd_vrt = []
ctrl_pts = []
db = b/n_span/2
for i in range(n_span):
    p1 = 0.25 * c + db * math.tan(sweep) * i, db*i, 0
    p2 = 0.25 * c + db * math.tan(sweep) * (i+1), db*(i+1), 0
    p3 = x_end, p1[1], 0
    p4 = x_end, p2[1], 0

    cpt = (p1[0] + 0.5*c + p2[0] + 0.5*c)/2, (p1[1] + p2[1])/2, 0

    bnd_vrt.append([p1, p2])
    ctrl_pts.append(cpt)

    plt.plot((p3[0], p1[0], p2[0], p4[0]), (p3[1], p1[1], p2[1], p4[1]), 'k')
    plt.plot(cpt[0], cpt[1], 'k.')

for i in range(n_span):
    p1 = 0.25 * c + db * math.tan(sweep) * i, db*-i, 0
    p2 = 0.25 * c + db * math.tan(sweep) * (i+1), db*(-i-1), 0
    p3 = x_end, p1[1], 0
    p4 = x_end, p2[1], 0

    cpt = (p1[0] + 0.5*c + p2[0] + 0.5*c)/2, (p1[1] + p2[1])/2, 0

    bnd_vrt.append([p1, p2])
    ctrl_pts.append(cpt)

    plt.plot((p3[0], p1[0], p2[0], p4[0]), (p3[1], p1[1], p2[1], p4[1]), 'k')
    plt.plot(cpt[0], cpt[1], 'k.')

gam_inflce_mat = np.zeros((len(ctrl_pts), len(bnd_vrt)))
norm_v_ref = np.zeros(len(ctrl_pts))

for m in range(len(ctrl_pts)):
    xm = ctrl_pts[m][0]
    ym = ctrl_pts[m][1]
    zm = ctrl_pts[m][2]
    for n in range(len(bnd_vrt)):
        # Boundary vortex must go from -y to +y
        idx_1 = np.argmin((bnd_vrt[n][0][1], bnd_vrt[n][1][1]))
        if idx_1 == 0:
            idx_2 = 1
        elif idx_1 == 1:
            idx_2 = 0
        x1n = bnd_vrt[n][idx_1][0]
        x2n = bnd_vrt[n][idx_2][0]
        y1n = bnd_vrt[n][idx_1][1]
        y2n = bnd_vrt[n][idx_2][1]
        z1n = bnd_vrt[n][idx_1][2]
        z2n = bnd_vrt[n][idx_2][2]

        if y1n > y2n:
            Exception('Bounded vortex direction erro')

        # Segment AB
        fac1_num = (ym-y1n)*(zm-z2n)-(ym-y2n)*(zm-z1n), -1*((xm-x1n)*(zm-z2n)-(xm-x2n)*(zm-z1n)), (xm-x1n)*(ym-y2n) - (xm-x2n)*(ym-y1n)
        fac1_den = fac1_num[0]**2 + fac1_num[1]**2 + fac1_num[2]**2
        fac1 = np.divide(fac1_num, fac1_den)
        fac2_num1 = (x2n-x1n)*(xm-x1n) + (y2n-y1n)*(ym-y1n) + (z2n-z1n)*(zm-z1n)
        fac2_den1 = math.sqrt((xm-x1n)**2 + (ym-y1n)**2 + (zm-z1n)**2)  
        fac2_num2 = (x2n-x1n)*(xm-x2n) + (y2n-y1n)*(ym-y2n) + (z2n-z1n)*(zm-z2n)
        fac2_den2 = math.sqrt((xm-x2n)**2 + (ym-y2n)**2 + (zm-z2n)**2)  
        fac2 = fac2_num1 / fac2_den1 - fac2_num2 / fac2_den2
        c_ab = np.multiply(fac1, fac2 / math.pi / 4)                  

        # Segment A_inf
        fac1_num = 0, (zm-z1n), (y1n-ym)
        fac1_den = fac1_num[1]**2 + fac1_num[2]**2
        fac1 = np.divide(fac1_num, fac1_den)
        fac2_num = xm - x1n
        fac2_den = math.sqrt((xm-x1n)**2 + (ym-y1n)**2 + (zm-z1n)**2)
        fac2 = 1 + fac2_num  / fac2_den
        ca_inf = np.multiply(fac1, fac2/math.pi/4)

        # Segment B_inf
        fac1_num = 0, (zm-z2n), (y2n-ym)
        fac1_den = fac1_num[1]**2 + fac1_num[2]**2
        fac1 = np.divide(fac1_num, fac1_den)
        fac2_num = xm - x2n
        fac2_den = math.sqrt((xm-x2n)**2 + (ym-y2n)**2 + (zm-z2n)**2)
        fac2 = 1 + fac2_num  / fac2_den
        cb_inf = np.multiply(fac1, fac2/math.pi/-4)

        c_tot = np.add(np.add(c_ab, ca_inf), cb_inf)
        gam_inflce_mat[m, n] = c_tot[2]

    dzdx = 0
    norm_v_ref[m] = dzdx - u_ref * a

gam_sol = np.linalg.solve(gam_inflce_mat, np.transpose(norm_v_ref))
df = np.multiply(gam_sol, roh_ref*u_ref*db)
cl = np.sum(df)/0.5/roh_ref/u_ref**2/c/b
print(cl/a)


plt.axis('equal')
plt.show()
