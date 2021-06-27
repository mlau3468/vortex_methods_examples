import math
from typing import Pattern
import numpy as np
import matplotlib.pyplot as plt

dt = 1
alpha = 5 #angle of attack, degrees
particles = []
U = 1 # freestream velocity
tsteps = 300
c = 1 #chord
rho = 1

class Particle:
    def __init__(self, gamma):
        self.dist = U*dt/2#initial distnace from TE. Place at middle point
        self.gamma = gamma
    def step(self):
        self.dist = self.dist + U*dt

class Panel:
    def __init__(self):
        self.gamma = 0
        self.last_gamma = 0
        self.dgamma_dt = 0
    
    def new_gamma(self, gamma):
        self.last_gamma = gamma
        self.gamma = gamma
        self.dgamma_dt = (self.gamma-self.last_gamma)/dt


#vortex on quarter chord, Kutta condition assumed satisfied.
alpha = np.radians(alpha)
panel = Panel()
gammas = []
lifts = []
for i in range(tsteps):
    # Step particles
    for part in particles:
        part.step()

    # sum wake contributions:
    wake_vel = 0 # small angle, only count in z direction
    wake_gam = 0
    for j in range(len(particles)):
        wake_vel = wake_vel + particles[j].gamma/(2*math.pi*(c/4+( particles[j].dist)))
        wake_gam = wake_gam + particles[j].gamma
    
    RHS = np.zeros((2,1))
    RHS[0,0] = -U*alpha-wake_vel
    RHS[1,0] = -wake_gam
    A = np.zeros((2,2))
    A[0,0] = -1/(2*math.pi*c/2)
    A[0,1] =  1/(2*math.pi*(c/4+U*dt/2))
    A[1,0] = 1
    A[1,1] = 1

    sol = np.linalg.solve(A,RHS)
    panel.new_gamma(sol[0])
    new_part = Particle(sol[1])
    particles.append(new_part)

    # lift per unit span
    lift = rho*(U*panel.gamma+panel.dgamma_dt*c)
    lifts.append(lift)
    gammas.append(panel.gamma)

print('CL: {}'.format(lifts[-1]/c))

plt.figure()
plt.plot(range(tsteps), gammas)
plt.xlabel('time step')
plt.ylabel('panel gamma')

plt.figure()
plt.plot(range(tsteps), lifts)
plt.xlabel('time step')
plt.ylabel('lift per unit span')
plt.show()