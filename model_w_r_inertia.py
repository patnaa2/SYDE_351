from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from math import tan, cos, sin, pi
from scipy.integrate import odeint, simps, cumtrapz 

##############
## y0 = yk 
## y1 = theta
## y2 = px
## y3 = py
##############

def model(y, t):
    yk, theta, vx, vy = y
    
    # constants
    m = 0.68 # Mass of pendulum
    k = 300 # Stiffness of Spring
    Y = 0
    b = 0.5 # Torsional Resistance
    R = 0.8 # Friction
    L = 0.61  # Length of Pendulum
    g = 9.81 # Gravitional acceleration

    # in between terms
    disp = (yk - Y)

    d_yk = vy + ( (L**2) * tan(theta) * vx) 
    d_theta = vx / (L * cos(theta)) 
    d_vy = g + ( -R * d_yk - k * yk)

    # the derivative causality is resolved here, so adding some in between
    # terms for easier debugging
    e_21 = tan(theta) * (d_vy - g) # comes from the left side of bg
    e_24 = d_theta * b  # torsional resistance
    e_22 = d_theta * tan(theta) * vx / (12 * (cos(theta)**2))
    factor = 1 / (1 + (1 / ( 12 * (cos(theta)**2))))

    d_vx = factor * (e_21 - e_22 - e_24)
    return [d_yk, d_theta, d_vx, d_vy]

time = np.linspace(0.0, 10.0, 10000)
yinit = [0.5, pi/8, 0, 0]
y = odeint(model, yinit, time)
ped_y = cumtrapz(y[:,3], time, initial=0)
ped_x = cumtrapz(y[:,2], time, initial=0)
plt.figure(1)

plt.subplot(211)
plt.plot(time, y[:,0], label="yk")
plt.xlabel('t')
plt.ylabel('Displacement of Spring')
plt.legend()

plt.subplot(212)
plt.plot(time, ped_y, label="p_y")
plt.xlabel('t')
plt.ylabel('Displacement of Pendulum in Y')
plt.legend()
plt.show()

plt.figure(2)

plt.subplot(211)
plt.plot(time, ped_x, label="p_x")
plt.xlabel('t')
plt.ylabel('Displacement of Pendulum in X')
plt.legend()

plt.subplot(212)
plt.plot(time, y[:,1], label="theta")
plt.xlabel('t')
plt.ylabel('Angle of rotation')
plt.legend()
plt.show()
