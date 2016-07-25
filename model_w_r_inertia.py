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
    m = 0.408 # Mass of pendulum
    k = 100 # Stiffness of Spring
    b = 1.2 # Torsional Resistance
    R = 2 # Friction
    L = 0.61  # Length of Pendulum
    g = 9.81 # Gravitional acceleration
    Y = 0.5 

    # in between terms
    disp = (yk - Y)

    d_yk = vy + ((tan(theta) * vx)) 
    d_theta = vx / (L * cos(theta)) 
    d_vy = g + (( -R * d_yk - k * yk)/m)

    # the derivative causality is resolved here, so adding some in between
    # terms for easier debugging
    e_21 = tan(theta) * (d_vy - g) # comes from the left side of bg
    e_24 = d_theta * b  # torsional resistance
    e_22 = d_theta * tan(theta) * vx / (12 * (cos(theta)**2))
    factor = 1 / (1 + (1 / ( 12 * (cos(theta)**2))))

    d_vx = factor * (e_21 - e_22 - e_24)
    return [d_yk, d_theta, d_vx, d_vy]

time = np.linspace(0.0, 10.0, 10000)
yinit = [0, pi/2, 0, 0]
y = odeint(model, yinit, time)

# the state equations give us velocity
# integrate again to get displacement
# our variable of interest
ped_y = cumtrapz(y[:,3], time, initial=0)
ped_x = cumtrapz(y[:,2], time, initial=0)

plt.figure(1)

plt.subplot(211)
plt.plot(time, y[:,0], label="yk")
plt.xlabel('t')
plt.ylabel('Displacement of Spring')
plt.legend()

plt.subplot(212)
plt.plot(time, y[:,1], label="theta")
plt.xlabel('t')
plt.ylabel('Angle of rotation')
plt.legend()
plt.show()

plt.figure(2)

plt.subplot(211)
plt.plot(time, ped_x, label="yx")
plt.xlabel('t')
plt.ylabel('Displacement of Pendulum in X')
plt.legend()

plt.subplot(212)
plt.plot(time, ped_y, label="yp")
plt.xlabel('t')
plt.ylabel('Displacement of Pendulum in Y')
plt.legend()
plt.show()
