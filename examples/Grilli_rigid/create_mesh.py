
import numpy as np


M = 500 #x divisions
N = 10 #y divisions
dx = 0.0125 #grid spacing


x0=0 # initial x coordinate
x = np.arange(x0, x0 + M*dx, dx) 


h0 = 0.33 #water depth
theta = 35*np.pi/180; #slope angle
x_s = 0.6; #slope finish


depth = np.empty((N, M))
depth[:, :] = h0
depth[:, x < x_s] = (-0.0901 + x * np.tan(theta))[x < x_s]
np.savetxt("depth.txt", depth, fmt="%.6f", delimiter=" ")
