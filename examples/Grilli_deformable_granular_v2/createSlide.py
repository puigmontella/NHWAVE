import numpy as np

M = 600
N = 10
dx = 0.01

SlideT = 0.084
SlideL = 0.1213
SlideX0 = 0.1987

# cell-center coordinates (matches typical xc(i))
x = (np.arange(M) + 0.5) * dx

x1 = SlideX0
x3 = SlideX0 + SlideL
x2 = x1 + 0.95 * (x3 - x1)

zt = np.zeros(M)

m1 = (x >= x1) & (x <= x2)
zt[m1] = SlideT * (x[m1] - x1) / (x2 - x1)

m2 = (x > x2) & (x <= x3)
zt[m2] = SlideT * (x3 - x[m2]) / (x3 - x2)

slide = np.tile(zt, (N, 1))  # uniform in y

np.savetxt("slide_combine1.txt", slide, fmt="%.6f", delimiter=" ")
print("Wrote slide_combine1.txt:", slide.shape, "max=", slide.max())
