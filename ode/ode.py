import matplotlib.pyplot
from mpl_toolkits.mplot3d import Axes3D

def X(x, y, s):
    return s * (y - x)
def Y(x, y, z, r):
    return (r - z) * x - y
def Z(x, y, z, b):
    return x * y - b * z

def runge_kutta1(x, y, z, h):
    x_new = x
    y_new = y
    z_new = z
    x_new += X(x, y, s) * h
    y_new += Y(x, y, z, r) * h
    z_new += Z(x, y, z, b) * h

    return (x_new, y_new, z_new)

def runge_kutta2(x, y, z, h):
    x_new = x
    y_new = y
    z_new = z
    x_new = X(x, y, s)
    y_new = Y(x, y, z, r)
    z_new = Z(x, y, z, b)

    k2 = X((x + x_new * h * c2), (y + y_new * h * c2), s)
    l2 = Y((x + x_new * h * c2), (y + y_new * h * c2), (z + z_new * h * c2), r)
    m2 = Z((x + x_new * h * c2), (y + y_new * h * c2), (z + z_new * h * c2), b)

    x += (x_new * b1 + b2 * k2) * h
    y += (y_new * b1 + b2 * l2) * h
    z += (z_new * b1 + b2 * m2) * h

    return (x, y, z)

s, r, b = 10, 28, 8 / 3
x_rk = [1]
y_rk = [1]
z_rk = [1]
x_runge_kutta1 = [1]
y_runge_kutta1 = [1]
z_runge_kutta1 = [1]
c2 = 0.5
b2 = 1 / (2 * c2)
b1 = 1 - 1 / (2 * c2)
h = 0.01

for i in range(1000):
    x = x_runge_kutta1[i]
    y = y_runge_kutta1[i]
    z = z_runge_kutta1[i]
    position = runge_kutta1(x, y, z, h)
    x_runge_kutta1.append(position[0])
    y_runge_kutta1.append(position[1])
    z_runge_kutta1.append(position[2])
    
for i in range(1000):
    x = x_rk[i]
    y = y_rk[i]
    z = z_rk[i]
    position = runge_kutta2(x, y, z, h)
    x_rk.append(position[0])
    y_rk.append(position[1])
    z_rk.append(position[2])

fig = matplotlib.pyplot.figure()
ax_rk = fig.gca(projection='3d')
ax_eu = fig.gca(projection='3d')
ax_eu.plot(x_runge_kutta1, y_runge_kutta1, z_runge_kutta1)
ax_rk.plot(x_rk, y_rk, z_rk)

matplotlib.pyplot.show()
