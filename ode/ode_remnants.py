import matplotlib.pyplot
import math
import numpy as np

A = 1
B = 3
C = 3
coefs_runkut = (1 / 6, 1 / 3, 1 / 3, 1 / 6)

def real_value(x):
    return np.array(
        [math.exp(math.sin(x ** 2)), math.exp(B * math.sin(x ** 2)), C * math.sin(x ** 2) + A, math.cos(x ** 2)])

def der_y(num, x, y):
    if num == 0:
        return 2 * x * (y[1]) ** (1 / B) * y[3]
    if num == 1:
        return 2 * B * x * math.exp(B / C * (y[2] - A)) * y[3]
    if num == 2:
        return 2 * C * x * y[3]
    if num == 3:
        return -2 * x * math.log(y[0])
    return 0

def y0():
    return np.array(list(map(float, [1, 1, A, 1])))

def get_k_runkut(num, x, y, h):
    if num == 0:
        return np.array([der_y(i, x, y) for i in range(4)])
    if num == 1:
        return np.array([der_y(i, x + 1 / 2 * h, y + 1 / 2 * h * get_k_runkut(0, x, y, h)) for i in range(4)])
    if num == 2:
        return np.array([der_y(i, x + 1 / 2 * h, y + 1 / 2 * h * get_k_runkut(1, x, y, h)) for i in range(4)])
    if num == 3:
        return np.array([der_y(i, x + h, y + h * get_k_runkut(2, x, y, h)) for i in range(4)])

def runge_kutta4_iteration(x, y, h):
    new_y = y.copy()
    for i in range(4):
        new_y += h * coefs_runkut[i] * get_k_runkut(i, x, y, h)
    return new_y

def runge_kutta4(x_start, y_start, x_end, h):
    x = x_start
    y = y_start.copy()
    while (x <= x_end):
        y = runge_kutta4_iteration(x, y, h)
        x += h
    return y

def get_opt_step(method, x_start, y_start, x_end, method_p, h_start=0.01, tol=0.0001):
    p = method_p
    h = h_start
    yn = method(x_start, y_start, x_end, h)
    yn2 = method(x_start, y_start, x_end, h / 2)
    Rn = (yn2 - yn) / (1 - 2 ** (-p))
    return h * (tol / np.linalg.norm(Rn)) ** (1 / p)

def auto_step(method_iteration, x_start, y_start, x_end, method_p, h_start=0.01, tol=0.001):
    p = method_p
    x_arr = [x_start]
    y = y_start.copy()
    y_arr = []
    h_arr = [h_start]
    delta_arr = []
    while (x_arr[-1] <= x_end):
        h = h_start
        x = x_arr[-1]
        y = method_iteration(x, y, h)
        while (True):
            yh2 = method_iteration(x, y, h / 2)
            yh = method_iteration(x + h / 2, yh2, h / 2)
            rn = (yh - y) / (1 - 2 ** (-p))
            norm_rn = np.linalg.norm(rn)
            if norm_rn > tol * 2 ** p:
                h = h / 2
                y = yh2
            elif tol < norm_rn <= tol * 2 ** p:
                x_arr.append(x + h)
                h = h / 2
                h_arr.append(h)
                y = yh
                y_arr.append(y)
                break
            elif tol / (2 ** (p + 1)) <= norm_rn <= tol:
                x_arr.append(x + h)
                h_arr.append(h)
                y_arr.append(y)
                break
            else:
                h = max(2 * h, max(h_arr))
                y = method_iteration(x, y, h)
        delta_arr.append(math.fabs(np.linalg.norm(real_value(x) - y) / (norm_rn)))
        print(x)
        x += h

    return x_arr, y_arr, h_arr, delta_arr

x_start = 0
x_end = 5
p_runkut = 4
h_start = 0.01
y_start = y0()
h_opt = get_opt_step(runge_kutta4, x_start, y_start, x_end, p_runkut, h_start)
y_opt = runge_kutta4(x_start, y_start, x_end, h_opt)
x_arr, y_arr, h_arr, delta_arr = auto_step(runge_kutta4_iteration, x_start, y_start, x_end,p_runkut, h_start)
print("Оптимальный шаг:", h_opt)
print("Полученное значение:", y_opt)
print("Точное значение:", real_value(x_end))
print("Полученное значение (авто):", y_arr[-1])
matplotlib.pyplot.plot(x_arr, h_arr)
matplotlib.pyplot.show()
matplotlib.pyplot.plot(x_arr[1:], delta_arr)
matplotlib.pyplot.show()
