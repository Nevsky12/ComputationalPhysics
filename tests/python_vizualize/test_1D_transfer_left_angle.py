import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import plot
import matplotlib.animation as animation
import pandas as pd

import os
import subprocess

path = "../../cmake-build-debug/tests/"

proc = subprocess.Popen([path + "./run_test_1D_transfer_left_angle"])
proc.wait()

fig, ax = plt.subplots()
ax = plt.axes(xlim=(0, 20), ylim=(-10, 10))

T = 18
dt1 = 0.3
dt2 = 0.5
dt3 = 0.505
h1 = 0.5
N = 41

df1 = pd.read_csv("CFL1")
df2 = pd.read_csv("CFL2")
df3 = pd.read_csv("CFL3")

x1 = np.array(df1["x"])
y1 = np.array(df1["y"])

x2 = np.array(df2["x"])
y2 = np.array(df2["y"])

x3 = np.array(df3["x"])
y3 = np.array(df3["y"])

line, = ax.plot([], [], lw=3)


def init():
    line.set_data([], [])
    return line,


def animate1(i):
    x = x1[i * N: N + i * N]
    y = y1[i * N: N + i * N]
    print(x)
    print(y)
    line.set_data(x, y)
    return line,

def animate2(i):
    x = x2[i * N: N + i * N]
    y = y2[i * N: N + i * N]
    print(x)
    print(y)
    line.set_data(x, y)
    return line,

def animate3(i):
    x = x3[i * N: N + i * N]
    y = y3[i * N: N + i * N]
    print(x)
    print(y)
    line.set_data(x, y)
    return line,


# create animation using the animate() function
# myAnimation1 = animation.FuncAnimation(fig, animate1, init_func=init, frames=10, \
#                                       interval=100, blit=True, repeat=True)
# myAnimation2 = animation.FuncAnimation(fig, animate2, init_func=init, frames=10, \
#                                        interval=100, blit=True, repeat=True)
myAnimation3 = animation.FuncAnimation(fig, animate3, init_func=init, frames=10, \
                                       interval=100, blit=True, repeat=True)

plt.show()
