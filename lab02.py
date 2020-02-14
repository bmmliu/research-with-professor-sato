# importing the required module

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import math
from scipy.interpolate import make_interp_spline, BSpline
from scipy.ndimage.filters import gaussian_filter1d



# giving a title to my graph
dt = 0.1
pcl = 200  # pacing cycle length (ms)
itr = 10  # the number of beats
tmax = pcl * itr  # total time


# convert to integer
tnmax = int(tmax/dt)
num = 400
# setting the x - coordinates
v_2d = np.zeros((tnmax, num))
v_test_2d = np.zeros((tmax, num))
h_2d = np.zeros((tnmax, num))
f_2d = np.zeros((tnmax, num))
t_set = []
z_set = []
stim = [0] * num
tmp = [0] * num

v_t = []
for i in range(tnmax):
    for j in range(num):
        v_2d[i][j] = 0
        h_2d[i][j] = 1.0
        f_2d[i][j] = 0.9

for i in range(num):
    stim[i] = 0
    z_set.append(i)

pcln = int(pcl/dt)
durn = int(1.0/dt)
for i in range (tnmax):
    t_set.append(i * 0.1)

tn = 0
while (tn < tnmax - 1):

#stimlate the end of the cable
    if tn % pcln < durn:
        for i in range(5):
            stim[i] = 0.3
    else:
        for i in range(5):
            stim[i] = 0
    for i in range(num - 1):
        tauso = 15
        taufi = 0.8
        tauh1 = 4.8
        tauh2 = 10.0
        tausi = 4.0
        tauf1 = 100
        tauf2 = 30
        minf = pow((v_2d[tn][i] / 0.2), 6) / (1 + pow((v_2d[tn][i] / 0.2), 6))
        hinf = 1 / (1 + pow((v_2d[tn][i] / 0.1), 6))
        dinf = pow((v_2d[tn][i] / 0.4), 4) / (1 + pow((v_2d[tn][i] / 0.4), 4))
        finf = 1 / (1 + pow((v_2d[tn][i] / 0.1), 4))
        tauh = tauh1 + tauh2 * math.exp(-20 * pow((v_2d[tn][i] - 0.1), 2))
        tauf = tauf2 + (tauf1 - tauf2) * v_2d[tn][i] * v_2d[tn][i] * v_2d[tn][i]
        jfi = h_2d[tn][i] * minf * (v_2d[tn][i] - 1.3) / taufi
        jsi = f_2d[tn][i] * dinf * (v_2d[tn][i] - 1.4) / tausi
        jso = (1 - math.exp(-4 * v_2d[tn][i])) / tauso
        ion = -(jfi + jsi + jso - stim[i])
        dh = (hinf - h_2d[tn][i]) / tauh
        df = (finf - f_2d[tn][i]) / tauf
        v_2d[tn + 1][i] = v_2d[tn][i] + ion * dt
        h_2d[tn + 1][i] = h_2d[tn][i] + dh * dt
        f_2d[tn + 1][i] = f_2d[tn][i] + df * dt


    dfu = 0.0005
    dx = 0.015
    v_2d[tn + 1][0] = v_2d[tn + 1][2]
    v_2d[tn + 1][num - 1] = v_2d[tn + 1][num - 3]

    i = 1
    while i < num - 1:
        tmp[i] = v_2d[tn + 1][i] + (v_2d[tn + 1][i - 1] + v_2d[tn + 1][i + 1] - 2 * v_2d[tn + 1][i]) * dfu * dt / (dx * dx)
        test1 = tmp[i]
        i = i + 1
    i = 1
    while i < num - 1:
        v_2d[tn + 1][i] = tmp[i]
        test = v_2d[tn + 1][i]
        i = i + 1

    tn = tn + 1

for i in range(tmax):
    for j in range(num):
        v_test_2d[i][j] = v_2d[i * 10][j]

apd1 = []
apd2 = []

startpoint = 7 * pcl
for j in range(num):
    points = []
    i = startpoint - 1
    while (i < 10 * pcl):
        if (v_test_2d[i][j] > 0.1):
            points.append(i)
            break
        i = i + 1
    startpoint = points[0]
    while (i < 10 * pcl):
        if (v_test_2d[i][j] < 0.1):
            points.append(i)
            break
        i = i + 1
    while (i < 10 * pcl):
        if (v_test_2d[i][j] > 0.1):
            points.append(i)
            break
        i = i + 1
    while (i < 10 * pcl):
        if (v_test_2d[i][j] < 0.1):
            points.append(i)
            break
        i = i + 1

    apd1.append(points[1] - points[0] )
    apd2.append(points[3] - points[2])


hf = plt.figure()
ha = hf.add_subplot(111, projection='3d')

X, Y = np.meshgrid(z_set, t_set)

surf = ha.plot_surface(Y,X, v_2d, cmap='PuOr')
ha.set_zlim(0, 1.3)
hf.colorbar(surf, shrink=0.5, aspect=5)
#ha.view_init(azim=0, elev=90)

#plt.xlabel('Space')
#plt.ylabel('Apd')

#plt.title('APD vs Space')
#plt.ylim([70, 130])
#plt.plot(z_set, apd1, label = "apdn")
#plt.plot(z_set, apd2,  label = "apdn+1")
#plt.legend()

plt.show()

