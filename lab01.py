# importing the required module
import matplotlib.pyplot as plt
import numpy as np
import math


# naming the x axis
plt.xlabel('Potential (mV)')
# giving a title to my graph


# setting the x - coordinates
arr_di = []
arr_n = []
arr_n1 = []
pcl_set = []

# initial value
v = 0.0  # membrane voltage
h = 1.0  # Na channel inactivation gate
f = 0.9  # Ca channel inactivation gate
# constants
tauso = 15.0
taufi = 0.8
tauh1 = 4.8
tauh2 = 10.0

tausi = 4.0
tauf1 = 100.0
tauf2 = 30.0

dt = 0.1
tn = 0.0

pcl = 200  # pacing cycle length (ms)
itr = 10  # the number of beats
tmax = pcl * itr  # total time
t = np.arange(0, tmax, 0.1)
# convert to integer
tnmax = int(tmax/dt)
pcln = int(pcl/dt)
durn = int(1.0/dt)

def get_k (xi, yi):
    stim = 0
    if (((xi * 10) % pcln) < durn):
        stim = 0.3

    minf = (math.pow((yi / 0.2), 6)) / (1 + (math.pow((yi / 0.2), 6)))
    hinf = 1 / (1 + (math.pow((yi / 0.1), 6)))
    dinf = (math.pow((yi / 0.4), 4)) / (1 + (math.pow((yi / 0.4), 4)))
    finf = 1 / (1 + (math.pow((yi / 0.1), 4)))

    tauh = tauh1 + tauh2 * math.exp(-20 * (math.pow((yi - 0.1), 2)))
    tauf = tauf2 + (tauf1 - tauf2) * math.pow(yi, 3)

    jfi = h * minf * (yi - 1.3) / taufi
    jsi = f * dinf * (yi - 1.4) / tausi
    jso = (1 - math.exp(-4 * yi)) / tauso
    ion = -(jfi + jsi + jso - stim)
    return ion


def get_h_k(x1 ,x2):
    hinf = 1 / (1 + (math.pow((x1 / 0.1), 6)))
    tauh = tauh1 + tauh2 * math.exp(-20 * (math.pow((x1 - 0.1), 2)))
    hk = (hinf - x2) / tauh
    return hk

def get_f_k(x1 ,x2):
    finf = 1 / (1 + (math.pow((x1 / 0.1), 4)))
    tauf = tauf2 + (tauf1 - tauf2) * math.pow(x1, 3)
    fk = (finf - x2) / tauf
    return fk



def get_graph():
    global v, h, f, tn, dt
    stim = 0
    if ((tn % pcln) < durn):
        stim = 0.3
    minf = (math.pow((v / 0.2), 6)) / (1 + (math.pow((v / 0.2), 6)))
    hinf = 1 / (1 + (math.pow((v / 0.1), 6)))
    dinf = (math.pow((v / 0.4),  4)) / (1 + (math.pow((v / 0.4), 4)))
    finf = 1 / (1 + (math.pow((v / 0.1), 4)))

    tauh = tauh1 + tauh2 * math.exp(-20 * (math.pow((v - 0.1), 2)))
    tauf = tauf2 + (tauf1 - tauf2) * math.pow(v, 3)

    jfi = h * minf * (v - 1.3) / taufi
    jsi = f * dinf * (v - 1.4) / tausi
    jso = (1 - math.exp(-4 * v)) / tauso
    ion = -(jfi + jsi + jso - stim)
    kv1 = ion
    kv2 = get_k((tn / 10 + 0.05), v + (kv1 / 2))
    kv3 = get_k((tn / 10 + 0.05), v + (kv2 / 2))
    kv4 = get_k((tn / 10 + 0.1), v + kv3)

    v = v + dt / 6 * (kv1 + 2 * kv2 + 2 * kv3 + kv4)
    # update variables
    kh1 = (hinf - h) / tauh
    v_23 = v - dt / 6 * (kv1 + 2 * kv2 + 2 * kv3 + kv4) + 0.05 * kv2
    kh2 = get_h_k(v_23, h + (kh1 / 2))
    kh3 = get_h_k(v_23, h + (kh2 / 2))
    kh4 = get_h_k(v, h + kh3)
    h = h + dt / 6 * (kh1 + 2 * kh2 + 2 * kh3 + kh4)

    kf1 = (finf - f) / tauf
    v_23 = v - dt / 6 * (kv1 + 2 * kv2 + 2 * kv3 + kv4) + 0.05 * kv2
    kf2 = get_f_k(v_23, f + (kf1 / 2))
    kf3 = get_f_k(v_23, f + (kf2 / 2))
    kf4 = get_f_k(v, f + kf3)
    f = f + dt / 6 * (kf1 + 2 * kf2 + 2 * kf3 + kf4)
    #f = f + dt * (finf - f) / tauf
    tn = tn + 1







def get_graph_v_t(t):

    get_graph()
    return v

def get_graph_h_t(t):

    get_graph()
    return h

def get_graph_f_t(t):

    get_graph()

    return f

def get_apd(pcl_l):
    global v, h, f, tn, tmax, pcln, pcl
    pcl = pcl_l
    tmax = pcl * itr  # total time
    pcln = int(pcl / dt)
    f1 = t / t * 0.1

    fig = plt.figure()
    line1 = fig.add_subplot(1, 1, 1)
    line1.plot(t, f1)
    v_graph = np.vectorize(get_graph_v_t)
    line2 = fig.add_subplot(1, 1, 1)
    line2.plot(t, v_graph(t))
    idx = np.argwhere(np.diff(np.sign(f1 - v_graph(t)))).flatten()

    i0 = idx[2] - idx[1]
    i1 = idx[3] - idx[2]
    i2 = idx[4] - idx[3]
    i3 = idx[5] - idx[4]
    a = [i0, i1, i2, i3]
    b = sorted(a)

    apd1 = b[3] / 10
    if (pcl < 139.7):
        apd2 = b[1] / 10
    else:
        apd2 = b[2] / 10
    di = pcl - apd1
    arr_di.append(di)
    arr_n.append(apd1)
    arr_n1.append(apd2)
    pcl_set.append(pcl)

    plt.close(fig)

    v = 0.0  # membrane voltage
    h = 1.0  # Na channel inactivation gate
    f = 0.9
    tn = 0.0



def get_di(pcl_l):
    global v, h, f, tn, tmax, pcln, pcl, arr_di, arr_n1
    pcl = pcl_l
    tmax = pcl * itr  # total time
    pcln = int(pcl / dt)
    f1 = t / t * 0.1

    fig = plt.figure()
    line1 = fig.add_subplot(1, 1, 1)
    line1.plot(t, f1)
    v_graph = np.vectorize(get_graph_v_t)
    line2 = fig.add_subplot(1, 1, 1)
    line2.plot(t, v_graph(t))
    idx = np.argwhere(np.diff(np.sign(f1 - v_graph(t)))).flatten()
    #plt.show()

    i0 = idx[2] - idx[1]
    i1 = idx[3] - idx[2]
    i2 = idx[4] - idx[3]
    i3 = idx[5] - idx[4]
    a = [i0, i1, i2, i3]
    b = sorted(a)

    apd1 = b[3] / 10
    if (pcl < 139.7):
        apd2 = b[1] / 10
    else:
        apd2 = b[2] / 10
    di = pcl - apd1
    arr_di.append(di)
    arr_n1.append(apd2)
    plt.close(fig)

    v = 0.0  # membrane voltage
    h = 1.0  # Na channel inactivation gate
    f = 0.9
    tn = 0.0
    return apd2

def main():
    mode = int(input('Choose mode\n'))
    global pcl, arr_n1, arr_di
    if mode == 1:
        plt.xlabel('Time')
        plt.ylabel('V')
        plt.title('V vs Time')
        v_graph = np.vectorize(get_graph_v_t)
        plt.plot(t, v_graph(t))
        plt.show()

    if mode == 2:
        plt.xlabel('Time')
        plt.ylabel('V')
        plt.title('H vs Time')
        h_graph = np.vectorize(get_graph_h_t)
        plt.plot(t, h_graph(t))
        plt.show()
    if mode == 3:
        plt.xlabel('Time')
        plt.ylabel('V')
        plt.title('F vs Time')
        f_graph = np.vectorize(get_graph_f_t)
        plt.plot(t, f_graph(t))
        plt.show()
    if mode == 4:
        plt.xlabel('PCL')
        plt.ylabel('APD')
        plt.title('APD vs PCL')
        pcl = np.arange(135.0, 165, 0.1)
        pcl_apd1_graph = np.vectorize(get_apd)
        plt.plot(pcl, pcl_apd1_graph(pcl))
        plt.clf()
        plt.plot(pcl_set, arr_n)
        plt.plot(pcl_set, arr_n1)

        plt.show()
    if mode == 5:
        plt.xlabel('APD(n+1)')
        plt.ylabel('DI(n)')
        plt.title('APD(n+1) vs DI(n)')
        pcl = np.arange(135.0, 165, 0.1)
        pcl_di_graph = np.vectorize(get_di)
        plt.plot(pcl, pcl_di_graph(pcl))
        plt.clf()
        plt.plot(arr_di, arr_n1)

        plt.show()

main()



