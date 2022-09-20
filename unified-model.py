from itertools import count
from sqlite3 import Row
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import sqrt

def nmosid(vd, vg, vt0, vdsat, kprime, lmbda):

    dl = 0
    dw = 0
    gamma = 0.4
    vs=0
    vb=0 
    W=1800 
    L=180
    # Fermi potential in p-substrate (page 89)
    phy_f = -0.3

    # See Rabaey, Eq. 3.27, and the explanation.
    # dl and dw model difference of effective (physical) transistor dimensions
    # from drawn dimensions
    Leff = L - dl
    Weff = W - dw

    vgs = vg - vs
    vds = vd - vs
    vsb = vs - vb

    assert vds >= 0 # convention, assumption in modeling equations
    assert vsb > 2 * phy_f # or source-bulk diode is in forward

    # compute threshold voltage, depending on vsb
    vt = vt0 + gamma * (sqrt (abs((-2 * phy_f + vsb))) - sqrt (abs (-2 * phy_f)))
    
    vgt = vgs - vt
    # print(vgs)
    if vgt < 0:
        return 0 # device is off

    # device transconductance 
    k = Weff/Leff * kprime

    # vmin depends on the operating region (lin, sat, velocity sat)
    vmin = min(vgt, vd, vdsat)

    id = k * (vgt * vmin - vmin * vmin / 2) * (1 + lmbda * vds)

    return id

data = pd.read_excel('LTspicedata.xls')

arr = np.array([0, 0, 0, 0, 0, 0])
unifiedarr = np.array([0, 0, 0, 0, 0, 0])

length = 901
i = 1
while i < length:
    row = np.array([data["vds"][i+1], data["Id(M1)"][1 + i], data["Id(M1)"][1 + i + (length+1)], data["Id(M1)"][1 + i + 2*(length+1)], data["Id(M1)"][1 + i + 3*(length+1)], data["Id(M1)"][1 + i + 4*(length+1)]])
    arr = np.vstack([arr, row])
    i += 1

global min_error 
min_error = 9999999999
global error
t_vt0 = 0.22
t_kprime = 0.00016

global min_kprime
global min_lambda
global min_vt0
global min_vdsat

while t_vt0 < 0.24:
    t_vdsat = 0.45
    while t_vdsat < 0.48:
        t_lambda = 0.2
        while t_lambda < 0.250:
            i = 1
            globalerror = []
            globalerror.clear()
            while i < len(arr):
                vd = arr[i][0]
                error = 0
                # print("nmosid ", nmosid(vd, 0.6, t_vt0, t_vdsat, t_kprime, t_lambda), "arr[", i, "][1]", arr[i][1], "error ", error)
                error += abs(nmosid(vd, 0.6, t_vt0, t_vdsat, t_kprime, t_lambda) - arr[i][1])
                # print("nmosid ", nmosid(vd, 0.9, t_vt0, t_vdsat, t_kprime, t_lambda), "arr[", i, "][1]", arr[i][1], "error ", error)
                error += abs(nmosid(vd, 0.9, t_vt0, t_vdsat, t_kprime, t_lambda) - arr[i][2])
                # print("nmosid ", nmosid(vd, 1.2, t_vt0, t_vdsat, t_kprime, t_lambda), "arr[", i, "][1]", arr[i][1], "error ", error)
                error += abs(nmosid(vd, 1.2, t_vt0, t_vdsat, t_kprime, t_lambda) - arr[i][3])
                # print("nmosid ", nmosid(vd, 1.5, t_vt0, t_vdsat, t_kprime, t_lambda), "arr[", i, "][1]", arr[i][1], "error ", error)
                error += abs(nmosid(vd, 1.5, t_vt0, t_vdsat, t_kprime, t_lambda) - arr[i][4])
                # print("nmosid ", nmosid(vd, 1.8, t_vt0, t_vdsat, t_kprime, t_lambda), "arr[", i, "][1]", arr[i][1], "error ", error)
                error += abs(nmosid(vd, 1.8, t_vt0, t_vdsat, t_kprime, t_lambda) - arr[i][5])
                # print("nmosid ", nmosid(vd, 1.8, t_vt0, t_vdsat, t_kprime, t_lambda), "arr[", i, "][1]", arr[i][1], "error ", error)
                i += 1
                # print(error)
                globalerror.append(error)
                #print(sum(globalerror), t_vt0, t_vdsat, t_lambda)
            # print(sum(globalerror), t_vt0, t_vdsat, t_lambda)
            # print(sum(globalerror))
            #print(globalerror)
            
            if(sum(globalerror) < min_error):
                min_kprime = t_kprime
                min_lambda = t_lambda
                min_vt0 = t_vt0
                min_vdsat = t_vdsat
                min_error = sum(globalerror)
                print("min_kprime", min_kprime, "min_lambda", min_lambda,"min_vt0", min_vt0,"min_vdsat", min_vdsat,"min_error", min_error)
            t_lambda += 0.001
        t_vdsat += 0.001
    t_vt0 += 0.001

print("min_kprime ", min_kprime, "\nmin_lambda", min_lambda, "\nmin_vt0", min_vt0, "\nmin_vdsat", min_vdsat)

# vd = arr[i][0]
# row = np.array([vd, nmosid(vd, 0.6), nmosid(vd, 0.9), nmosid(vd, 1.2), nmosid(vd, 1.5), nmosid(vd, 1.8)])
# unifiedarr = np.vstack([unifiedarr, row])


# plot = unifiedarr.transpose()
# plot = plot * 1000
# plot[0] = plot[0] / 1000


# plt.plot (plot[0], plot[1])
# plt.plot (plot[0], plot[2])
# plt.plot (plot[0], plot[3])
# plt.plot (plot[0], plot[4])
# plt.plot (plot[0], plot[5])
# plt.title('$I_d$ against $V_{ds}$ unified model')
# plt.xlabel('$V_{ds}$ (V)')
# plt.ylabel('$I_d$ (mA)')
# plt.xlim(0,1.8)
# plt.ylim(0, 1.5)
# plt.savefig("Id against Vds unified model.png", dpi=600)
# plt.show()
# exit()
