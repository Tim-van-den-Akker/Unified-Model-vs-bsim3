from itertools import count
from sqlite3 import Row
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import sqrt

def nmosid(vd, vg, vt0 = 0.237, vdsat = 0.468, kprime = 0.0001587, lmbda = 0.261):

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

i = 0
while i < len(arr):
    vd = arr[i][0]
    row = np.array([vd, nmosid(vd, 0.6), nmosid(vd, 0.9), nmosid(vd, 1.2), nmosid(vd, 1.5), nmosid(vd, 1.8)])
    unifiedarr = np.vstack([unifiedarr, row])
    i += 1


plot = arr.transpose()
plot = plot * 1000
plot[0] = plot[0] / 1000

plotu = unifiedarr.transpose()
plotu = plotu * 1000
plotu[0] = plotu[0] / 1000


plt.plot (plot[0], plot[1], label='0.6V BSIM')
plt.plot (plot[0], plot[2], label='0.9V BSIM')
plt.plot (plot[0], plot[3], label='1.2V BSIM')
plt.plot (plot[0], plot[4], label='1.5V BSIM')
plt.plot (plot[0], plot[5], label='1.8V BSIM')

plt.plot (plotu[0], plotu[1],'--', label='0.6V UM')
plt.plot (plotu[0], plotu[2],'--', label='0.9V UM')
plt.plot (plotu[0], plotu[3],'--', label='1.2V UM')
plt.plot (plotu[0], plotu[4],'--', label='1.5V UM')
plt.plot (plotu[0], plotu[5],'--', label='1.8V UM')
plt.legend(prop={'size' : 9})
plt.title('$I_d$ against $V_{ds}$ simulated Unified Model vs BSIM3 simulation')
plt.xlabel('$V_{ds}$ (V)')
plt.ylabel('$I_d$ (mA)')
plt.xlim(0,1.8)
plt.ylim(0, 1.5)
plt.savefig("Id against Vds unified model.png", dpi=600)
plt.show()
exit()