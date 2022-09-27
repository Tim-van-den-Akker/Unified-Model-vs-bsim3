from itertools import count
from sqlite3 import Row
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import sqrt

data = pd.read_excel('sqrtId_DId.xlsx')
IdArr = np.array([0, 0, 0, 0, 0, 0, 0, 0])
DIdArr = np.array([0, 0, 0, 0, 0, 0, 0, 0])

length = 901
i = 1

spiceName1 = "sqrt(Id(M1))"
spiceName2 = "D(sqrt(Id(M1)))"

while i < length:
    row = np.array([data["vgs"][i+1], data[spiceName1][1 + i], data[spiceName1][1 + i + (length+1)], data[spiceName1][1 + i + 2*(length+1)], data[spiceName1][1 + i + 3*(length+1)], data[spiceName1][1 + i + 4*(length+1)], data[spiceName1][1 + i + 5*(length+1)], data[spiceName1][1 + i + 6*(length+1)]])
    IdArr = np.vstack([IdArr, row])
    i += 1

i = 1
while i < length:
    row = np.array([data["vgs"][i+1], data[spiceName2][1 + i], data[spiceName2][1 + i + (length+1)], data[spiceName2][1 + i + 2*(length+1)], data[spiceName2][1 + i + 3*(length+1)], data[spiceName2][1 + i + 4*(length+1)], data[spiceName2][1 + i + 5*(length+1)], data[spiceName2][1 + i + 5*(length+1)]])
    DIdArr = np.vstack([DIdArr, row])
    i += 1

plot = IdArr.transpose()
plot = plot * 1000
plot[0] = plot[0] / 1000

plotu = DIdArr.transpose()
plotu = plotu * 1000
plotu[0] = plotu[0] / 1000

i = 2
while i < np.shape(plotu)[0] - 5:
    max = np.where(plotu[i] == np.amax(plotu[i]))[0]
    plt.axvline(plotu[0][max], linestyle = 'dotted')
    b = plot[i][max] - (plotu[i][max] * plotu[0][max])
    x = np.array(plot[0])
    n = 1
    y = np.array([0, b])
    while n < len(plot[0]):
        row = x[n], plotu[i][max] * x[n] + b
        y = np.vstack([y, row])
        n += 1
    y = y.transpose()
    plt.plot (plot[0], y[1], '-.')
    print(y)
    i += 1

plt.axvline(0.33579813339955766, linestyle = 'dotted')
# plt.axvline(0.3128452357912393, linestyle = 'dotted')
# plt.axvline(0.28993139829439424, linestyle = 'dotted')
# plt.axvline(0.26768077355900893, linestyle = 'dotted')
# plt.axvline(0.24605056577640688, linestyle = 'dotted')
# plt.axvline(0.22496920689347988, linestyle = 'dotted')
# plt.plot (plot[0], plot[1], label='0V')
plt.plot (plot[0], plot[2], label='0.3V')
# plt.plot (plot[0], plot[3], label='0.6V')
# plt.plot (plot[0], plot[4], label='0.9V')
# plt.plot (plot[0], plot[5], label='1.2V')
# plt.plot (plot[0], plot[6], label='1.5V')
# plt.plot (plot[0], plot[7], label='1.8V')

# plt.plot (plotu[0], plotu[1],'--', label='0V')
plt.plot (plotu[0], plotu[2],'--', label='0.3V')
# plt.plot (plotu[0], plotu[3],'--', label='0.6V')
# plt.plot (plotu[0], plotu[4],'--', label='0.9V')
# plt.plot (plotu[0], plotu[5],'--', label='1.2V')
# plt.plot (plotu[0], plotu[6],'--', label='1.5V')
# plt.plot (plotu[0], plotu[7],'--', label='1.8V')
plt.legend(prop={'size' : 9})
plt.title('$I_d$ against $V_{ds}$ simulated Unified Model vs BSIM3 simulation')
plt.xlabel('$V_{gs}$ (V)')
plt.ylabel('$\sqrt{I_d}$ ($\sqrt{mA}$)')
plt.ylim(0, 40)
plt.savefig("D(sqrt(Id))_sqrt(Id)", dpi=600)
# plt.show()
exit()

