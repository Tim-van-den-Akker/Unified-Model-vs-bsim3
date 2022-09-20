from itertools import count
from sqlite3 import Row
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import sqrt

data = pd.read_excel('LTspicedata.xls')

arr = np.array([0, 0, 0, 0, 0, 0])

length = 901
i = 1
while i < length:
    row = np.array([data["vds"][i+1], data["Id(M1)"][1 + i], data["Id(M1)"][1 + i + (length+1)], data["Id(M1)"][1 + i + 2*(length+1)], data["Id(M1)"][1 + i + 3*(length+1)], data["Id(M1)"][1 + i + 4*(length+1)]])
    arr = np.vstack([arr, row])
    i += 1

plot = arr.transpose()
plot = plot * 1000
plot[0] = plot[0] / 1000


plt.plot (plot[0], plot[1])
plt.plot (plot[0], plot[2])
plt.plot (plot[0], plot[3])
plt.plot (plot[0], plot[4])
plt.plot (plot[0], plot[5])
plt.title('$I_d$ against $V_{ds}$')
plt.xlabel('$V_{ds}$ (V)')
plt.ylabel('$I_d$ (mA)')
plt.savefig("Id against Vds.png", dpi=600)
plt.show()
exit()
