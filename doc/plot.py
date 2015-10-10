#!/usr/bin/env python

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import gridspec

# example data
x = np.array([2,4,8,16,32,64,128,256,512,1024])
n2y = np.array([
3.844791,
6.534534,
11.920509,
22.424794,
44.526984,
80.354958,
169.092560,
397.745738,
1478.529697,
5587.880485
])
n2err = np.array([
19.325180,
26.571852,
24.873007,
30.824656,
46.399137,
67.960798,
89.407254,
130.292613,
278.698871,
787.487320
])

fjy = np.array([
1.06291,
1.67782,
3.53681,
11.5995,
48.4113,
44.0311,
85.8257,
204.854,
538.002,
1236.26
])
fjerr = np.array([
0.858589,
0.851119,
0.98387,
2.2595,
7.44948,
6.00628,
6.24451,
19.1661,
52.984,
86.091
])

def fpoly1(x, a, b): return a*x+b
def fpoly2(x, a, b, c): return a*x*x + b*x + c

par1, cov1 = curve_fit(fpoly1, x[:5], n2y[:5])
par2, cov2 = curve_fit(fpoly2, x[4:], n2y[4:])

fig = plt.figure(figsize=(8, 6))
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0)

ax1 = plt.subplot(gs[0])
ax1.loglog(basex=2,nonposy='clip')
ax1.errorbar(x, n2y, yerr=n2err, fmt='bo', label='N2jet Java')
ax1.errorbar(x, fjy, yerr=fjerr, fmt='ms', label='FastJet C++')
ax1.plot(x[:5], fpoly1(x[:5], *par1), 'r-',
         label="fit $%.1fN+%.1f$" % (np.exp(par1[1]),par1[0]))
ax1.plot(x[4:], fpoly2(x[4:], *par2), 'g-',
         label="fit $%.4fN^{2}+%.2fN+%.0f$" % (par2[0],par2[1],par2[2]))
ax1.set_ylim([5e-1, 1e4])
ax1.set_ylabel('Mean single event clusterization time [mks]')
plt.legend(loc=2)

ax2 = plt.subplot(gs[1], sharex=ax1)
ratio = n2y/fjy
avg = np.mean(ratio)
print(avg)
ax2.scatter(x, ratio)
ax2.set_ylim([0, 5.5])
ax2.set_yticks(np.arange(0, 6, 1))
ax2.axhline(avg, linestyle='--')
ax2.axhline(1, linestyle='--', color='black')
ax2.set_ylabel('$t_\mathtt{N2}/t_\mathtt{FJ}$', fontsize=18, labelpad=15)

plt.xlim([1.7, 1300])
for tick in ax2.xaxis.get_major_ticks():
    tick.label.set_fontsize(14)
plt.xlabel('Number of particles per event, $N$', fontsize=15)
plt.setp(ax1.get_xticklabels(), visible=False)

plt.text(0.5, 0.51, 'average',
         transform=ax2.transAxes, fontsize=14, color='blue',
         verticalalignment='bottom', horizontalalignment='center')

plt.tight_layout()
plt.savefig('benchmark.pdf')
