#!/usr/bin/env python

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

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
0.989182,
1.58708,
3.4148,
10.9166,
45.645,
41.4099,
81.3343,
192.879,
504.688,
1163.61
])
fjerr = np.array([
0.628493,
0.624047,
0.542491,
1.24241,
6.22873,
2.62773,
4.05763,
6.11079
])

def fpow(x, a, b): return np.exp(b)*pow(x,a)
def fpoly1(x, a, b): return a*x+b
#def fpoly01(x, a): return a*x
def fpoly2(x, a, b, c): return a*x*x + b*x + c
#def fpoly3(x, a, b, c, d): return a*x*x*x + b*x*x + c*x + d

par1, cov1 = curve_fit(fpoly1, np.log(x), np.log(n2y))
#par2, cov2 = curve_fit(fpoly10, x, y)
par3, cov3 = curve_fit(fpoly2, x, n2y)
#par4, cov4 = curve_fit(fpoly3, x, y)

fig, ax = plt.subplots()
ax.errorbar(x, n2y, yerr=n2err, fmt='bo', label='N2jet Java')
ax.errorbar(x, fjy, yerr=fjerr, fmt='mo', label='FastJet')
plt.plot(x, fpow(x, *par1), 'r-',
         label="fit $%.3fN^{%.3f}$" % (np.exp(par1[1]),par1[0]))
#plt.plot(x, fpoly01(x, *par2), 'y-',
#         label="fit $%.3fn$" % (par2[0]))
plt.plot(x, fpoly2(x, *par3), 'g-',
         label="fit $%.3fN^{2}+%.3fN+%.3f$" % (par3[0],par3[1],par3[2]))
#plt.plot(x, fpoly3(x, *par4), 'y-',
#         label="fit $%.3en^{3}+%.3fn^{2}+%.3fn+%.3f$" % (par4[0],par4[1],par4[2],par4[3]))
plt.loglog(basex=2,nonposy='clip')
ax.set_title('Benchmark antikt 0.4 $N$ 100000')
plt.xlabel('number of particles per event, $N$')
plt.ylabel('mean single event clusterization time, mks')

plt.legend(loc=2)

plt.savefig('benchmark.pdf')

