#!/usr/bin/env python

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# example data
x = np.array([2,4,8,16,32,64,128,256,512,1024])
n2y = np.array([5.650047,10.567299,20.708395,39.470656,82.984448,
                150.885494,241.011962,533.382033,1568.262288,6002.008349])
n2err = np.array([14.063481,16.620472,20.416913,28.153894,40.749753,
                  56.940555,77.336775,113.822303,246.451930,614.800570])

fjy = np.array([1.08526,1.60225,3.40594,10.476,42.1891,43.5207,85.1275,201.4,521.304,1172.13])
fjerr = np.array([0.512875,0.532326,0.526778,1.12427,5.30958,2.11707,3.00389,5.91203,11.1181,17.9807])

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

