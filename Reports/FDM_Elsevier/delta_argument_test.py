import numpy as np
import matplotlib.pyplot as plt; plt.ion()
from scipy.optimize import curve_fit


d_type = list(zip(["a","b","c","d"], [float]*4))

f = lambda x, a,b,c,d: a*np.sin(b*x + c) + d
l = lambda x, a,b: a + b*x
def fit(function, data, pguess):
    popt, pcov = curve_fit(function, s, data, p0 = pguess)
    perr = np.sqrt(np.diag(pcov))
    return (popt, perr)


n = 250

ba = 1e-6
freq1 = 30
freq2 = 31
ph1 = 0
ph2 = np.pi/16
z1 = 1e-12
z2 = -1e-12
SQUID_meas_error = 1e-12 #local measurement error of the beams relative offsets from the design orbit

s = np.linspace(0, 150, n)
err = np.random.normal(0, SQUID_mean_error, n) 
y1 = f(s, ba, freq1, ph1, z1) + err
y2 = f(s, ba, freq2, ph2, z2) + err
D = y1 - y2

y1par = fit(f,y1,[ba,freq1,ph1,0])
y2par = fit(f,y2,[ba,freq1,ph1,0])
Dpar = fit(l, D, [0, 0])

print(Dpar[0][0], '+-', Dpar[1][0], '(', 2*ba/np.sqrt(n), ')')

if False:
    f, ax = plt.subplots(1,1)
    ax.plot(s, D, 'o-')
    a = Dpar[0][0]
    b = Dpar[0][1]; b=0
    sa = Dpar[1][0]
    ax.plot(s, a + b*s, '--r')
    ax.plot(s, a+sa + b*s, '--r')
    ax.plot(s, a-sa + b*s, '--r')
    ax.set_xlabel('s [m]')
    ax.set_ylabel('BPM measurement')
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
