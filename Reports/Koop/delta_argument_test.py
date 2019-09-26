# this test checks if it is possible to determine the mean
# CO offsets of counter-circulating beams (z1 and z2 are of opposite sign)
# if a number of BPMs are deistibuted along the beamline, and each
# can measure the local vertical beam position with a precision SQUID_meas_error
# the beams betatron oscillate vertically with an amplitude ba >> SQUID_meas_error
# The beam CO vertical offsets are determined via model fitting to BPM data.
# Multiple trials are carried out, CO offsets are histogrammed for both beams
# What we look for is whether the <z1> and <z2> are statistically distinguishable

import numpy as np
import matplotlib.pyplot as plt; plt.ion()
from scipy.optimize import curve_fit
from collections import namedtuple
import sys
sys.path.append('/Users/alexaksentyev/REPOS/COSYINF/analysis')
import analysis as ana

plt.ion()
plt.rc('text', usetex=True)
plt.rc('font', size=16, family='normal')
form = lambda ax: ax.ticklabel_format(style='sci', scilimits=(0,0), axis='both')

Pair = namedtuple('Pair', 'one two')
st_type = list(zip(["est","se"], [float]*2))

f = lambda x, a,w,p,d: a*np.sin(w*x + p) + d
l = lambda x, a,b: a + b*x
def fit(function, data, pguess):
    popt, pcov = curve_fit(function, data[0], data[1], p0 = pguess)
    perr = np.sqrt(np.diag(pcov))
    return np.array(list(zip(popt, perr)),dtype=st_type)

ntrl = 100 # number of test trials
ntrn = int(1e0)

# BPM parameters
nbpm = 25*ntrn # number per beamline
SQUID_meas_error = 1e-12 #local BPM measurement error

# betatron oscillation parameters
ba = 1e-6 # amplitude
freq1 = 30 # frequency
freq2 = freq1 + 0.074
ph1 = 0 # phase
ph2 = np.pi/16
z1 = 1e-12 # closed orbit offset
z2 = -1e-12

s = np.tile(np.linspace(0, 150*ntrn, nbpm), (ntrl,1)) # beamline coordinate
err1 = np.random.normal(0, SQUID_meas_error, (ntrl,nbpm)) # measurement error distribution
err2 = np.random.normal(0, SQUID_meas_error, (ntrl,nbpm)) # measurement error distribution
y1 = f(s, ba, freq1, ph1, z1) + err1 # measured vertical position for the primary beam
y2 = f(s, ba, freq2, ph2, z2) + err2 # same for the secondary beam

# fitting vertical beam position measurements with offset-sine function
y1par = ana.fit_matrix(s,y1,f,ini_guess={'a':ba,'w':freq1,'p':ph1,'d':0})
y2par = ana.fit_matrix(s,y2,f,ini_guess={'a':ba,'w':freq2,'p':ph2,'d':0})


zmean = Pair(y1par['d'][:,0], y2par['d'][:,0])
zse = Pair(y1par['d'][:,1], y2par['d'][:,1])

z_diff = zmean.one - zmean.two
z_se = np.sqrt(zse.one**2 + zse.two**2)
test_score = z_diff/z_se

fig, ax = plt.subplots(1,1)
ax.hist(zmean.one, label='one (z1={:4.2e})'.format(z1))
ax.hist(zmean.two, label='two (z2={:4.2e})'.format(z2))
ax.axvline(z1,linestyle='--', color='k')
ax.axvline(z2,linestyle='--', color='k')
ax.legend()
