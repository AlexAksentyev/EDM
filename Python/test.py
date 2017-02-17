from ggplot import *
import numpy as np
import matplotlib.pyplot as plt
import CBunch
import CSignal

b1 = CBunch.Bunch(Npart=3, wDist="norm")
at = numpy.arange(0,2000, .95*np.pi/b1.Synch["wFreq"])
#at = np.linspace(0,2000, 1e3)

s1 = CSignal.Signal(b1, at)


s1.addNoise(3e-2)
ggplot(s1.Signal, aes(x="Time",y="ValNs")) + geom_point(color="red", size=5) + geom_line(aes(x="Time", y="Val"), color="gray", size=.5) + theme_bw()

f = lambda t, l,w,p: 1e3*np.exp(l*t)*np.sin(w*t + p)
s1.fit(f, [-1e-4, s1.Bunch.Synch["wFreq"], 0])

s1.Spectrum()

plt.semilogy(s1.PSD.wFreq, s1.PSD.Pow)

