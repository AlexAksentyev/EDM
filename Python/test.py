from ggplot import ggplot, aes, geom_point, geom_line, theme_bw
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import CBunch
import CSignal
from importlib import reload

b1 = CBunch.Bunch(Npart=1e4)
at = np.arange(0,5, .05*np.pi/b1.Synch["wFreq"])
#at = np.linspace(0,2000, 1e3)

s1 = CSignal.Signal(b1, at)

s1.addNoise(3e-2)
ggplot(s1.Signal, aes(x="Time",y="ValNs")) + geom_point(color="red", size=5) + geom_line(aes(x="Time", y="Val"), color="gray", size=.5) + theme_bw() + geom_point(aes(x="Time",y="Val"), data=s1.specPts, col="blue")

#f = lambda t, l,w,p: 1e3*np.exp(l*t)*np.sin(w*t + p)
#s1.fit(f, [-1e-4, s1.Bunch.Synch["wFreq"], 0])

#s1.Spectrum()
#plt.figure(figsize=(20,20))
#plt.semilogy(s1.PSD.wFreq, s1.PSD.Pow)

s1.findPts()

pts=s1.specPts
