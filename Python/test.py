from ggplot import *
import numpy
import CBunch
import CSignal

b1 = CBunch.Bunch(Npart=10, wDist="norm")
at = numpy.arange(0,2000, .95*numpy.pi/b1.Synch["wFreq"])

s1 = CSignal.Signal(b1, at)
ggplot(s1.Signal, aes(x="Time",y="Val")) + geom_line(linetype="dotted", color="blue", size=.5) + theme_bw()
