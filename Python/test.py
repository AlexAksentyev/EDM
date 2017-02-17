from ggplot import *
import numpy


b1 = Bunch(Npart=1e5)
at = numpy.linspace(0,2000,1e4)
#at = numpy.arange(0,2000, .95*numpy.pi/b1.Synch["wFreq"])
b1.project(at)
ggplot(b1.Pproj, aes(x="Time",y="Val")) + geom_line(linetype="dotted", color="red", size=.5) + theme_bw()