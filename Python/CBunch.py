import pandas
import numpy

def rnorm(size=1000, mean=0, sd=1):
        """Generates Normally distributed random numbers.
        """
        return(numpy.random.normal(mean, sd, size))

class Bunch:
    """Contains fields and methods of a bunch.
    """
    _G = 700
    _SD = dict(Phi=0, Dgamma=0)
    Synch = dict(wFreq=3, Phi=0)
    EnsPS = None
    
    def __init__(self, Npart=1000, SDdy=1e-3, SDphi=1e-2):
        self._SD["Phi"] = SDphi
        self._SD["Dgamma"] = SDdy
        
        self.EnsPS = pandas.DataFrame({
            "Phi": rnorm(Npart, self.Synch["Phi"], SDphi),
            "wFreq": self.Synch["wFreq"] + self._G*rnorm(Npart, sd=SDdy)**2
        })
        
    def Phase(self, at):
        return(numpy.outer(self.PS.wFreq,at) + numpy.outer(self.PS.Phi,1))
    
    def _Signal(self, at):
        return( numpy.sum( numpy.sin(self.Phase(at)), axis=0) )
        
    def project(self, at):
        self.Pproj = pandas.DataFrame({"Time":at, "Val":self._Signal(at)})


