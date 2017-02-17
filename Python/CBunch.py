import pandas
import numpy
import multiprocessing
from multiprocessing import pool

def rnorm(size=1000, mean=0, sd=1):
        """Generates Normally distributed random numbers.
        """
        return(numpy.random.normal(mean, sd, size))

class Bunch:
    """Contains fields and methods of a bunch.
    """
    _G = 700
    _SD = dict(Phi=0, Dgamma=0)
    def _polProj(self, at):
        return numpy.sum(numpy.sin(self.EnsPS.wFreq*at+self.EnsPS.Phi))
    
    Synch = dict(wFreq=3, Phi=0)
    EnsPS = None
    
    def __init__(self, Npart=1000, SDdy=1e-3, SDphi=1e-2, wDist="phys"):
        self._SD["Phi"] = SDphi
        self._SD["Dgamma"] = SDdy
        
        switch = lambda x: {
                "phys" : self.Synch["wFreq"] + self._G*rnorm(Npart, sd=SDdy)**2,
                "norm" : rnorm(Npart, self.Synch["wFreq"], self._G*SDdy**2)
            }.get(x, None)
        
        self.EnsPS = pandas.DataFrame({
            "Phi": rnorm(Npart, self.Synch["Phi"], SDphi),
            "wFreq": switch(wDist)
        })
        
    def Phase(self, at):
        pass
        #return(numpy.outer(self.PS.wFreq,at) + numpy.outer(self.PS.Phi,1))
        
    def project(self, at):
        p = pool.Pool(processes=multiprocessing.cpu_count())
        res = pandas.DataFrame({"Time":at, "Val":p.map(self._polProj,at)})
        p.close()
        return res


