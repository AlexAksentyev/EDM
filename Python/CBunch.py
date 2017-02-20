import pandas
import numpy as np
import multiprocessing as mp

def rnorm(size=1000, mean=0, sd=1):
        """Generates Normally distributed random numbers.
        """
        return(np.random.normal(mean, sd, size))

class Bunch:
    """Contains fields and methods of a bunch.
    """
    _G = 700
    _SD = dict(Phi=0, Dgamma=0)
    def _polProj(self, at):
        #p = mp.ProcessingPool(mp.cpu_count())
        #def fn(sub, x):
        #    import numpy as np
        #    return np.sin(sub.wFreq*x + sub.Phi)
        #    
        #df_g = self.EnsPS.groupby(list(range(len(self.EnsPS))))
        #sin_lst = p.map(fn, [group for name, group in df_g], [at]*len(df_g))
        #p.close()
        #p.join()
        #return np.sum(sin_lst)
        return np.sum(np.sin(self.EnsPS.wFreq*at+self.EnsPS.Phi))
    
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
        p = mp.Pool(processes=mp.cpu_count())
        res = pandas.DataFrame({"Time":at, "Val":p.map(self._polProj, at)})
        p.close()
        p.join()
        return res


