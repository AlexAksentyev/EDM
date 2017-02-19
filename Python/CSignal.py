import pandas
import numpy as np
import multiprocessing
import pathos.multiprocessing as mp
from scipy.optimize import curve_fit
from scipy.signal import welch

class Signal:
    """Contains fields and methods of a signal
    
    A signal is a property of a bunch. 
    """
    
    Bunch=None
    Signal=None
    specPts=None
    Model=None
    PSD=None
    __NSmpl=None
    
    def __init__(self, bunch, smpl_pts):
        self.Bunch = bunch
        self.Signal = self.Bunch.project(smpl_pts)
        self._NSmpl = 1
        
    def split(self):
        pass
    
    def sample(self, at):
        pass
    
    def addNoise(self, rerr):
        self.Signal = self.Signal.assign(Noise = self.Signal.Val*np.random.normal(0, rerr, len(self.Signal)))
        self.Signal = self.Signal.assign(ValNs = self.Signal.Val+self.Signal.Noise)
        
    def fit(self, FUN, Guess):
        self.Model = curve_fit(FUN, self.Signal.Time, self.Signal.ValNs, p0=Guess)
        return {"Estimate" : self.Model[0], "SE" : np.sqrt(np.diag(self.Model[1]))}
        
    def Spectrum(self):
        fs = 1/(self.Signal.Time[2] - self.Signal.Time[1])
        psd = welch(self.Signal.ValNs, fs)
        self.PSD = pandas.DataFrame({
            "Freq": psd[0],
            "Pow" : psd[1]
        })
        self.PSD = self.PSD.assign(wFreq = self.PSD.Freq*2*np.pi)

    def findPts(self, what="Node", wguess=None):
        
        switch = lambda x: {"Node": 1, "Envelope": -1}.get(x, 1)
        
        k = switch(what) #this'll determine the direction of optimization
        fn = lambda x: k*self.Bunch.project(x).Val**2
        def finder(e): # work on here
            from scipy.optimize import minimize
            x0 = e.Time
            dx = np.pi/self.Bunch.Synch["wFreq"]/2
            bnd = x0 + (-dx, +dx)
            x0 = minimize(fn, x0, bounds=((bnd[0], bnd[1]),)).x
            #df1 = pandas.DataFrame({
            #    "N": pts0.ix[i].N,
            #    "Side": pts0.ix[i].Side,
            #    "Which": "Optim"
            #}, index = [i])
            #return pandas.concat([df1, self.Bunch.project(x0)], axis=1, join_axes=[df1.index])
            return x0
      
        pts0 = self._NullSpecPts(what, wguess)
        pts0_l = [pts0.ix[i] for i in range(len(pts0))]
        pool = mp.Pool(processes=multiprocessing.cpu_count())
        pts1 = pool.map(finder, pts0_l)
        
        p.close()
        p.join()
      
        #self.specPts = [pts0, pts1].sort_values("Time")
        return pts1
        
    def _NullSpecPts(self, what="Node", wref=None):
        if wref is not None:
            w0 = wref 
        else: 
            w0 = self.Bunch.Synch["wFreq"]
        
        p0 = self.Bunch.Synch["Phi"]
        
        _dum = lambda Time: int(np.floor((w0*Time+p0-np.pi/2)/2/np.pi))
      
        Nstt = _dum(self.Signal.Time[0])
        Ntot = _dum(self.Signal.Time[len(self.Signal)-1])
      
        switch = lambda x: {"Node": 0, "Envelope": np.pi/2}.get(x, 0)
      
        d = switch(what)
        tnu = [(2*np.pi*x -p0+d)/w0  for x in range(Nstt,Ntot+1)]
        tnu = [x for x in tnu if not x < 0]
        tnd = [x+np.pi/w0 for x in tnu]
      
        df1 = pandas.DataFrame({
            "N" : list(range(1,len(tnu)+1)) + list(range(1,len(tnd)+1)),
            "Side": ["U"]*len(tnu)+ ["D"]*len(tnd),
            "Which": "Null"
        })
        
        return pandas.concat([df1, self.Bunch.project(tnu+tnd)], axis=1)
