import pandas
import numpy as np
import multiprocessing
from multiprocessing import pool
from scipy.optimize import curve_fit

class Signal:
    """Contains fields and methods of a signal
    
    A signal is a property of a bunch. 
    """
    
    Bunch=None
    Signal=None
    specPts=None
    Model=None
    __NSmpl=None
    
    def __init__(self, bunch, smpl_pts):
        self.Bunch = bunch
        self.Signal = self.Bunch.project(smpl_pts)
        self._NSmpl = 1
        
    def fit(self, FUN, Guess):
        self.Model = curve_fit(FUN, self.Signal.Time, self.Signal.ValNs, p0=Guess)
        return {"Estimate" : self.Model[0], "SE" : np.sqrt(np.diag(self.Model[1]))}
        
    def addNoise(self, rerr):
        self.Signal = self.Signal.assign(Noise = self.Signal.Val*np.random.normal(0, rerr, len(self.Signal)))
        self.Signal = self.Signal.assign(ValNs = self.Signal.Val+self.Signal.Noise)
        
