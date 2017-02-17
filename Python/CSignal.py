import pandas
import numpy
import multiprocessing
from multiprocessing import pool


class Signal:
    """Contains fields and methods of a signal
    
    A signal is a property of a bunch. 
    """
    
    Bunch=None
    Signal=None
    specPts=None
    ModelCoef=None
    __NSmpl=None
    
    def __init__(self, bunch, smpl_pts):
        self.Bunch = bunch
        self.Signal = self.Bunch.project(smpl_pts)
        self._NSmpl = 1
        

