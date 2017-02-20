import pathos.multiprocessing as mp
from scipy.optimize import minimize
import pandas
import scipy.optimize
import numpy as np
import CBunch
import CSignal

b1 = CBunch.Bunch(Npart=10)
s1 = CSignal.Signal(b1, np.linspace(0,5, 10))

pts0 = s1._NullSpecPts()
pts0 = pts0.groupby(list(range(len(pts0))))
