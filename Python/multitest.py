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
pts0_l = [pts0.ix[i] for i in range(len(pts0))]

fn = lambda x: s1.Bunch.project(x).Val**2
def finder(e, fn):
    import numpy as np
    x0 = e.Time
    dx = .5
    bnd = x0 + (-dx, dx)
    #x0 = minimize(fn, x0, bounds = ((bnd[0],bnd[1]),)).x
    return fn(x0)

p = mp.Pool(2)

p.map(finder, fn, l)
