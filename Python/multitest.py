import pathos.multiprocessing as mp



pool = mp.ProcessingPool(3)


fn = lambda x,y,z: x**2+y-z

pool.map(fn, list(range(5)), list(range(5)), [10]*5)
