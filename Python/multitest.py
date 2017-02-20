import pathos.multiprocessing as mp



pool = mp.Pool(3)


fn = lambda x,y,x: x**2+y-z

pool.map(fn, list(range(5)), list(range(5)), [10]*5)