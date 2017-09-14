from CParticle import *
from CElement import *
from importlib  import reload
from ggplot import *
import pandas

theme_bw()

state = [1e-3, -1e-3, 0, 0, 0, 1e-4, 0, 0, 1, 0]
ds = Drift(.25)
fquad = MQuad(5,.831)
dquad = MQuad(5,-.86)
fsext = MSext(5,1.023)
dsext = MSext(5,-1.34)
p = Particle(state)

p.track([fquad, ds, dquad],500)

x = [p.fState[i][0] for i in p.fState]
y = [p.fState[i][1] for i in p.fState]
t = [p.fState[i][2] for i in p.fState]
dW = [p.fState[i][5] for i in p.fState]
Sx = [p.fState[i][6] for i in p.fState]
Sy = [p.fState[i][7] for i in p.fState]
Ss = [p.fState[i][8] for i in p.fState]
H = [p.fState[i][9] for i in p.fState]

df = pandas.DataFrame({'x':x,'y':y,'Sx':Sx,'Sy':Sy,'t':t,'H':H,'dW':dW})
df = pandas.melt(df, id_vars=['t','H'])
ggplot(df.loc[df['variable'].isin(['x','y','Sx','Sy'])],aes(x='t',y='value'))+\
    geom_line() + facet_wrap('variable',scales='free')



