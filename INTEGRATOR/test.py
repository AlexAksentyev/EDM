from importlib  import reload
from ggplot import ggplot, aes, geom_line, theme_bw, facet_wrap
import pandas as PDS
import CParticle as PCL
import CElement as ENT

theme_bw()

state = [1e-3, -1e-3, 0, 0, 0, 1e-4, 0, 0, 1, 0]
ds_25 = ENT.Drift(.25)
ds_15 = ENT.Drift(.15)
ds2_2 = ENT.Drift(2.2)
fquad = ENT.MQuad(5,.831)
dquad = ENT.MQuad(5,-.86)
fsext = ENT.MSext(5,1.023)
dsext = ENT.MSext(5,-1.34)

p = PCL.Particle(state)

FODO = [fquad, ds_25, dquad, ds_25]

R = 7.55
B0 = p.Pc(p.fKinEn0)*1e6/R
mdip = ENT.MDipole(1.8, R, (0,B0,0))
WA = ENT.Wien(1.808,R,.05,120e5,B0)


p.track([mdip]*3,40,FWD=False)

x = [p.fStateLog[i][0] for i in p.fStateLog]
y = [p.fStateLog[i][1] for i in p.fStateLog]
t = [p.fStateLog[i][2] for i in p.fStateLog]
dW = [p.fStateLog[i][5] for i in p.fStateLog]
Sx = [p.fStateLog[i][6] for i in p.fStateLog]
Sy = [p.fStateLog[i][7] for i in p.fStateLog]
Ss = [p.fStateLog[i][8] for i in p.fStateLog]
H = [p.fStateLog[i][9] for i in p.fStateLog]

df = PDS.DataFrame({'x':x,'y':y,'Sx':Sx,'Sy':Sy,'t':t,'H':H,'dW':dW})
df = PDS.melt(df, id_vars=['t','H'])
ggplot(df.loc[df['variable'].isin(['x','y','Sx','Sy'])],aes(x='t',y='value'))+\
    geom_line() + facet_wrap('variable',scales='free')


