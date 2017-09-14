from importlib  import reload
import matplotlib.pyplot as plt
import CParticle
import CElement



state = [1e-3, -1e-3, 0, 0, 0, 1e-4, 0, 0, 1, 0]
d1 = Drift(500)
p = Particle(state)

p.track([d1],5)

x = [p.fState[i][0] for i in p.fState]
y = [p.fState[i][1] for i in p.fState]
t = [p.fState[i][2] for i in p.fState]
Sx = [p.fState[i][6] for i in p.fState]
Sy = [p.fState[i][7] for i in p.fState]
Ss = [p.fState[i][8] for i in p.fState]
H = [p.fState[i][9] for i in p.fState]


plt.figure(1)
plt.subplot(221); plt.plot(t,x, 'r--'); plt.ylabel('x')
plt.subplot(222); plt.plot(t,y, 'r--'); plt.ylabel('y')
plt.subplot(223); plt.plot(t,Sx, 'r--'); plt.ylabel('Sx')
plt.subplot(224); plt.plot(t,Sy, 'r--'); plt.ylabel('Sy')

