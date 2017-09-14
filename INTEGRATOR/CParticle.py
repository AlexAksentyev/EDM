from scipy.integrate import odeint
import numpy as np



class Particle:
    
    _Stats = {}
    
    _EZERO = 1.602176462e-19 # Coulomb
    _clight = 2.99792458e8 # m/s
    
    _fIniState = None
    fState = None 
    
    fMass0 = 1876.5592 # deuteron mass in MeV
    fKinEn0 = 270.005183 # deuteron magic energy
    fG = -.142987
    
    fGamma0 = None # reference particle's
    fBeta0 = None  # gamma, beta
    
    def __init__(self, State0):
        
        self._fIniState = list(State0)
        self.fState = {0:list(State0)}
        
        self.fGamma0, self.fBeta0 = self._GammaBeta(self.fKinEn0)
    
    def _GammaBeta(self, NRG):
        gamma = NRG / self.fMass0 + 1
        beta = np.sqrt(gamma**2-1)/gamma
        return (gamma, beta)
    
    def _RHS(self, state, at, element):
        x,y,s,px,py,dEn,Sx,Sy,Ss,H = state.values() # px, py are normalized to P0c for consistency with the other vars
        
        KinEn = self.fKinEn0*(1+dEn) # dEn = (En - En0) / En0
        lPC = lambda KNRG:  np.sqrt((self.fMass0 + KNRG)**2 - self.fMass0**2)
        
        Pc = lPC(KinEn) # momentum
        P0c = lPC(self.fKinEn0) # reference momentum
        
        Px,Py = [P0c*x for x in (px,py)] # turn px,py back to MeVs
        Ps = np.sqrt(Pc**2 - Px**2 - Py**2)
        
        Ex,Ey,Es = element.EField(state)
        Bx,By,Bs = element.BField(state)
        
        kappa = element.fCurve
        hs = 1 + kappa*x # don't exactly understand this
        H = Pc*hs/Ps # or this
        
        xp,yp = [x * hs/Ps for x in (Px,Py)]
        Wp = (Ex*xp +Ey*yp +Es) * 1e-6 # Kinetic energy prime (in MeV)
        gammap = Wp/self.fMass0 # gamma prime
        
        gamma,beta = self._GammaBeta(KinEn)
        v = beta*self._clight
        
        m0 = self.fMass0/self._clight**2
        q = self._EZERO
        
        ## I don't understand the following formulas
        betap = (Wp*(self.fMass0)**2)/((KinEn+self.fMass0)**2*np.sqrt(KinEn**2+2*KinEn*self.fMass0))
        tp = H/v
        
        D = (q/(m0*hs))*(xp*By-yp*Bx+H*Es/v)-((gamma*v)/(H*hs))*3*kappa*xp # what's this?
        xpp=((-H*D)/(gamma*v))*xp+(H/(Pc*1e6))*(H*Ex/v+yp*Bs-hs*By)+kappa*hs
        ypp=((-H*D)/(gamma*v))*yp+(H/(Pc*1e6))*(H*Ey/v+hs*Bx-xp*Bs)
        
        Pxp = Px*(betap/beta - gammap/gamma)+Pc*xpp/H-Px*((Px*xpp)/(Pc*H)+(Py*ypp)/(Pc*H)+(hs*kappa*xp)/(H**2))
        Pyp = Py*(betap/beta - gammap/gamma)+Pc*ypp/H-Py*((Px*xpp)/(Pc*H)+(Py*ypp)/(Pc*H)+(hs*kappa*xp)/(H**2))
        
        t5 = tp
        t6 =  t5* (q / (gamma * m0 * self.fMass0)) * (self.fG + 1/(1 + gamma))
        sp1 = t5*(-q / (gamma*m0))*(1 + self.fG * gamma)
        sp2 = t5*( q / (gamma*m0**2 * self.fMass0)) * (self.fG/(1 + gamma))*(Px*Bx+Py*By+Ps*Bs)
        
        
        Sxp =      kappa * Ss + t6 * ((Ps * Ex - Px * Es) * Ss - (Px * Ey - Py * Ex) * Sy) + (sp1*By+sp2*Py)*Ss-(sp1*Bs+sp2*Ps)*Sy
        Syp =                   t6 * ((Px * Ey - Py * Ex) * Sx - (Py * Es - Ps * Ey) * Ss) + (sp1*Bs+sp2*Ps)*Sx-(sp1*Bx+sp2*Px)*Ss
        Ssp = (-1)*kappa * Sx + t6 * ((Py * Es - Ps * Ey) * Sy - (Ps * Ex - Px * Es) * Sx) + (sp1*Bx+sp2*Px)*Sy-(sp1*By+sp2*Py)*Sx
        
        DX = [xp, yp, tp, Pxp/P0c, Pyp/P0c, Wp/self.fKinEn0, Sxp, Syp, Ssp, H]
        
        return DX
    
    def track(self, ElementSeq, ntimes):
        brks = 101
        Xtmp = self._fIniState
        for n in range(1,ntimes+1):
            for i in range(len(ElementSeq)):
                element = ElementSeq[i]
                at = np.linspace(0, element.fLength, brks)
                
                Xtmp = element.frontKick(Xtmp)
                Xtmp = odeint(self._RHS, Xtmp, at, args=(element,))[brks-1]
                Xtmp = element.rearKick(Xtmp)
            self.fState.update({n:Xtmp})
        
        


