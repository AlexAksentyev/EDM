import numpy as NP




class Element:
    
    fCurve = None
    fLength = None
    
    def __init__(self, Curve, Length):
        self.fCurve = Curve
        self.fLength = Length
    
    def EField(self,arg):
        return (0,0,0)
    
    def BField(self,arg):
        return (0,0,0)

    def frontKick(self,particle):
        pass # do nothing
    
    def rearKick(self,particle):
        pass # do nothing

class Drift(Element):
    """ drift space
    """
    
    def __init__(self, Length):
        Element.__init__(self, 0, Length)
        

class MQuad(Element):
    """ magnetic quadrupole
    """
    
    _fGrad = None
    
    def __init__(self, Length, Grad):
        Element.__init__(self, 0, Length)
        self._fGrad = Grad
        
    def BField(self, arg):
        x,y = arg[0:2]
        return (-self._fGrad*x, -self._fGrad*y,0)
        

class MDipole(Element):
    """ bending magnetic dipole (horizontally/vertically bending);
    define _BField as a tuple;
    could also be used as a solenoid, if _BField = (0,0,Bz)
    """
    
    _BField = None
    
    def __init__(self, Length, R, BField):
        Element.__init__(self, 1/R, Length)
        self.setBField(BField)
        
    def setBField(self,BField):
        self._BField = BField[:]
        
    def getBField(self):
        return self._BField[:]
    
    @staticmethod    
    def computeBStrength(particle, R):
        return particle.Pc(particle.fKinEn0)*1e6/R
    
    @staticmethod
    def computeRaduis(particle, BField):
        return particle.Pc(particle.fKinEn0)/(particle.CLIGHT()*1e-6*BField)
        
        
    def BField(self, arg):
        return self._BField
        

class MSext(Element):
    """ magnetic sextupole
    """
    
    _fGrad = None
    
    def __init__(self, Length, Grad):
        Element.__init__(self, 0, Length)
        self._fGrad = Grad
        
    def BField(self, arg):
        x,y=arg[0:2]
        return (self._fGrad*x*y,.5*self._fGrad*(x**2 - y**2), 0)
        
class Wien(Element):
    """ wien filter ?!?! 
    The magnetic field definition is weird
    """
    
    _R = None
    _Volt = None
    _BField = None
    
    def __init__(self, Length, R, Hgap, Voltage, BField):
        Element.__init__(self, 1/R, Length)
        self._R = [R, R - Hgap, R**2/(R-Hgap)]
        self._Volt = Voltage
        self._BField = BField
        
    @staticmethod
    def computeVoltage(particle, R, Hgap):
        gamma,beta = particle.GammaBeta(particle.fKinEn0)
        v = beta*particle.CLIGHT()
        R = [R, R-Hgap, R**2/(R-Hgap)]
        E0 = - particle.Pc(particle.fKinEn0)*v/R[0] * 1e6
        return (E0 * R[0] * NP.log(R[2] / R[1])) / (-2)
    
    @staticmethod
    def computeBStrength(particle, R, Hgap): # no idea about the formulas here
        x = particle[0]
        gamma,beta = particle.GammaBeta(particle.fKinEn0)
        v = beta*particle.CLIGHT()
        R = [R, R-Hgap, R**2/(R-Hgap)]
        E0 = - particle.Pc(particle.fKinEn0)*v/R[0] * 1e6
        qm = particle.EZERO()/particle.fMass0
        k = qm * ((2 - 3*beta**2 - .75*beta**4)/beta**2 + 1/gamma)
        
        B0 = -E0/v
        k = k*1.18
        return 0.018935*(-B0)*(-1/R + k*B0*v)*x
        
    def EField(self, arg):
        x = arg[0]
        Ex = -2*self._Volt/(NP.log(self._R[2]/self._R[1])*(self._R[0]+x))
        return (Ex, 0, 0)
    
    def BField(self, arg):        
        return (0, self._BField, 0)
    
    def frontKick(self, particle):
        x=particle._fState[0]
        R = self._R[0]
        R1 = self._R[1]
        R2 = self._R[2]
        V = self._Volt
        u = -V + 2*V*NP.log((R+x)/R1)/NP.log(R2/R1)
        Xk = particle.getState()
        Xk[5] -= u*1e-6/particle.fKinEn0
        particle.setState(Xk)
        
    def rearKick(self, particle):
        x=particle._fState[0]
        R = self._R[0]
        R1 = self._R[1]
        R2 = self._R[2]
        V = self._Volt
        u = -V + 2*V*NP.log((R+x)/R1)/NP.log(R2/R1)
        Xk = particle.getState()
        Xk[5] += u*1e-6/particle.fKinEn0
        particle.setState(Xk)