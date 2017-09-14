




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

    def frontKick(self,arg):
        return arg
    
    def rearKick(self,arg):
        return arg

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
        self._BField = BField
        
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
        
