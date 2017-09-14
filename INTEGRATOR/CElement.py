




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
    
    
    