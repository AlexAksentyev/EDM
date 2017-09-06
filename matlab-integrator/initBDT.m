% creation of global variables
global m q G c MeV W0 options
MeV = 1.60217646e-13;
c = 299792458;

options = odeset('AbsTol',1e-5, 'RelTol',1e-5, 'MaxStep', 0.1);