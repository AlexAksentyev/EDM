% creation of global variables
global MeV c DEGRAD AMUMEV AMU options
MeV = 1.60217646e-13;
c = 299792458;
DEGRAD = 180/pi; 
AMUMEV = 931.4940136;
AMU = 1.66053873e-27;

options = odeset('AbsTol',1e-5, 'RelTol',1e-5, 'MaxStep', 0.1);