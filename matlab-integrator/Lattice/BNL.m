clc; clear;

initBDT; setD;

G = -.14;
h = .05;
R = 9.207;
Ebp = -12e6; %V/m
Lbp = pi/16 *(m*c*c/MeV)*g^3*b^2*G*1e6/(Ebp*(G+1));

Lq = 5;
GSFP = 0;
GSDP = 0;

drift_25 = @()emdrift(.25);
drift_15 = @()emdrift(.15);
drift2_2 = @()emdrift(2.2);
bpm = @(15);

dquadSS1 = @()mquadr(Lq, -.86);
fquadSS1 = @()mquadr(Lq, .831);

fquadSS2 = @()mquadr(Lq, .831);
dquadSS2 = @()mquadr(Lq, -.86);

fquadA = @()mquadr(Lq,  1.364);
dquadrA = @()mquadr(Lq,-1.023);

fsexA = @()msext(.15, GSFP);
dsexA = @()msext(.15, GSDP);

wienA = @()wien(Lbp, j, R, h);




