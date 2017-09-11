% to start work use followig commands:
% >> initBDT;   %initialization of necessary global variables
% >> setProton; %for proton
%
% available elements:
%
% · drift(length)
% //---------------------------------------------------------------------
% MAGNETIC FIELDS:
% · mdipole(length, B0, R)
%     represents a horizontally bending dipole B_y=B0
%     if B0 is equal to j then B0 = p/(R*q)
%
% · mquadr(length, g)
%     represents a normal quadrupole B_x=-g*y, B_y=-g*x
%
% · msquadr(length, g)
%     represents a skew quadrupole B_x=g*x, B_y=-g*y
%
% · msext(length, g)
%     represents a normal sextupole B_x=gp*x*y, B_y=0.5*gp*(x^2-y^2)
%
% · msol(length, g)
%     represents a normal solenoid B_z=Bz
%
% · mextractor(length, B0)
%     represent a horizontally bending dipole B_y=B0 with respect to the
%     straight direction
% //---------------------------------------------------------------------
% ELECTROSTATIC FIELDS:
% //--------------------------------------------------------------------