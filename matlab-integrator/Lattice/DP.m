clc; clear

initBDT; setD;

%global MeV AMU AMUMEV DEGRAD c m q W0

v0 = W2vel(W0*MeV);
g0 = W2gamma(W0*MeV);

p0 = (m/AMU*AMUMEV*v0/c*g0); %in MeV

L=1.8; B0=.46;
R=p0/B0/c*1e6;
PHI=L/R/DEGRAD;

h=0.05;

R1 = @()ecyldeflm(L, j, R, h);

lattice = [{R1}];

X0 = [2e-3 zeros(1,7) 1 0]; n=100;
X = turn(lattice, X0, n);
x = X(1, :);
y = X(2, :);
sy = X(8,:);
sx = X(7,:);


subplot(2,2,1); plot([1:1:n+1], x); title('x');
subplot(2,2,2); plot([1:1:n+1], y); title('y');
subplot(2,2,3); plot([1:1:n+1], sx); title('Sx');
subplot(2,2,4); plot([1:1:n+1], sy); title('Sy');