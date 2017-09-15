clc; clear

disp("Running DP\n");

initBDT; setD;

%global MeV AMU AMUMEV DEGRAD c m q W0

v0 = W2vel(W0*MeV);
g0 = W2gamma(W0*MeV);

p0 = (m/AMU*AMUMEV*v0/c*g0); %in MeV

L=1.8; B0=.46;
R=p0/B0/c*1e6;
PHI=L/R/DEGRAD;

h=0.05;

R1 = @()mdipole(L, j, R);

lattice = [{R1}];

n=100
X0 = [1e-3*ones(1,3) -2e-3:2e-3:2e-3; [-2e-3:2e-3:2e-3] 1e-3*ones(1,3); zeros(6,6); ones(1,6); zeros(1,6)];
X = cell(6)
for i=1:4
  disp(num2str(i));
  X{i} = turn(lattice, X0(:,i), n);
end

x=[X{1}(1,:); X{2}(1,:); X{3}(1,:); X{4}(1,:)]
y=[X{1}(2,:); X{2}(2,:); X{3}(2,:); X{4}(2,:)]
sx=[X{1}(7,:); X{2}(7,:); X{3}(7,:); X{4}(7,:)]
sy=[X{1}(8,:); X{2}(8,:); X{3}(8,:); X{4}(8,:)]



for i=1:4; subplot(2,2,1); plot(x(i,:)); hold all; title('x'); end; legend("p1","p2","p3","p4");
for i=1:4; subplot(2,2,2); plot(y(i,:)); hold all; title('y'); end; legend("p1","p2","p3","p4");
for i=1:4; subplot(2,2,3); plot(sx(i,:)); hold all; title('Sx'); end; legend("p1","p2","p3","p4");
for i=1:4; subplot(2,2,4); plot(sy(i,:)); hold all; title('Sy'); end; legend("p1","p2","p3","p4");

