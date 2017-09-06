% DX = emfunction(X, kappa)
% calculates right-hand part of equation of spin-orbit motion
% x = [x, y, t, px, py, W, Sx, Sy, Ss] 
function DX = emfunction(s, X, kappa, fEx, fEy, fEs, fBx, fBy, fBs)

global m q G c MeV W0
v0 = W2vel(W0*MeV);
g0 = W2gamma(W0*MeV);

W =W0*MeV*(1+X(6));
g = W2gamma(W);
v = W2vel(W);
p0 = (m*v0*g0);
p = (m*v*g);


x =X(1);
y =X(2);
t =X(3);
px=X(4)*p0;
py=X(5)*p0;
Sx=X(7);
Sy=X(8);
Ss=X(9);

hs = 1+kappa*x;

Ex = fEx(t, x,y,s);
Ey = fEy(t, x,y,s);
Es = fEs(t, x,y,s);
Bx = fBx(t, x,y,s);
By = fBy(t, x,y,s);
Bs = fBs(t, x,y,s);

ps = sqrt(p^2 - px^2 - py^2);
xp = (px*hs)/ps;
yp = (py*hs)/ps;
Wp = q*(Ex*xp +Ey*yp +Es);

H = (p*hs)/ps;
vp = W2velprime(W, Wp);
tp = H/v;

gp = Wp2gammaprime(Wp);

D = (q/(m*hs))*(xp*By-yp*Bx+H*Es/v)-((g*v)/(H*hs))*3*kappa*xp;

xpp=((-H*D)/(g*v))*xp+((q*H)/(p))*(H*Ex/v+yp*Bs-hs*By)+kappa*hs;
ypp=((-H*D)/(g*v))*yp+((q*H)/(p))*(H*Ey/v+hs*Bx-xp*Bs);


pxp = px*(vp/v - gp/g)+p*xpp/H-px*((px*xpp)/(p*H)+(py*ypp)/(p*H)+(hs*kappa*xp)/(H^2));
pyp = py*(vp/v - gp/g)+p*ypp/H-py*((px*xpp)/(p*H)+(py*ypp)/(p*H)+(hs*kappa*xp)/(H^2));

t5 = H/v;
t6 =  t5* (q / (g * m * m * c * c)) * (G + 1/(1 + g));
sp1 = t5*(-q / (g*m))*(1 + G * g);
sp2 = t5*( q / (g*m * m * m*c*c)) * (G/(1 + g))*(px*Bx+py*By+ps*Bs);


Sxp =      kappa * Ss + t6 * ((ps * Ex - px * Es) * Ss - (px * Ey - py * Ex) * Sy) + (sp1*By+sp2*py)*Ss-(sp1*Bs+sp2*ps)*Sy;
Syp =                   t6 * ((px * Ey - py * Ex) * Sx - (py * Es - ps * Ey) * Ss) + (sp1*Bs+sp2*ps)*Sx-(sp1*Bx+sp2*px)*Ss;
Ssp = (-1)*kappa * Sx + t6 * ((py * Es - ps * Ey) * Sy - (ps * Ex - px * Es) * Sx) + (sp1*Bx+sp2*px)*Sy-(sp1*By+sp2*py)*Sx;

DX = [ xp; yp; tp; pxp/p0; pyp/p0; Wp/(W0*MeV); Sxp; Syp; Ssp; H];

end