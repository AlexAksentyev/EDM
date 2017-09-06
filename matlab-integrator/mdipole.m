% · mdipole(length, B0, R)
%     represents a horizontally bending dipole B_y=B0
%     if B0 is equal to j then B0 = p/(R*q)
function [f, L, left, right] = mdipole(length, B0, R)

    global m q G c MeV W0
    
    if(B0==j)
        v = W2vel(W0*MeV);
        g = W2gamma(W0*MeV);

        B0 = (m * g * v) / (R * q);
    end
    if(B0==-j)
        v = W2vel(W0*MeV);
        g = W2gamma(W0*MeV);

        B0 = -(m * g * v) / (R * q);
    end
    f = @(s, X)emfunction(s, X, 1/R, @(t,x,y,s)Ex(t,x,y,s), @(t,x,y,s)Ey(t,x,y,s), @(t,x,y,s)Es(t,x,y,s), @(t,x,y,s)Bx(t,x,y,s), @(t,x,y,s)By(t,x,y,s, B0), @(t,x,y,s)Bs(t,x,y,s));
    L = length;
    
    left = @(X)leftKick(X);
    right = @(X)rightKick(X);
    
end

function Ex = Ex(t,x, y, s)
    Ex = 0;
end
function Ey = Ey(t,x, y, s)
    Ey = 0;
end
function Es = Es(t,x, y, s)
    Es = 0;
end
function Bx = Bx(t,x, y, s)
    Bx = 0;
end
function By = By(t,x, y, s, B0)
    By = B0;
end
function Bs = Bs(t,x, y, s)
    Bs = 0;
end

function Xk = leftKick(X)
    Xk = X;
end
function Xk = rightKick(X)
    Xk = X;
end