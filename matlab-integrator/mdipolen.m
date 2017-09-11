function [f, L, left, right] = mdipole(length, B0, R)

    global m q G c MeV W0
    
    if(B0==j)
        v = W2vel(W0*MeV);
        g = W2gamma(W0*MeV);

        B0 = (m * g * v) / (R * q);
    end
    
    f = @(s, X)emfunction(s, X, 1/R, @(t,x,y,s)Ex(t,x,y,s), @(t,x,y,s)Ey(t,x,y,s), @(t,x,y,s)Es(t,x,y,s), @(t,x,y,s)Bx(t,x,y,s), @(t,x,y,s)By(t,x,y,s, B0, R, v), @(t,x,y,s)Bs(t,x,y,s));
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
function By = By(t,x, y, s, B0, R, v0)
    global m q c MeV W0
    
    kw = 1.000;
    
    v = W2vel(kw*W0*MeV);
    b = v/c;
    g = 1/sqrt(1-(v*v)/(c*c));
    
    k = (q/m)*((2-3*b^2-0.75*b^4)/(v*v) + 1/(g*c*c));

    E0 = -B0*v0;
    B = -E0/v;

    k = k * 1.18;
    
    By = B*(1 + (-1/R + k*B*v)*x +(1/R)*(1/R-3*k*B*v)*x*x);
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