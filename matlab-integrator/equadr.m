function [f, L, left, right] = equadr(length, V, h)

    global m q G c MeV W0
    
    g = 2*V/(h*h);
    
    f = @(s, X)emfunction(s, X, 0, @(t,x,y,s)Ex(t,x,y,s,g), @(t,x,y,s)Ey(t,x,y,s,g), @(t,x,y,s)Es(t,x,y,s), @(t,x,y,s)Bx(t,x,y,s), @(t,x,y,s)By(t,x,y,s), @(t,x,y,s)Bs(t,x,y,s));
    L = length;
    
    left = @(X)leftKick(X, V, h, q, W0*MeV);
    right = @(X)rightKick(X, V, h, q, W0*MeV);
    
end

function Ex = Ex(t,x, y, s, g)
    Ex = -g*x;
end
function Ey = Ey(t,x, y, s, g)
    Ey = g*x;
end
function Es = Es(t,x, y, s)
    Es = 0;
end
function Bx = Bx(t,x, y, s)
    Bx = 0;
end
function By = By(t,x, y, s)
    By = 0;
end
function Bs = Bs(t,x, y, s)
    Bs = 0;
end

function Xk = leftKick(X, V, h, q, W0)
    x=X(1);
    y=X(2);
    u = V*(x*x-y*y)/(h*h);
    Xk = X;
    Xk(6) = Xk(6) - q *u/W0;
end
function Xk = rightKick(X, V, h, q, W0)
    x=X(1);
    y=X(2);
    u = V*(x*x-y*y)/(h*h);
    Xk = X;
    Xk(6) = Xk(6) + q *u/W0;
end