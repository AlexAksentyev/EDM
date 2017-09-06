function [TIME, y] = srk4(f, T, h,  X)
t0 = T(1);
t1 = T(2);
t=t0;
x=X;
y=[x'];
TIME = [t0];
N=100;
eps = 1e-6;
k1 = zeros(length(X), 1);
k2 = zeros(length(X), 1);

while (t <= t1)
    [k1, k2] = ImplicitSolution(f, t, x, k1, k2, h, N, eps);
    x = x + 0.5 * h * (k1 + k2);
    t = t + h;
    y = [y; x'];
    TIME = [TIME; t];
end


end

function [k1, k2] = ImplicitSolution(f, t, x, k1_old, k2_old, h, N, eps)
    coef_GL_4 = 0.5 / sqrt(3);
    k1=k1_old;
    k2 = k2_old;
    count = 0;
    while (count==0 || (count <= N && (norm(k1 - k1_old) > eps || norm(k2 - k2_old) > eps)))
        
        k1_old = k1;
        k2_old = k2;
        
        curs = t + (0.5 + coef_GL_4) * h;
        tmp1 = x + h * (0.25 * k1_old + (0.25 - coef_GL_4) * k2_old);
        tmp2 = x + h * ((0.25 + coef_GL_4) * k1_old + 0.25 * k2_old);
        k1 = f(curs, tmp1);
        k2 = f(curs, tmp2);
        count=count+1;
    end

end