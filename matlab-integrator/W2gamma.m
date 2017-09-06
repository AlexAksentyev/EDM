% g = W2gamma(W)
% converts Kinetic Energy W to Lorentz factor g
function g = W2gamma(W)
    global m c
    g = W/(m*c*c)+1;
end