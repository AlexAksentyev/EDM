% gp = Wp2gammaprime(W)
% converts Kinetic Energy prime Wp to Lorentz factor prime gp
function gp = Wp2gammaprime(Wp)
    global m c
    gp = Wp/(m*c*c);
end