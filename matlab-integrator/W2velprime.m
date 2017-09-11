% vp = W2velprime(W, Wp)
% converts Kinetic Energy W and its prime to velocity prime vp
function vp = W2velprime(W, Wp)
    global m c
    vp = (Wp*c*(m*c*c)^2)/((W+m*c*c)^2*sqrt(W*W+2*W*m*c*c));
end
