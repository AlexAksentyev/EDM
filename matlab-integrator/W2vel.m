% v = W2vel(W)
% converts Kinetic Energy W to velocity v
function v = W2vel(W)
    global m c
    v = (c*sqrt(W*W+2*W*m*c*c))/(W+m*c*c);
end
