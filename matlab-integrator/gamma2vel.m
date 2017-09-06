% v = gamma2vel(g)
% converts gamma factor g to velocity v
function v = gamma2vel(g)
    global c
    v = c*sqrt(1-1/(g*g));
end
