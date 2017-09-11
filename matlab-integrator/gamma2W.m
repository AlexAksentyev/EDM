% W = gamma2vel(g)
% converts gamma factor g to energy W in MeV
function W = gamma2W(g)
    global c m MeV W0
    
    v = gamma2vel(g);
    W = m*c*c*(g-1)/MeV;
end
