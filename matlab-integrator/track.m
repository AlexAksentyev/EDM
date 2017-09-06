function X = track(els, X0, n)
global options
%options = odeset('AbsTol',1e-5, 'RelTol',1e-5, 'MaxStep', 0.1);
X=[];
Xtmp = X0;
for k=1:n
    for i=1:length(els)
        elf = els{i};
        [element, L, leftcorrection, rightcorrection] = elf();
        Xtmp = leftcorrection(Xtmp);
        [t, x] = ode45(element, [0 L], Xtmp, options);
        Xtmp = x(length(x), :);
        Xtmp=Xtmp';
        Xtmp = rightcorrection(Xtmp);
        x=x';
        X = [X x];
    end
end