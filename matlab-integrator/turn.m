function X = turn(els, X0, n)
global options
options = odeset('AbsTol',1e-5, 'RelTol',1e-5, 'MaxStep', 1, 'InitialStep', 1);
X = zeros(length(X0), n+1);
X(:,1)=X0;

Xtmp = X0;

perc = fix(n/10);
N=1;
for k=1:n
    for i=1:length(els)
        elf = els{i}; 
        [element, L, leftcorrection, rightcorrection] = elf(); 
        Xtmp = leftcorrection(Xtmp);
        [t, x] = ode45(element, [0 L], Xtmp, options);
        %x
        Xtmp = x(end, :);
        Xtmp=Xtmp';
        Xtmp = rightcorrection(Xtmp);
        %x=x';
        %X = [X Xtmp];
    end
    X(:, k+1)=Xtmp;
     if(k == perc*N)
         fprintf('%i0/100\n', N)
         N=N+1;
     end
end