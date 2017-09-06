function X = demo(X0, n)

global W0 MeV c m q

R = 24.6642277;
h=0.05;


R1 = @()ecyldeflEy(9.68561959, j, R, h, 10);

lattice = [{R1}];

disp('1) cw (x0 --> x1)...')
X = turn(lattice, X0, n);
x1 = X(:, end);
sy1 = X(8,:);

disp('2) ccw (x0 --> x2)...')
X = turnCCW(lattice, X0, n);
x2 = X(:, end);
sy2 = X(8,:);

disp('3) ccw as inverse for 1): x1 --> x0...')
X = turnCCW(lattice, x1, n);
x0i = X(:, end);
sy0i = X(8,:);



%plot([1:1:n], sy1, 'b.', [1:1:n], sy2, 'r.', [1:1:n], sy0i, 'g.', 'MarkerSize', 2)
plot([1:1:n+1], sy1+sy2, 'b-', [1:1:n+1], sy1-sy0i(end:-1:1), 'r-', 'MarkerSize', 2)


X=[];
end