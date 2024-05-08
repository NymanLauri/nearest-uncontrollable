maxiter = 50;
maxtime = 100;

n = 8;

rng(1,"twister")

% %Real input
% A = randn(n);
% b = randn(n,1);

%Complex input
A = randn(n) + 1i*randn(n) ;
b = randn(n,1) + 1i*randn(n,1);


[S,t,distance,time_seconds,Q,infotable] = nearest_uncontrollable(A, b, maxiter, maxtime);

S_tri = Q'*S*Q;
t_tri = Q'*t;

