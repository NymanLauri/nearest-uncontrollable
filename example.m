options.maxiter = 50;
options.maxtime = 100;
options.verbosity = 2;

n = 10;

rng(1,"twister")

% %Real input
% A = randn(n);
% b = randn(n,1);

%Complex input
A = randn(n) + 1i*randn(n) ;
b = randn(n,1) + 1i*randn(n,1);

[S,t,distance,time_seconds,Q,infotable] = nearest_uncontrollable(A, b, options);

% The constructed uncontrollable pencil is [t S+xI]

% Construct the canonical form
S_tri = Q'*S*Q;
t_tri = Q'*t;

% % The constant part of the canonical form
% [t_tri S_tri]

