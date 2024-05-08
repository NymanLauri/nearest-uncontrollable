% addpath("./Overton")
% addpath("./Overton/bfgs0_91")

addpath("./software_rob_sc/fo_dist_uncont")

n = 8;

rng(1,"twister")

%Real input
A = randn(n);
b = randn(n,1);

% %Complex input
% A = randn(n) + 1i*randn(n) ;
% b = randn(n,1) + 1i*randn(n,1);

[f,z,tol] = dist_uncont_hybrid(A, b, options)