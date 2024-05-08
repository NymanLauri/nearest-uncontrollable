% addpath("./Overton")
% addpath("./Overton/bfgs0_91")

% addpath("./software_rob_sc/fo_dist_uncont")
% addpath("./software_rob_sc/fo_dist_uncont/private_folder")

n = 10;

% rng(1,"twister")

%Real input
A = randn(n);
b = randn(n,1);

% %Complex input
% A = randn(n) + 1i*randn(n) ;
% b = randn(n,1) + 1i*randn(n,1);

options = struct();
options.eig_method = 0;
options.method = 2;
options.print_dtl = 1;
options.tol = 0.001;

[f,z,tol] = dist_uncont_hybrid(A, b, options)

% svd([A-A^0*z b])

% ----------------
% 
% maxiter = 50;
% maxtime = 100;
% [S,t,distance,time_seconds,Q,infotable] = nearest_uncontrollable(A, b, maxiter, maxtime);
% 
% % The constructed uncontrollable pencil is [t S+xI]
% 
% % Construct the canonical form
% S_tri = Q'*S*Q;
% t_tri = Q'*t;
% 
% f
% distance(end)

