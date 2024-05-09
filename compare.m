% addpath("./Overton")
% addpath("./Overton/bfgs0_91")

% addpath("./software_rob_sc/fo_dist_uncont")
% addpath("./software_rob_sc/fo_dist_uncont/private_folder")

rng(1,"twister")

options = struct();
options.eig_method = 0;
options.method = 1;
options.print_dtl = 1;
options.tol = 1e-8;
% options.prtlevel = 1;

options_riemann.maxiter = 50;
options_riemann.maxtime = 100;
options_riemann.verbosity = 1;

list_sizes = [15];
n_sample = 10;

d_riemann = zeros(list_sizes(end), n_sample);
d_them = zeros(list_sizes(end), n_sample);

t_riemann = zeros(list_sizes(end), n_sample);
t_them = zeros(list_sizes(end), n_sample);

for j = 1:n_sample
    j
    for m = list_sizes

        % Complex input
        A = randn(m) + 1i*randn(m) ;
        b = randn(m,1) + 1i*randn(m,1);
              
        tic
        [f,z,tol] = dist_uncont_hybrid(A, b, options);
        t1 = toc;
        
        tic
        [S,t,distance,Q] = outer(A, b, options_riemann);
        t2 = toc;

        d_riemann(m,j) = distance(end);
        d_them(m,j) = f;
    
        t_riemann(m,j) = t2;
        t_them(m,j) = t1;
    end
end

[mean(t_riemann(m,:),2) mean(t_them(m,:),2)]
[mean(d_riemann(m,:),2) mean(d_them(m,:),2)]



