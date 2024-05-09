function [S,t,distance,Q] = outer(A, b, options)

tol = 1e-8;

n = length(A);

sinvals = svd([b A]);
up_bound = min(sinvals);
low_bound = 0;
eigs_tol = 10^-7;
est_normAK = 10^4*max(sinvals);

iter = 0;
eig_method = 0;
print_dtl = 2;

z=0;

while up_bound-low_bound > tol
    iter = iter + 1;
  
    [U,S,V] = svd([b A-A^0*z]);
    S(n,n) = 0;
    bA_trun = U*S*V';
    bb = bA_trun(:,1);
    AA = bA_trun(:,2:end) + A^0*z;

    krylov_matrix = zeros(n,n);
    for i=0:n-1
        krylov_matrix(:,i+1) = AA^i*bb;
    end
    [Q,R] = qr(krylov_matrix);
    
    x0 = Q;
   
    [S,t,distance,time_seconds,Q,infotable] = nearest_uncontrollable(A, b, options, x0);
    
    delta1 = max(distance(end) - 0.45*tol,tol);
    delta2 = max(distance(end)- 0.9*tol,0);

    [test1_test2,zt,iter2] = perform_test_hybrid(A,b,delta1,delta2, ...
                                            eig_method,print_dtl, ...
                                            est_normAK, eigs_tol, iter);
    
    up_bound = distance(end);
    if test1_test2
        up_bound = delta1;
        z = zt;                
    else
        low_bound = delta2;
    end 

    if print_dtl
        fprintf('iter:%d (%3.10f %3.10f] z:%3.10f+%3.10fi\n', iter, low_bound, ...
                                        up_bound, real(z),imag(z));
    end

end



end