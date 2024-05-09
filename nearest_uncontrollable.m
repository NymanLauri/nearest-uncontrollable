function [S,t,distance,time_seconds,Q,infotable] = nearest_uncontrollable(A, b, options, x0)
% Computes a (locally) nearest complex uncontrollable pencil [t S+xI]
% to the rectangular pencil [b A+xI]. If the pencil [b A+xI] and the 
% starting point x0 are both real, the output [t S+xI] will be real.
% 
% Input:
%   b, A (vector of size n, matrix of size (n,n))
%       the constituents of the pencil [b A+xI] in consideration
%   maxiter (integer)
%        maximum amount of outer iterations, used in Manopt's
%        "trustregions" solver
%   timemax_seconds (integer)
%        maximum amount of time (in seconds) the algorithm can run before
%        the iteration is stopped
%   x0 (matrix of size (n,n))
%       initial value of the optimization variable (unitary matrix)
%       
%
% Output:
%   t, S (vector of size n, matrix of size (n,n))
%       the constituents of the nearest uncontrollable pencil [t S+xI]
%   distance (vector)
%       the distance to the constructed uncontrollable pencil after each 
%       iteration
%   time_seconds (vector)
%       elapsed time after each iteration 
%   Q (matrix of size (n,n))
%       final value of the optimization variable (unitary matrix)
%   infotable (table)
%       contains various types of information for diagnostic purposes
%
% Requirement: Manopt needs to be imported

n = length(A);

if not(exist('options', 'var'))
    options = struct();
end
if not(exist('x0', 'var'))
    x0 = [];
end
if isreal([A,b]) && isreal(x0)
    problem.M = stiefelfactory(n, n, 1); 
else
    problem.M = stiefelcomplexfactory(n, n, 1); 
end

% % Rescale the pencil to be of norm 100
% P_norm = norm([A B], 'f')*1e-2;
% A = A / P_norm;
% B = B / P_norm;

problem.cost = @cost;
% The code uses the Euclidean gradient. Projection to 
% the tangent space of U(n) is handled automatically (see 
% stiefelcomplexfactory documentation)
problem.egrad = @egrad;
% Euclidean Hessian. Projection is handled automatically.
problem.ehess = @ehess;

default.tolgradnorm = 1e-10;
default.maxiter = 100;
default.maxtime = 100;
default.verbosity = 1; % 2 = Default; 0 = No output; 

options = mergeOptions(default, options);

[Q, xcost, info, ~] = trustregions(problem, x0, options);

infotable = struct2table(info);
distance = sqrt(infotable.cost);
time_seconds = infotable.time;

% Construct the nearest singular pencil
t = triu(Q'*b);
S = triu(Q'*A*Q,-1);

[~,k] = min(diag(abs([t S])));

if k==1
    t(1) = 0;
else
    S(k,k-1) = 0;
end

t = Q*t;
S = Q*S*Q';

% Squared distance to the uncontrollable pencil. Should be equal to xcost
assert(abs(xcost -  norm([b A] - [t S],'f')^2) < 1e-10)

% % Rescale back
% T = T*P_norm;
% S = S*P_norm;
% distance = distance*P_norm;

% ---------

function f = cost(Q)

    t = Q'*b;
    S = Q'*A*Q;
    f = norm(tril([t S],-1),'fro')^2 + min(diag(abs([t S])).^2);
    
end

function g = egrad(Q)

    tS = Q'*[b A]*blkdiag(1,Q);

    L = tril(tS,-1);

    % Add "min" part of the objective function
    [~,k] = min(abs(diag(tS)).^2);
    L(k,k) = tS(k,k);

    g = 2* [b A]*blkdiag(1,Q)*L' + 2* (Q'*A)' * L(:,2:end);
    
end

function H = ehess(Q, d)

    tS = Q'*[b A]*blkdiag(1,Q);

    L = tril(ones(size(tS)),-1);

    % Add "min" part of the objective function
    [~,k] = min(abs(diag(tS)).^2);
    L(k,k) = 1;

    DL = L.*(d'*[b A]*blkdiag(1,Q) + Q'*[b A]*blkdiag(0,d));

    H = [b A]*blkdiag(0,d)*(L.*tS)' + [b A]*blkdiag(1,Q)*DL' ...
        + (d'*A)' * (L(:,2:end).*tS(:,2:end)) + (Q'*A)' * DL(:,2:end);
    H = 2*H;
   
end

end