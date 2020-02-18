function [ mmse_ub, sX_opt, alpha_opt ] = kl_upper_bound( s0, sN, kl_eps )
%
% INPUT
% s0:       Eigenvalues of reference covariance matrix, vector of length K
% sN:       Eigenvalues of noise covariance matrix, vecor of length K
% kl_eps:   Radius of feasible KL divergence ball
%
% OUTPUT
% mmse_ub:      Upper bound on the MMSE
% sX_opt:       Eigenvalues of the covariance matrix that attains the bound
% alpha_opt:    Optimal Lagrange multiplier

% initialize
s0 = s0(:);
sN = sN(:);

if ~(length(s0) == length(sN) && all(s0 > 0) && all(sN > 0))
    error('s0 and sN must be positive vectors of the same length.');
end

% catch KL div = 0 case
if eps == 0
    alpha_opt = 0;
else
    % objective function for root finding
    func = @(a) get_KL_div(s0, sN, a) - kl_eps;
    
    % calculate bounds on a
    sX_bound = 2*(log(2)+2*kl_eps)*s0;
    a_bound = (sX_bound./sN+1).^2.*(1./s0 - 1./sX_bound)+1;
    
    % find root
    alpha_opt = fzero(func, [0, min(a_bound)]);
end

% get optimal variances
[~, sX_opt] = get_KL_div(s0, sN, alpha_opt);

% calculate bound
mmse_ub = sum(sX_opt.*sN./(sX_opt+sN));

end

