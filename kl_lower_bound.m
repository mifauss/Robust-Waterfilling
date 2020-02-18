function [mmse_lb, gamma_opt] = kl_lower_bound(S0, SN, kl_eps)
%
% INPUT
% S0:       Reference covariance matrix, needs to be a matrix of size K x
%           K, or a vector of length K
% SN:       Noise covariance matrices, needs to be a matrix of size K x K
%           or a vector of length K
% kl_eps:   Radius of feasible KL divergence ball
%
% OUTPUT
% mmse_lb:      Lower bound on the MMSE
% gamma_opt:    Optimal Lagrange multiplier

% catch vector case
if any(size(S0) == 1)
    S0 = diag(S0);
end
if any(size(SN) == 1)
    SN = diag(SN);
end

% check positive definiteness
[~, p0] = chol(S0);
[~, pN] = chol(SN);

% calculate eigenvalues
if all(size(S0) == size(SN)) && p0 == 0 && pN == 0
    Xi = S0*((S0+SN)\SN);
    xi = eig(Xi);
else
    error('S0 and SN must be positive definite matrices of the same size or positive vectors of the same length.');
end

% solve for gamma
func = @(gamma) sum(log(1+gamma*xi) - gamma.*(xi./(gamma.*xi+1))) - 2*kl_eps;
gamma_max = 1;
while func(gamma_max) <= 0
    gamma_max = 2*gamma_max;
end
gamma_opt = fzero(func, [0 gamma_max]);

% calculate bound
mmse_lb = sum(xi./(gamma_opt*xi+1));