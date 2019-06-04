%% KL divergence - Uniform
function KL_div = get_KL_div_U(K)

KL_div = K/2 + K/2*log(1/(1+K/2)) + log(gamma(K/2+1));

end