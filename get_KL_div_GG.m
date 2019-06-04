%% KL divergence - GG
function KL_div = get_KL_div_GG(a1, p1, a2, p2)

KL_div = log(p1*a2*gamma(1/p2)) - log(p2*a1*gamma(1/p1)) + (a1/a2)^p2 * gamma((p2+1)/p1) / gamma(1/p1) - 1/p1;

end