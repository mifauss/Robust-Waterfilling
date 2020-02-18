function [KL_div, sX] = get_KL_div(s0, sN, a)

sX = zeros(size(s0));
if a==0
    sX = s0;
    KL_div = 0;
else
    for i=1:length(s0)
        % find positive real roots of polynomial p
        p = [1, -(s0(i)-2*sN(i)), -sN(i)*(2*s0(i)+a*s0(i)*sN(i)-sN(i)), -s0(i)*sN(i)^2];
        r = roots(p);
        idx = and(imag(r) == 0, real(r) >= 0);
        sX(i) = r(idx);
    end
    snr = sX./s0;
    KL_div = 0.5*(sum(snr-1-log(snr)));
end