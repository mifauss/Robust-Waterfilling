function [Ec] = apprx_mi(c,g,p)
% approximated mutual information
% all inputs are pn x k martix, pn: # points, k: # channel

[pn,k] = size(p);

lmmsei = @(ci,x) ci./(1+ci.*x);
Iic = @(ci,snr) 1/2*integral(@(x)lmmsei(ci,x),0,ci*snr);

Ec = zeros(1,pn);
for n = 1:pn
    for i = 1:k
        Ec(n) = Ec(n) + Iic(c(n,i),g(n,i)*p(n,i))/k;
    end
end


end