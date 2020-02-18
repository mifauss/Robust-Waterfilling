function R = sum_mi_iidGG(betav, h, pv)
% sum of mutual informations 

[pn,k] = size(pv);
lmmsei = @(h,x) h./(1 + h.*x);

R = zeros(1,pn);
for n = 1:pn
    n
    for i = 1:k
        if betav(n)==2
        % Gaussian
        R(1,n) = R(1,n)+1/2 * integral(@(x)lmmsei(h(i),x),0,h(i)*pv(n,i)) /k;
        else
        R(1,n) = R(1,n) + miGG(betav(n),h(i),pv(n,i))/k;
        end
    end
end


end