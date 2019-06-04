function R = sum_mi_iidGG(betav, h, pv)
% sum of mutual informations 

[pn,k] = size(pv);

R = zeros(1,pn);
for n = 1:pn
    n
    for i = 1:k
        R(1,n) = R(1,n) + miGG(betav(n),h(i),pv(n,i))/k;
    end
end


end