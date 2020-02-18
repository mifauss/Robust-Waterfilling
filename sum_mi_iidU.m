function R = sum_mi_iidU(h, p)
% mutual information of the parallel channels
% h: pn x k, channel gains h_i
% p: pn x k, power p_i * P

[pn,k] = size(p);

R = zeros(1,pn);
for n = 1:pn
    n
    for i = 1:k
    R(n) = R(n) + miUnif(h(n,i),p(n,i))/k;
    end
end


end