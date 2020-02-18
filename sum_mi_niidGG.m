function R = sum_mi_niidGG(beta, h, p)
% mutual information of the parallel channels
% beta: 1 x k, parameter of GG dist
% h: pn x k, channel gains h_i
% p: pn x k, power p_i * P

[pn,k] = size(p);
lmmsei = @(h,x) h./(1 + h.*x);

R = zeros(1,pn);
for n = 1:pn
    n
    for i = 1:(k-1)
        if beta(i)==2
        % Gaussian
        R(n) = R(n)+1/2 * integral(@(x)lmmsei(h(n,i),x),0,h(n,i)*p(n,i)) /k;
        else
        % GG
        R(n) = R(n) + miGG(beta(i),h(n,i),p(n,i))/k;
        end
    end
    % Uni
    i = k;
    R(n) = R(n) + miUnif(h(n,i),p(n,i))/k;
end


end