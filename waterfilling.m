function [ eta,pv ] = waterfilling( g, a, b )
%   waterfilling
%   g: 1xk, channel gains gamma_i
%   a and bc: coefficient of MMSE bound $b \lmmse_a(snr)$


[pn,k] = size(g);

v = g.*a.*b;
[Vs,I] = sort(v,2,'descend');

for j = 1:k, Bs(:,j) = b(:,I(:,j)); end
for j = 1:k, Hs(:,j) = g(:,I(:,j)).*a(:,I(:,j)); end

Hext = [Hs -1];

eta = inf;
i = 1;
while i<=k
    etai = (sum(Bs(1:i)))/(1+sum(1./Hs(1:i)));
    if etai > Hext(i+1) && etai <= Hs(i)
        eta = etai;
        i = k+1;
    else
        i = i+1;
    end
end

pv = max(0, b./eta - 1./(g.*a));

end


