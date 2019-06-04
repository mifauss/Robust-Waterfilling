function [ eta,pv ] = waterfilling( g )
%   waterfilling
%   g: 1xk, channel gains gamma_i


[pn,k] = size(g);

gs = sort(g,2,'descend');

Hext = [gs 0];

eta = inf;
for i = 1:k
    etai = i/(1+sum(1./gs(1:i)));
    if etai > Hext(i+1) && etai <= gs(i)
        eta = etai;
    end
end

fiv = 1/eta - 1./g;
pv = max(0,fiv);

end


