% fixed non-identical GG input, [beta_i] = [1,2,4,infty]
% varying identical channel gains, h_i=h varying from 0.5 to 2

clc; clear all;

k = 4;      % number of channels, i = 1,...,k
s0 = 1; % E[si^2]=1
P = 10; % power constraint
betav = [1, 2, 4];  % vector of $\beta_i$, i=1,...,k-1, and $\beta_k = \infty$
hv = 10:10:100;    % channel gain scaler
[px,pn] = size(hv);

%% Coefficients of MMSE bounds 
wv = ones(pn,k);% LMMSE bound 1/(1+snr)

% Generalized Gaussian distribution
alphaf = @(beta) sqrt(gamma(1/beta)/gamma(3/beta)); % s.t. $E[S_i^2]=1$

% for i=1,...,k-1, Si ~ G(alpha,beta) s.t. E[Si^2]=1
for bn = 1:(k-1)
    
    beta = betav(bn);
    alpha = alphaf(beta);
    
    % Reference distribution G(alpha0,2)
    beta0 = 2; 
    alpha0 = alphaf(beta0);
    
    % KL ball radius: epsilon
    epsG = get_KL_div_GG(alpha, beta, alpha0, beta0);
        
    % mmse lower bound
    lv(:,bn) = s0*ones(pn,1);   % coe
    clv(:,bn) = -real(lambertw(0,-exp(-(1+2*epsG))))*ones(pn,1);
    
    % mmse upper bound
    sX_max = -real(lambertw(-1,-exp(-(1+2*epsG)))*s0);
    %     ub = sX_max*sN/(sX_max + sN);
    uv(:,bn) = sX_max*ones(pn,1);

end

% Uniform: Sk ~ U_K(R)
K = 1;	% dimension
R = 1;  % radius, s.t. E[Si^2]=1

% Reference distribution N(0,s0),s0 = 1
epsU = get_KL_div_U(K); % KL ball radius: epsilon

% mmse lower bound
lv(:,k) = s0*ones(pn,1);   % coe
cl = -real(lambertw(0,-exp(-(1+2*epsU))));
clv(:,k) = cl*ones(pn,1);

% mmse upper bound
sX_max = -real(lambertw(-1,-exp(-(1+2*epsU/K)))*s0);
uv(:,k) = sX_max*ones(pn,1);


%% RPA
for cnt = 1:pn
    h = hv(cnt).*ones(1,k); % identical channel gains, $h_i = h$
    g = P*h.^2;             % $\gamma_i$: strength of the channels
    gv(cnt,:) = g;
    
    % Mismatched waterfilling: using LMMSE bound, coe=1
    [etaw(cnt),pw(cnt,:)] = waterfilling(g,wv(cnt,:),ones(1,k));
    
    % RPA
    [etal(cnt),pl(cnt,:)] = waterfilling(g,lv(cnt,:),clv(cnt,:));
    [etalm(cnt),plm(cnt,:)] = waterfilling(g,lv(cnt,:),ones(1,k));
    [etau(cnt),pu(cnt,:)] = waterfilling(g,uv(cnt,:),ones(1,k));
    
end

    
%% plot quotient of water-levels and c_W v.s. beta
cWwl = max(wv./clv,[],2); 
cWul = max(uv./clv,[],2);

figure
hold on
plot(hv,cWwl,'-k','LineWidth',1.2)
plot(hv,cWul,'--k','LineWidth',1.2)
plot(hv,etaw./etal,'-k')
plot(hv,etau./etal,'--k')
xlabel('beta');
ylabel('Quotient between water levels');
legend('c_W for eta_w/eta_l','c_W for eta_u/eta_l','eta_w/eta_l','eta_u/eta_l');
legend('location','south');


%% plot ralative rate w.r.t. Rw, the achievable rate of legacy waterfilling
Rw = sum_mi_niidGG(betav,hv'*ones(1,k),P*pw);
Rl = sum_mi_niidGG(betav,hv'*ones(1,k),P*pl);
Ru = sum_mi_niidGG(betav,hv'*ones(1,k),P*pu);

figure
hold on
plot(hv,Rw./Rw,'--k')
plot(hv,Rl./Rw,'-k')
plot(hv,Ru./Rw,'-.k')
xlabel('Channel gains h_i');
ylabel('Relative achievable rate r(p)');
legend('r(p_w)','r(p_l)','r(p_u)');
legend('location','south');


%% plot c_R and difference between R(p) v.s. h_i
crul = log(cWul)/2;    % using the upper and lower KL-ball bounds
crwl = log(cWwl)/2;    % using LMMSE bound and lower KL-ball bound

figure
hold on
plot(hv,crul','-k','LineWidth',1.2)
plot(hv,crwl','-.k','LineWidth',1.2)
plot(hv,abs(Ru-Rl),'-k')
plot(hv,abs(Rw-Rl),'-.k')
xlabel('beta');
ylabel('Difference between achievable rates');
legend('c_R for R_{p_u}-R(p_l)','c_R for R_{p_w}-R(p_l)','|R_{u}-R(p_l)|','|R_{w}-R(p_l)|');
