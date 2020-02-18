% i.i.d. GG input, beta_i = beta varying
% fixed non-identical channel gains

clc; clear all;

k = 4;      % number of channels, i = 1,...,k
P = 10;	% power constraint: $\sum_i \var(X_i) \le P$
s0 = 1;     % $E[S_i^2]=1$, where $S_i$ is normalized $X_i$
nv = -1:0.2:2;       % vector of $n_i$'s, where $n_i = \log_2 \beta_i$
betav = 2.^nv;          % vector of $\beta$, in iid case, for all $i$, $\beta_i = \beta$
[px, pn] = size(betav); % pn is the length of betav
h = [8 1 2 7];  % vector of channel gain $h_i$
g = P*h.^2;     % vector of channel strength $\gamma_i$

%% Coefficients of MMSE bounds for Generalized Gaussian distribution
alphaf = @(beta) sqrt(gamma(1/beta)/gamma(3/beta)); % s.t. $E[S_i^2]=1$
wv = ones(pn,k); % coef vector of LMMSE bound

sN = 1;
snr = s0/sN;

for cnt = 1:pn
    
    beta = betav(cnt);
    alpha = alphaf(beta);

    % Reference distribution G(alpha0,2)
    beta0 = 2;
    alpha0 = alphaf(beta0);

    % KL ball radius: epsilon
    eps = get_KL_div_GG(alpha, beta, alpha0, beta0);
    
    % Fisher information
    J = ((beta^2) * gamma(3/beta) * gamma(2-1/beta)) / (s0 * gamma(1/beta)^2);
    
    % mmse lower bound
    lv(cnt,:) = ones(1,k);
    cl = -real(lambertw(0,-exp(-(1+2*eps))));
    clv(cnt,:) = cl*ones(1,k);
    lb(cnt,1) = cl/(snr + 1); % bound 
    
    % mmse upper bound
    u = -real(lambertw(-1,-exp(-(1+2*eps)))*s0);
    uv(cnt,:) = u*ones(1,k);
    ub(cnt,1) = u/(u*snr + 1); % bound
    
    % Cramer Rao bound
    if beta > 0.5
        c = 1/J;%1/(snr+J); 
    else
        c = 0;
    end
    cv(cnt,:) = c*ones(1,k);
    cb(cnt,1) = c/(c*snr + 1); % bound
  
    % LMMSE bound
    wb(cnt,1) = 1/(1+snr);

    % Mismatched waterfilling: using LMMSE bound, coe=1
    [etaw(cnt),pw(cnt,:)] = waterfilling(g,wv(cnt,:),ones(1,k));
    
    % RPA with Cramer-Rao bound
    [etac(cnt),pc(cnt,:)] = waterfilling(g,cv(cnt,:),ones(1,k));

    % RPA
    [etal(cnt),pl(cnt,:)] = waterfilling(g,lv(cnt,:),clv(cnt,:));
    [etalm(cnt),plm(cnt,:)] = waterfilling(g,lv(cnt,:),ones(1,k));
    [etau(cnt),pu(cnt,:)] = waterfilling(g,uv(cnt,:),ones(1,k));
    
end


%% plot quotient of water-levels and c_W v.s. beta
cWwl = max(wv./clv,[],2); % upper bound of the difference between water levels
cWul = max(uv./clv,[],2);
cWwc = max(wv./cv,[],2);

figure
hold on
plot(betav,cWwl,'-k','LineWidth',2)
plot(betav,cWul,'--k','LineWidth',2)
plot(betav,cWwc,'-.k','LineWidth',2)
plot(betav,etaw./etal,'-k')
plot(betav,etau./etal,'--k')
plot(betav,etaw./etac,'-.k')
xlabel('beta');
ylabel('Quotient between water levels');
legend('c_W for eta_w/eta_l','eta_w/eta_l','c_W for eta_u/eta_l','eta_u/eta_l','c_W for eta_w/eta_c','eta_w/eta_c');
legend('location','south');


%% plot ralative rate w.r.t. Rw, the achievable rate of legacy waterfilling
Rw = sum_mi_iidGG(betav,h,P*pw); % computing the achievable rate via Monte Carlo simulations
Rc = sum_mi_iidGG(betav,h,P*pc);
Rl = sum_mi_iidGG(betav,h,P*pl);
Ru = sum_mi_iidGG(betav,h,P*pu);


figure
hold on
plot(betav,Rw./Rw,':k')
plot(betav,Rc./Rw,'--k')
plot(betav,Rl./Rw,'-k')
plot(betav,Ru./Rw,'-.k')
xlabel('beta');
ylabel('Relative achievable rate r(p)');
legend('r(p_w)','r(p_c)','r(p_l)','r(p_u)');
legend('location','south');

%% plot c_R and difference between R(p) v.s. beta
crul = log2(max(uv./clv,[],2))/2;    % using the upper and lower KL-ball bounds
crwl = log2(max(wv./clv,[],2))/2;    % using LMMSE bound and lower KL-ball bound
crwc = log2(max(wv./cv,[],2))/2;    % using LMMSE bound and Cramer-Rao bound

% an upper bound of the optimal rate
Rub = apprx_mi(wv, ones(pn,1)*g, pw);    

figure
hold on
plot(betav,crul','-.k','LineWidth',1.2)
plot(betav,crwl','-k','LineWidth',1.2)
plot(betav,crwc','--k','LineWidth',1.2)
plot(betav,Ru-Rl,'-.k')
plot(betav,Rw-Rl,'-.k')
plot(betav,Rw-Rc,'--k')
xlabel('beta');
ylabel('Difference between achievable rates');
legend('lu','wl','wc','|R_{p_u}-R(p_l)|','|R_{p_w}-R(p_l)|','|R_{p_w}-R(p_c)|');
legend('location','south');
