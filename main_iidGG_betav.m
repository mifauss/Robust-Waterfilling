% i.i.d. GG input, beta_i = beta varying
% fixed non-identical channel gains

clc; clear all;

k = 4;      % number of channels, i = 1,...,k
P = 100;	% power constraint: $\sum_i \var(X_i) \le P$
s0 = 1;     % $E[S_i^2]=1$, where $S_i$ is normalized $X_i$
nv = -1:0.25:2;       % vector of $n_i$'s, where $n_i = \log_2 \beta_i$
betav = 2.^nv;          % vector of $\beta$, in iid case, for all $i$, $\beta_i = \beta$
[px, pn] = size(betav); % pn is the length of betav
h = [8 1 2 7];  % vector of channel gain $h_i$
g = P*h.^2;     % vector of channel strength $\gamma_i$

%% Coefficients of MMSE bounds for Generalized Gaussian distribution
alphaf = @(beta) sqrt(gamma(1/beta)/gamma(3/beta)); % s.t. $E[S_i^2]=1$
wv = ones(pn,1); % coef vector of LMMSE bound

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
    sX_min = -real(lambertw(0,-exp(-(1+2*eps)))*s0);
    l = sX_min;   % coe
    lv(cnt,1)=l;
    
    % mmse upper bound
    sX_max = -real(lambertw(-1,-exp(-(1+2*eps)))*s0);
    u = sX_max;
    uv(cnt,1)=u;
    
    % Cramer Rao bound
    if beta > 0.5
        c = 1/J;%1/(snr+J); 
    else
        c = 0;
    end
    cv(cnt,1) = c;
  
    % Mismatched waterfilling: using LMMSE bound, coe=1
    [etaw(cnt),pw(cnt,:)] = waterfilling(g);
    
    % RPA with Cramer-Rao bound
    [etac(cnt),pc(cnt,:)] = waterfilling(g.*c);

    % RPA
    [etal(cnt),pl(cnt,:)] = waterfilling(g.*l);
    [etau(cnt),pu(cnt,:)] = waterfilling(g.*u);
    
end

%% plot coefficients of MMSE bounds

figure
hold on
plot(betav,wv,':k')
plot(betav,cv,'--k')
plot(betav,lv,'-k')
plot(betav,uv,'-.k')
xlabel('beta');
ylabel('MMSE bound coefficient');
legend('LMMSE bound','Cramer-Rao bound','the lower MMSE bound','the upper MMSE bound');
legend('location','north');


%% plot water level & floor
al = {'MM';'CR';'RL';'RU'};
for cnt=1:pn
    n = nv(cnt);
    if n == 0 
        wlw = 1/etaw(cnt)*ones(1,k); %water level
        wlc = 1/ etac(cnt)*ones(1,k);
        wll = 1/ etal(cnt)*ones(1,k);
        wlu = 1/ etau(cnt)*ones(1,k);
        fw = 1/etaw(cnt) - pw(cnt,:); %floor surface
        fc = 1/etac(cnt) - pc(cnt,:);
        fl = 1/etal(cnt) - pl(cnt,:);
        fu = 1/etau(cnt) - pu(cnt,:);
        figure
        b(1,:) = bar([wlw; wlc; wll; wlu],1,'FaceColor',[1 1 1],'EdgeColor',[0 .9 .9],'LineWidth',1.5,'Displayname','power p_i');
        hold on
        b(2,:) = bar([fw; fc; fl; fu],1,'FaceColor',[.6 .6 .6],'Displayname','noise 1/(gamma_i x_i)');
        legend(b(:,1),'location','northwest')
        set(gca,'xticklabel',al)
    end
end

%% plot water-level v.s. beta
figure
hold on
plot(betav,1./etaw,':k')
plot(betav,1./etac,'--k')
plot(betav,1./etal,'-k')
plot(betav,1./etau,'-.k')
xlabel('beta');
ylabel('Water level 1/eta');
legend('1/eta_w','1/eta_c','1/eta_l','1/eta_u');
legend('location','south');

%% plot difference of water-levels and c_W v.s. beta
cw = sum(((1./lv - 1./uv)*ones(1,k))./(ones(pn,1)*g),2); % upper bound of the difference between water levels

figure
hold on
plot(betav,cw,'-k')
plot(betav,1./etal-1./etaw,'-.k')
plot(betav,1./etal-1./etau,'--k')
plot(betav,max(1./etac-1./etaw,0),':k')
xlabel('beta');
ylabel('Difference between water levels');
legend('c_W','1./eta_l-1/eta_w','1/eta_l-1/eta_u','1./eta_c-1/eta_w');
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
crul = log2(uv./lv)/2;    % using the upper and lower KL-ball bounds
crwl = log2(wv./lv)/2;    % using LMMSE bound and lower KL-ball bound
crwc = log2(wv./cv)/2;    % using LMMSE bound and Cramer-Rao bound

% an upper bound of the optimal rate
Rub = apprx_mi(wv*ones(1,k), ones(pn,1)*g, pw);    

figure
hold on
plot(betav,crwl,'-k')
plot(betav,Rub-Rw,'-.k')
plot(betav,Rub-Rl,'-.k')
xlabel('beta');
ylabel('Difference between achievable rates R_{ub}-R(p)');
legend('c_R','R_{ub}-R(p_w)','R_{ub}-R(p_l)');
legend('location','south');


figure
hold on
plot(betav,crul','-.k')
plot(betav,crwl','-k')
plot(betav,crwc','--k')
xlabel('beta');
ylabel('c_R');
legend('lu','wl','wc');

