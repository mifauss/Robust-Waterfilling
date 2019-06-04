% fixed non-identical GG input, [beta_i] = [1,2,4,infty]
% varying identical channel gains, h_i=h varying from 0.5 to 2

clc; clear all;

k = 4;      % number of channels, i = 1,...,k
s0 = 1; % E[si^2]=1
P = 10; % power constraint
betav = [1, 2, 4];  % vector of $\beta_i$, i=1,...,k-1, and $\beta_k = \infty$
hv = 0.5:0.05:2;    % channel gain scaler

%% Coefficients of MMSE bounds 
wv = ones(1,k);% LMMSE bound 1/(1+snr)

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
    sX_min = -real(lambertw(0,-exp(-(1+2*epsG)))*s0);
    %     lb = sX_min*sN/(sX_min + sN); % bound
    lv(bn) = sX_min;   % coe
    
    % mmse upper bound
    sX_max = -real(lambertw(-1,-exp(-(1+2*epsG)))*s0);
    %     ub = sX_max*sN/(sX_max + sN);
    uv(bn) = sX_max;

end

% Uniform: Sk ~ U_K(R)
K = 1;	% dimension
R = 1;  % radius, s.t. E[Si^2]=1

% Reference distribution N(0,s0),s0 = 1
epsU = get_KL_div_U(K); % KL ball radius: epsilon

% mmse lower bound
sX_min = -real(lambertw(0,-exp(-(1+2*epsU/K)))*s0);
lv(k) = sX_min;   % coe

% mmse upper bound
sX_max = -real(lambertw(-1,-exp(-(1+2*epsU/K)))*s0);
uv(k) = sX_max;


%% RPA
[px,pn] = size(hv);
for cnt = 1:pn
    h = hv(cnt).*ones(1,k); % identical channel gains, $h_i = h$
    g = P*h.^2;             % $\gamma_i$: strength of the channels
    gv(cnt,:) = g;
    
    % RPA
    [etaw(cnt),pw(cnt,:)] = waterfilling(g.*wv);% Mismatched waterfilling: using LMMSE bound
    [etal(cnt),pl(cnt,:)] = waterfilling(g.*lv);
    [etau(cnt),pu(cnt,:)] = waterfilling(g.*uv);
    
end

%% plot water level & floor
al = {'MM';'RL';'RU'};
for cnt=1:pn
    if cnt == 1 || cnt == pn
        wlw = 1/etaw(cnt)*ones(1,k); %water level
        wll = 1/ etal(cnt)*ones(1,k);
        wlu = 1/ etau(cnt)*ones(1,k);
        fw = 1/etaw(cnt) - pw(cnt,:); %floor surface
        fl = 1/etal(cnt) - pl(cnt,:);
        fu = 1/etau(cnt) - pu(cnt,:);
        figure
        b(1,:) = bar([wlw; wll; wlu],1,'FaceColor',[1 1 1],'EdgeColor',[0 .9 .9],'LineWidth',1.5,'Displayname','power p_i');
        hold on
        b(2,:) = bar([fw; fl; fu],1,'FaceColor',[.6 .6 .6],'Displayname','noise 1/(gamma_i x_i)');
        legend(b(:,1),'location','northwest')
        set(gca,'xticklabel',al)
    end
end

%% plot water-level v.s. beta
figure
hold on
plot(nv,1./etaw,':k')
plot(nv,1./etal,'-k')
plot(nv,1./etau,'-.k')
xlabel('Channel gains h_i');
ylabel('Water level 1/eta');
legend('1/eta_w','1/eta_l','1/eta_u');
legend('location','south');
    
%% plot difference of water-levels and c_W v.s. h_i
cw = sum((ones(pn,1)*(1./lv - 1./uv))./gv,2);% upper bound of the difference between water levels

figure
hold on
plot(hv,cw,'-r')
plot(hv,1./etal-1./etaw,':k')
plot(hv,1./etal-1./etau,'-.k')
xlabel('h_i');
ylabel('Difference between water levels');
legend('c_W','1./eta_l-1/eta_w','1/eta_l-1/eta_u');
legend('location','south');

%% plot ralative rate w.r.t. Rw, the achievable rate of legacy waterfilling
Rw = sum_mi_niidGG(betav,hv,P*pw);
Rl = sum_mi_niidGG(betav,hv,P*pl);
Ru = sum_mi_niidGG(betav,hv,P*pu);

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
crul = 1/2/k*sum(log2(uv./lv))*ones(1,pn);
crwl = 1/2/k*sum(log2(wv./lv))*ones(1,pn);

figure
hold on
plot(hv,crul','-.k')
plot(hv,crwl','-k')
xlabel('beta');
ylabel('c_R');
legend('lu','wl');

% an upper bound of the optimal rate
Rub = apprx_mi(ones(pn,1)*wv, gv, pw);

figure
hold on
plot(hv,crwl,'-k')
plot(hv,Rub-Rl,'--k')
xlabel('h_i');
ylabel('Difference between achievable rates R_{ub}-R(p)');
legend('c_R','R_{ub}-R(p_l)');
legend('location','south');


