% % % The MATLAB CODE is used for the following paper: Hieu V. Nguyen, Van-Dinh Nguyen, Octavia A. Dobre, Diep N. Nguyen, Eryk Dutkiewicz, and Oh-Soon Shin, 
% % % "Joint Power Control and User Association for NOMA-Based Full-Duplex Systems,"
% % % IEEE Transactions on Communications, vol. 67, no. 11, pp. 8037-8055, Nov. 2019.


function [ GetBreak, OptimalValue, OptimalValue_preStep, UplinkRate_PerUser, DownlinkRate_PerUser, W_current, p_current, addVar_current, PenaltyChain ] = GetInitialization( Pbs, P, Hu, Hd, G_SI, G_lk, rho, Rate_Threshold, max_min, eta, InitInfo, InitPoint, zone3 )
%GETINITIALIZATION Summary of this function goes here
%   Detailed explanation goes here

if (nargin<nargin('GetInitialization'))
    zone3 = 0;
end

Init = InitInfo{1};
Link = InitInfo{2};

N = size(G_SI,1);
L = size(G_lk,1);
if (zone3)
    K = floor(size(G_lk,2)/3);
else
K = floor(size(G_lk,2)/2);
end
sigma_K = sqrt(10^(-174/10)*10^(7))*ones(1, K);

if (max_min)
    FesiPoint = 0.5;
else
    FesiPoint = 0.2;
end

% if (Fixed_timegroup_assignment==2)
%     G=1;
% end

OptimalValue_preStep = 0;
OptimalValue = 0;

n=0;
GetBreak = 0;

if (Init==1 || (Init==2 && strcmp(Link,'ExSearch')))

W1_current = 10*(rand(N*K,1)+1i*rand(N*K,1));
W2_current = 10*(rand(N*K,1)+1i*rand(N*K,1));
W3_current = 10*(rand(N*K,1)+1i*rand(N*K,1));
    
p_current = 2*rand(L, 1);%zeros(L, G);
if (~isempty(strfind(Link,'UL')))
%     addVar_current = 
    if (zone3)
        alpha_current = InitPoint{1,3}{1,1}(1:K,:);
        alpha_current2 = InitPoint{1,3}{1,1}(K+1:2*K,:);
        alpha_current3 = InitPoint{1,3}{1,1}(2*K+1:3*K,:);
    else
        alpha_current = InitPoint{3}{1};
    end
elseif (strcmp(Link,'ExSearch'))
    alpha = InitPoint{1}
    alpha_idx = alpha + (alpha>0).*([0:(K-1)]*K);
    alpha_idx = alpha_idx(find(alpha_idx>0));
    alpha_current = zeros(K,K);
    alpha_current(alpha_idx) = 1;
else
    alpha_current = rand(K,K);%0.5*ones(K,K);%
    if (zone3)
    alpha_current2 = rand(K,K);
    alpha_current3 = alpha_current'*alpha_current2;
    end
end

if (strcmp(Link,'ExSearch'))
    beta_order = InitPoint{2};
    if (sum(beta_order)==0)
        beta_current = 1-eye(L);
    else
    beta_current = zeros(L,L);
    for l = 1:1:L
        for ll = l+1:1:L
            beta_current(beta_order(l),beta_order(ll)) = 1;
        end
    end
    end
else
    beta_current = rand(L,L).*(1-eye(L));%0.5*ones(L,L);%
end
% beta_current = zeros(L,L);
% for l = 1:1:L
%     for ll = 1:1:L
%         if (l<ll)
%             beta_current(l,ll) = 1-10^(-5);
%         else
%             beta_current(l,ll) = 10^(-5);
%         end
%     end
% end

if (zone3)
    omega_current = zeros(K, 3);
else
omega_current = zeros(K, 2);
end

if (zone3)
    for zone = 1:3
        if (zone==1)
            W1_current_mat = reshape(W1_current, N, K);
            for k = 1:1:K
                phi_k = sum(abs(Hd(:,k)'*W1_current_mat(:,[1:(k-1) (k+1):K])).^2) + ...
                        sum((1-alpha_current(k,:)).*abs(Hd(:,k)'*reshape(W2_current, N, K)).^2) + ...
                        sum((1-alpha_current2(k,:)).*abs(Hd(:,k)'*reshape(W3_current, N, K)).^2) + ...
                        sum(p_current.^2.*G_lk(:,k).^2) + sigma_K(k)^2;
                omega_current(k,zone) = real(phi_k/(abs(Hd(:,k)'*W1_current_mat(:,k)).^2));
            end
        elseif (zone==2)
            W2_current_mat = reshape(W2_current, N, K);
            for j = 1:1:K
                varphi_j = sum(abs(Hd(:,K+j)'*reshape(W1_current, N, K)).^2) + ...
                           sum(abs(Hd(:,K+j)'*W2_current_mat(:,[1:(j-1) (j+1):K])).^2) + ...
                           sum((1-alpha_current3(j,:)).*abs(Hd(:,K+j)'*reshape(W3_current, N, K)).^2) + ...
                           sum(p_current.^2.*G_lk(:,K+j).^2) + sigma_K(j)^2;
                minSINR = real((abs(Hd(:,K+j)'*W2_current_mat(:,j)).^2)/varphi_j);
                for k = 1:1:K
                    psi_kj = sum(abs(Hd(:,k)'*reshape(W1_current, N, K)).^2) + ...
                           sum(abs(Hd(:,k)'*W2_current_mat(:,[1:(j-1) (j+1):K])).^2) + ...
                           sum((1-alpha_current2(k,:)).*abs(Hd(:,k)'*reshape(W3_current, N, K)).^2) + ...
                           sum(p_current.^2.*G_lk(:,k).^2) + sigma_K(k)^2;
                    minSINR = min(minSINR, real((abs(Hd(:,k)'*W2_current_mat(:,j)).^2)/((alpha_current(k,j)+10^(-3))*psi_kj) ));
                end
                omega_current(j,zone) = 1/minSINR;
            end
        else
            W3_current_mat = reshape(W3_current, N, K);
            for j = 1:1:K
                varphi_j = sum(abs(Hd(:,2*K+j)'*reshape(W1_current, N, K)).^2) + ...
                           sum(abs(Hd(:,2*K+j)'*reshape(W2_current, N, K)).^2) + ...
                           sum(abs(Hd(:,2*K+j)'*W3_current_mat(:,[1:(j-1) (j+1):K])).^2) + ...
                           sum(p_current.^2.*G_lk(:,2*K+j).^2) + sigma_K(j)^2;
                minSINR = real((abs(Hd(:,K+j)'*W3_current_mat(:,j)).^2)/varphi_j);
                
                for k = 1:1:K
                    psi_kj = sum(abs(Hd(:,k)'*reshape(W1_current, N, K)).^2) + ...
                             sum(abs(Hd(:,k)'*reshape(W2_current, N, K)).^2) + ...
                             sum(abs(Hd(:,k)'*W3_current_mat(:,[1:(j-1) (j+1):K])).^2) + ...
                             sum(p_current.^2.*G_lk(:,k).^2) + sigma_K(k)^2;
                    minSINR = min(minSINR, real((abs(Hd(:,k)'*W3_current_mat(:,j)).^2)/((alpha_current2(k,j)+10^(-3))*psi_kj) ));
                end
                
                for k = 1:1:K
                    psi_kj = sum(abs(Hd(:,K+k)'*reshape(W1_current, N, K)).^2) + ...
                             sum(abs(Hd(:,K+k)'*reshape(W2_current, N, K)).^2) + ...
                             sum(abs(Hd(:,K+k)'*W3_current_mat(:,[1:(j-1) (j+1):K])).^2) + ...
                             sum(p_current.^2.*G_lk(:,K+k).^2) + sigma_K(k)^2;
                    minSINR = min(minSINR, real((abs(Hd(:,K+k)'*W3_current_mat(:,j)).^2)/((alpha_current3(k,j)+10^(-3))*psi_kj) ));
                end
                
                omega_current(j,zone) = 1/minSINR;
            end
        end
    end
else
    for zone = 1:2
        if (zone==1)
            W1_current_mat = reshape(W1_current, N, K);
            for k = 1:1:K
                phi_k = sum(abs(Hd(:,k)'*W1_current_mat(:,[1:(k-1) (k+1):K])).^2) + ...
                        sum((1-alpha_current(k,:)).*abs(Hd(:,k)'*reshape(W2_current, N, K)).^2) + ...
                        sum(p_current.^2.*G_lk(:,k).^2) + sigma_K(k)^2;
                omega_current(k,zone) = real(phi_k/(abs(Hd(:,k)'*W1_current_mat(:,k)).^2));
            end
        else
            W2_current_mat = reshape(W2_current, N, K);
            for j = 1:1:K
                varphi_j = sum(abs(Hd(:,K+j)'*reshape(W1_current, N, K)).^2) + ...
                           sum(abs(Hd(:,K+j)'*W2_current_mat(:,[1:(j-1) (j+1):K])).^2) + ...
                           sum(p_current.^2.*G_lk(:,K+j).^2) + sigma_K(j)^2;
                minSINR = real((abs(Hd(:,K+j)'*W2_current_mat(:,j)).^2)/varphi_j);
                for k = 1:1:K
                    psi_kj = sum(abs(Hd(:,k)'*reshape(W1_current, N, K)).^2) + ...
                           sum(abs(Hd(:,k)'*W2_current_mat(:,[1:(j-1) (j+1):K])).^2) + ...
                           sum(p_current.^2.*G_lk(:,k).^2) + sigma_K(k)^2;
                    minSINR = min(minSINR, real((abs(Hd(:,k)'*W2_current_mat(:,j)).^2)/((alpha_current(k,j)+10^(-3))*psi_kj) ));
                end
                omega_current(j,zone) = 1/minSINR;
            end
        end
    end
end
% omegap_current = 1./omega_current;

lambda_current = 1-alpha_current;
if (zone3)
    lambda_current2 = 1-alpha_current2;
    lambda_current3 = 1-alpha_current3;
end
% lambda_current = 1./alpha_current;

mu_current = zeros(K,K);
if (zone3)
    mu_current2 = zeros(K,K);
    mu_current3 = zeros(K,K);
end

for k = 1:1:K
    W2_current_mat = reshape(W2_current, N, K);
    for j = 1:1:K
        mu_current(k,j) = abs(Hd(:,k)'*W2_current_mat(:,j))^2;
        if (zone3)
            mu_current2(k,j) = abs(Hd(:,k)'*W3_current_mat(:,j))^2;
            mu_current3(k,j) = abs(Hd(:,K+k)'*W3_current_mat(:,j))^2;
        end
    end
end

nu_current = beta_current*diag(p_current.^2);

% u_current = zeros(L,1);
% for l = 1:1:L
%     u_current(l) = sum(beta_current(l,:));
% end

u_current = zeros(L,L);
for l = 1:1:L
    for ll = 1:1:L
        u_current(l,ll) = sum(beta_current(l,:))-sum(beta_current(ll,:));
    end
end

uu_current = sum(beta_current,2);
if (zone3)
    W_current = [W1_current; W2_current; W3_current];
    alpha_all = [alpha_current; alpha_current2; alpha_current3];
    lambda_all = [lambda_current; lambda_current2; lambda_current3];
    mu_all = [mu_current; mu_current2; mu_current3];
    addVar_current = {alpha_all, beta_current, omega_current, ...
                      lambda_all, mu_all, nu_current, u_current, uu_current};
else
    W_current = [W1_current; W2_current];
    addVar_current = {alpha_current, beta_current, omega_current, ...
                      lambda_current, mu_current, nu_current, u_current, uu_current};
end
% addVar_current = {alpha_current, beta_current, omega_current, ...
%                   lambda_current, mu_current, nu_current, u_current, uu_current};
% addVar_current = {alpha_current, beta_current, [omega_current; omegap_current], ...
%                   lambda_current, mu_current, nu_current, u_current};
else
    
   W_current = InitPoint{1}; 
   p_current = InitPoint{2};
   addVar_current = InitPoint{3};

end

MaxInitialIter = 50; 
OptimalValue_current = -1000;
Status = 'Infeasible';
PenaltyChain_alpha = [];
PenaltyChain_beta = [];
       
% while (~(OptimalValue>FesiPoint && isempty(findstr('Infeasible',Status))) && ~(Init==1 && isempty(findstr('Infeasible',Status)) && isempty(findstr('Maxim',Status)) ) )
while (~(OptimalValue>FesiPoint && isempty(findstr('Infeasible',Status))) && ~(Init==1 && isempty(findstr('Infeasible',Status)) ) )
% while (1)%(~(OptimalValue>FesiPoint && isempty(findstr('Infeasible',Status))) )
% while (~(OptimalValue>FesiPoint && isempty(findstr('Infeasible',Status))) && ~(Init==1 && isempty(findstr('Infeasible',Status)) && abs(OptimalValue)<1  ) )
    
    if (n>0)
        OptimalValue_preStep = OptimalValue;
    end

    disp(['---- For rho = ' num2str(10*log10(rho)) ' dB,  Init = ' num2str(Init) ' for ' InitInfo{2}]);
    disp(['Iteration: n = ' num2str(n+1)]);
    
%     [ OptimalValue, DownlinkRate_PerGroupPerUser, UplinkRate_PerGroupPerUser, RDown_current, RThDown_current, RUp_current, RTh_current, W_current, p_current, phi_current, time_current, w_tilde_current, p_bar_current, alpha_bar_current, alpha, beta_bar_current, beta] = Get_optSolutionPerIteration(Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, time_current, w_tilde_current, p_bar_current, RDown_current, RThDown_current, RUp_current, RTh_current, alpha_bar_current, beta_bar_current, rho, Fixed_timegroup_assignment, Rate_Threshold);
    
    if (~isempty(strfind(Link, 'PEN')))
        if (zone3)
            [ OptimalValue, UplinkRate_PerUser, DownlinkRate_PerUser, W_next, p_next, addVar_next, Status] = ... 
                Get_optSolutionPerIteration3_MultiZones(Pbs, P, Hu, Hd, G_SI, G_lk, rho, Rate_Threshold, W_current, p_current, addVar_current, InitInfo, max_min, eta, n);
        else
            [ OptimalValue, UplinkRate_PerUser, DownlinkRate_PerUser, W_next, p_next, addVar_next, Status] = ... 
                Get_optSolutionPerIteration3(Pbs, P, Hu, Hd, G_SI, G_lk, rho, Rate_Threshold, W_current, p_current, addVar_current, InitInfo, 0, eta, n);
        end
        PenaltyChain_alpha = [];%[PenaltyChain_alpha norm(vec(addVar_current{1}.^2-addVar_current{1}),'Inf')];
        PenaltyChain_beta = [];%[PenaltyChain_beta norm(vec(addVar_current{2}.^2-addVar_current{2}),'Inf')];
    else
        if (zone3)
            [ OptimalValue, UplinkRate_PerUser, DownlinkRate_PerUser, W_next, p_next, addVar_next, Status] = ... 
                Get_optSolutionPerIteration3_MultiZones(Pbs, P, Hu, Hd, G_SI, G_lk, rho, Rate_Threshold, W_current, p_current, addVar_current, InitInfo, max_min, eta, n);
        else
        [ OptimalValue, UplinkRate_PerUser, DownlinkRate_PerUser, W_next, p_next, addVar_next, Status] = ... 
        Get_optSolutionPerIteration2(Pbs, P, Hu, Hd, G_SI, G_lk, rho, Rate_Threshold, W_current, p_current, addVar_current, InitInfo, 0, eta);
        end
    end
    
    
    
    if (n>10 && abs(OptimalValue-OptimalValue_current)<0.01 && abs(OptimalValue+Rate_Threshold)<0.001 && abs(OptimalValue)>Rate_Threshold)
        GetBreak = 1;
        disp('######## Exit Initialization ########')
        break;
    end
    OptimalValue
    OptimalValue_current = OptimalValue;
    
    W_current = W_next;
    
    p_current = p_next;
    
%     if (~isempty(findstr('Success',Status)))
%         for k = 1:1:K
%             MVal = max(addVar_next{1}(k,:));
%             alpha_k = addVar_next{1}(k,:);
%             idx = find(alpha_k==MVal);
%             alpha_k = zeros(1,K)+10^(-5);
%             alpha_k(idx) = 1-10^(-5);
%             addVar_next{1}(k,:) = alpha_k;
%         end
%         addVar_next{1}
%     end
    addVar_current = addVar_next;
    
%     if (~isempty(findstr('Success',Status)))
%         break;
%     end
    
    n = n+1;
    
    if (n>MaxInitialIter || (round(OptimalValue,20)==0) || (abs(OptimalValue)>10^8))
%     if (n>MaxInitialIter || (round(OptimalValue,20)==0))
%     if (n>50)
        GetBreak = 1;
        break;
    end
    
end

PenaltyChain = {PenaltyChain_alpha, PenaltyChain_beta};

disp(['#Iterations: ' num2str(n) ]);
OptimalValue_preStep
OptimalValue 
UplinkRate_PerUser
DownlinkRate_PerUser


end

