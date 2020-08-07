function [ OptimalValue, UplinkRate_PerUser, DownlinkRate_PerUser, W, p, mu, vartheta ] = Get_optSolutionPerIteration2( Pbs, P, Hu_tilde, Hd_tilde, G_SI_tilde, G_lk, error_lk, rho, Fixed_timegroup_assignment, Rate_Threshold, Omega, W_current, p_current, mu_current, vartheta_current, Init, max_min, eta )
%GET_OPTSOLUTIONPERITERATION Summary of this function goes here
%   Detailed explanation goes here

if (nargin<18)
    eta = 1;
end

N = floor(size(G_SI_tilde,1)/2);
G = 2;
L = size(G_lk,1);
K = size(G_lk,2);


beta = [(Omega(1,1) & Omega(2,1)), (Omega(1,2) & Omega(2,2))];
chi = [~(Omega(1,1) | Omega(2,1)), ~(Omega(1,2) | Omega(2,2))];


% sigma_K = 0.01*ones(1, K);
sigma_K = sqrt(10^(-174/10)*10^(7))*ones(1, K);
% Rate_Threshold = Threshold;

% BarrierValue = 100;

W_next = sdpvar(2*N*K,G,'full','complex');


p_next = sdpvar(L, G, 'full','real');

if (Fixed_timegroup_assignment==1)
    
    mu_next = G*ones(1,G);
    
elseif (Fixed_timegroup_assignment==2)
    
    secondtimegroup = 10^-5;
    
    mu_next = [(1/(1-secondtimegroup)) (1/secondtimegroup)];
    
else

    mu_next = sdpvar(1, G, 'full', 'real');

end

Theta_next = sdpvar(K, G, 'full', 'real');

vartheta_next = sdpvar(K, G, 'full', 'real');

varrho = sdpvar(1, 1);

phi = sdpvar(1, 1);

phi_up = sdpvar(1, 1);

phi_down = sdpvar(1, 1);



obj = 0;
cons = [];

for j = 1:1:G
%     Wj_next = [];
           
     
    
     
        for i = 1:1:2
            if (Omega(i,j)==0)
                for k = 1:1:K
                    cons = [cons, cone(W_next(((k-1)*2*N+(i-1)*N+1):((k-1)*2*N+i*N),j), 10^(-10) ) ];
                end
            end
        end
        
        
%         if (Omega(i,j)==0)
%             Wj_next = [Wj_next; zeros(N*K,1)];
            
%         else
%             Wj_next = [Wj_next; sdpvar(N*K,1,'full','complex')];
%         end

%     W_next = [W_next Wj_next];
%     if (sum(Omega(:,j))==2)
%         for l = 1:1:L
%             cons = [cons, p_next(l,j) <= 10^(-10)];
%         end
%     end
end


%% add objective function and constraints for UPLINK rate in Eq. (23)

s_next = sdpvar(L, G, 'full', 'real'); % varphi_vec <= s_next


PerGroupUplinkRate = sdpvar(L,G, 'full', 'real');

tempsum_up = sdpvar(L,1,'full', 'real');


for l = 1:1:L
    
    temp_up = 0;
    
    for g = 1:1:G % g~j
    
    
			
		% get parameters in Eq. (24)

		[zeta, varpi, varphi_vec, uprate] = GetParameter_Uplink(Hu_tilde(:,l:L,g), G_SI_tilde(:,:,g), W_current(:,g), p_current(l:L,g), W_next(:,g), p_next(l:L,g), rho);

        
		% add first term in Eq. (23)
		

%         PerGroupUplinkRate(l,g) = real(zeta)/mu_current(g);

        temp_up = temp_up + real(zeta)/mu_current(g);
		
		
		% add second term in Eq. (23)
		

%         PerGroupUplinkRate(l,g) = PerGroupUplinkRate(l,g) + real(varpi)/mu_current(g)*p_next(l,g);

        temp_up = temp_up + real(varpi)/mu_current(g)*p_next(l,g);


		% add third term in Eq. (23)
		

%         PerGroupUplinkRate(l,g) = PerGroupUplinkRate(l,g) - s_next(l,g)/mu_current(g);

        temp_up = temp_up - s_next(l,g)/mu_current(g);
        
        
        
        % add fourth term in Eq. (23)
        
%         PerGroupUplinkRate(l,g) = PerGroupUplinkRate(l,g) - uprate/(mu_current(g)^2)*mu_next(g);

        temp_up = temp_up - uprate/(mu_current(g)^2)*mu_next(g);
        
        PerGroupUplinkRate(l,g) = real(zeta)/mu_current(g) + real(varpi)/mu_current(g)*p_next(l,g) - s_next(l,g)/mu_current(g) - uprate/(mu_current(g)^2)*mu_next(g);
		
		cons = [cons, cone([varphi_vec; 0.5*(s_next(l,g)-1)], 0.5*(s_next(l,g)+1) ) ];
        
        cons = [cons, s_next(l,g)>=10^(-20)];
        
%         if (~Init)

        cons = [cons, PerGroupUplinkRate(l,g)>=0];
        
%         end
        
    end
    

    if (Init || max_min==1)
        
%         cons = [cons, real(temp_up/Rate_Threshold)>= varrho];
%         cons = [cons, real(PerGroupUplinkRate(l,g)/Rate_Threshold)>= varrho];
        cons = [cons, real(sum(PerGroupUplinkRate(l,:))/Rate_Threshold)>= varrho];
        
    else
        
        if (max_min==0)

%             cons = [cons, real(temp_up/log(2))>=Rate_Threshold];
%             cons = [cons, temp_up>=tempsum_up(l)];
            
            cons = [cons, real(sum(PerGroupUplinkRate(l,:))/log(2))>=Rate_Threshold];
            cons = [cons, sum(PerGroupUplinkRate(l,:))>=tempsum_up(l)];

        else
            
%             cons = [cons, real(temp_up/Rate_Threshold)>= phi_up];
%             cons = [cons, real(temp_up/Rate_Threshold)>= phi];
            cons = [cons, real(sum(PerGroupUplinkRate(l,:))/Rate_Threshold)>= phi];
            
        end
        
    end
    
	
end


obj = obj + sum(tempsum_up);


%% add objective function and constraints for DOWNLINK rate in Eq. (31)

PerGroupDownlinkRate = sdpvar(K,G,'full','real');

tempsum_down = sdpvar(K,1,'full','real');

for k = 1:1:K
    
    temp_down = 0;
    
    for g = 1:1:G
          
		
		% get parameters in Eq. (16)

        [nu, xi, lambda] = GetParameter_Downlink(vartheta_current(k, g), mu_current(g));
		
		% add first term in Eq. (31)
		
%         PerGroupDownlinkRate(k,g) = real(nu);

        temp_down = temp_down + real(nu);
		
		
		% add second term in Eq. (31)
		

%         PerGroupDownlinkRate(k,g) = PerGroupDownlinkRate(k,g) + real(xi)*vartheta_next(k,g);

        temp_down = temp_down + real(xi)*vartheta_next(k,g);
        
        
        % add third term in Eq. (31)
        
%         PerGroupDownlinkRate(k,g) = PerGroupDownlinkRate(k,g) + real(lambda)*mu_next(g);

        temp_down = temp_down + real(lambda)*mu_next(g);
        
        PerGroupDownlinkRate(k,g) = real(nu) + real(xi)*vartheta_next(k,g) + real(lambda)*mu_next(g);
        
%         if (~Init)
        
        cons = [cons, PerGroupDownlinkRate(k,g)>=0];
        
%         end

		        
    end


    if (Init || max_min==1)
        
%         cons = [cons, real(temp_down/Rate_Threshold)>= varrho];
        cons = [cons, real(sum(PerGroupDownlinkRate(k,:))/Rate_Threshold)>= varrho];
        
    else
        
        if (max_min==0)

%             cons = [cons, real(temp_down/log(2))>=Rate_Threshold];
%             cons = [cons, temp_down>=tempsum_down(k)];
            
            cons = [cons, real(sum(PerGroupDownlinkRate(k,:))/log(2))>=Rate_Threshold];
            cons = [cons, sum(PerGroupDownlinkRate(k,:))>=tempsum_down(k)];
            
        else
            
%             cons = [cons, real(temp_down/Rate_Threshold)>= phi_down];
%             cons = [cons, real(temp_down/Rate_Threshold)>= eta*phi];
%             disp('hereeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee');

            cons = [cons, real(sum(PerGroupDownlinkRate(k,:))/Rate_Threshold)>= eta*phi];
            
        end
        
    end
    
    
end


obj = obj + sum(tempsum_down);



%% Add power constraints

for l = 1:1:L
    for g = 1:1:G
%         cons = [cons, real(vec(p_next))>=0]; %11f
        cons = [cons, real(p_next(l,g))>=10^(-20)]; %11f
    end
end

% Constraints 20 21

if (Fixed_timegroup_assignment~=2)
    
mu_bar_next = sdpvar(1, G, 'full', 'real'); 

for g = 1:1:G
    
    cons = [cons, mu_next(g)>=1+10^(-4)];
    
    cons = [cons, mu_next(g)<=10000];
    
    cons = [cons, cone([1 0.5*(mu_next(g)-mu_bar_next(g))], 0.5*(mu_next(g)+mu_bar_next(g))) ];
    
    cons = [cons, mu_bar_next(g)>=10^(-4)];
%     cons = [cons, mu_bar_next(g)<=1-10^(-2)];
    
end

cons = [cons, mu_bar_next(1)+mu_bar_next(2)<=1];

end

% Constraints 26 27 29 30

for g = 1:1:G
    
    if ( ((chi(1)) && (g==1)) || ((chi(2)) && (g==2)) )
        cons = [cons, vartheta_next(:,g)>=10^(-10)];
        continue
    end
    
    W_next_mat_all = reshape(W_next(:,g), 2*N, K);
    
    for k = 1:1:K
    
        W_next_mat = W_next_mat_all(:,[1:(k-1) (k+1):K]);  
        
        % Constraints 26 27
        
%         cons = [cons, real(Hd_tilde(:,k,g)'*W_next(((k-1)*2*N+1):(k*2*N),g)) >= 0 ];
%         cons = [cons, real(Hd_tilde(:,k,g)'*W_current(((k-1)*2*N+1):(k*2*N),g))*(2*real(Hd_tilde(:,k,g)'*W_next(((k-1)*2*N+1):(k*2*N),g))-(real(Hd_tilde(:,k,g)'*W_current(((k-1)*2*N+1):(k*2*N),g))))>=0]; %xxxxxxxxxxxxxxxxxxxxxxxxxxxx
        
%         cons = [cons, 2*real(Hd_tilde(:,k,g)'*W_current(((k-1)*2*N+1):(k*2*N),g))*real(Hd_tilde(:,k,g)'*W_next(((k-1)*2*N+1):(k*2*N),g))-(real(Hd_tilde(:,k,g)'*W_current(((k-1)*2*N+1):(k*2*N),g)))^2>=Theta_next(k,g) ];
        Theta_next(k,g) = 2*real(Hd_tilde(:,k,g)'*W_current(((k-1)*2*N+1):(k*2*N),g))*real(Hd_tilde(:,k,g)'*W_next(((k-1)*2*N+1):(k*2*N),g))-(real(Hd_tilde(:,k,g)'*W_current(((k-1)*2*N+1):(k*2*N),g)))^2;
        cons = [cons, Theta_next(k,g)>= 0];
        
        
        % Constraints 29 30
        
%         cons = [cons, vartheta_next(k,g)>0];
        
        Sum_vec = [(Hd_tilde(:,k,g)'*W_next_mat), (p_next(:,g).*G_lk(:,k))', (p_next(:,g).*sqrt(error_lk(:,k)))', sigma_K(k)];
		
        cons = [cons, cone([Sum_vec, 0.5*(Theta_next(k,g)-vartheta_next(k,g))], 0.5*(Theta_next(k,g)+vartheta_next(k,g)) )];
        
%         cons = [cons, vartheta_next(k,g)>=10^-(1.25)];
        cons = [cons, vartheta_next(k,g)>=10^(-10)];
        
%         cons = [cons, vartheta_next(k,g)<=10^(10)];
        
        
    end
end


%% Constraints 34a 34b

if (chi(2))
    PW1_next = sdpvar(1, 1, 'full', 'real');
    PW2_next = sdpvar(1, 1, 'full', 'real');
    cons = [cons, PW1_next>=10^(-20)];
    cons = [cons, PW2_next>=10^(-20)];
    cons = [cons, cone([W_next(:,1); 0.5*(PW1_next-1)], 0.5*(PW1_next+1))];
    cons = [cons, cone([W_next(:,2); 0.5*(PW2_next-1)], 0.5*(PW2_next+1))];
    cons = [cons, PW2_next <= 10^(-10)];
    cons = [cons, PW1_next + PW2_next <= Pbs];
    
elseif (chi(1))
%     PW1_next = sdpvar(1, 1, 'full', 'real');
%     PW2_next = sdpvar(1, 1, 'full', 'real');
%     cons = [cons, PW1_next>=10^(-10)];
%     cons = [cons, PW2_next>=10^(-10)];
%     cons = [cons, cone([W_next(:,1); 0.5*(PW1_next-1)], 0.5*(PW1_next+1))];
%     cons = [cons, cone([W_next(:,2); 0.5*(PW2_next-1)], 0.5*(PW2_next+1))];
%     cons = [cons, PW1_next <= 10^(-9)];
%     cons = [cons, PW1_next + PW2_next <= Pbs];

    
%     PW2_next = sdpvar(1, 1, 'full', 'real');
%     cons = [cons, PW2_next>=10^(-10)];
%     cons = [cons, cone(W_next(:,1), 10^(-10))];
%     cons = [cons, cone([W_next(:,2); 0.5*(PW2_next-1)], 0.5*(PW2_next+1))];
%     cons = [cons, PW2_next <= Pbs];
    
    cons = [cons, cone([W_next(:,1); W_next(:,2)], sqrt(Pbs))];
else
%     disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    mu_W1_next = sdpvar(1, 1, 'full', 'real');
    mu_W2_next = sdpvar(1, 1, 'full', 'real');
    cons = [cons, mu_W1_next>=10^(-10)];
    cons = [cons, mu_W2_next>=10^(-10)];
    
%     cons = [cons, cone([W_next(:,1); 0.5*(mu_W1_next-1)], 0.5*(mu_W1_next+1)) ];
    cons = [cons, cone([W_next(:,1); 0.5*(mu_W1_next-mu_next(1))], 0.5*(mu_W1_next+mu_next(1))) ];
    
    cons = [cons, cone([W_next(:,2); 0.5*(mu_W2_next-mu_next(2))], 0.5*(mu_W2_next+mu_next(2))) ];
    
%     cons = [cons, (mu_W1_next+mu_W2_next-2*real(W_current(:,1)'*W_next(:,1))/mu_current(2) + (norm(W_current(:,1))^2)/(mu_current(2)^2)*mu_next(2)) <=  Pbs ];
    cons = [cons, (mu_W1_next+mu_W2_next) <= Pbs ];
end


% 
% 
% cons = [cons, mu_tilde_next(2)>=10^(-10)];
% 
% 
% 
% cons = [cons, cone([W_next(:,2); 0.5*(mu_next(2)-mu_tilde_next(2))], 0.5*(mu_tilde_next(2)+mu_next(2))) ];
% 
% if (chi(2))
%     
%     mu_tilde2_next1 = sdpvar(1, 1, 'full', 'real');
%     cons = [cons, mu_tilde2_next1>=10^(-10)];
%     cons = [cons, cone([W_next(:,1); 0.5*(mu_next(2)-mu_tilde2_next1)], 0.5*(mu_tilde2_next1+mu_next(2)) ) ]; %%%%%
% else
%     mu_tilde2_next1 = 0;
% end
% 
% if (chi(1))
%     mu_tilde2_next2 = sdpvar(1, 1, 'full', 'real');
%     cons = [cons, mu_tilde2_next2>=10^(-10)];
%     cons = [cons, cone([W_next(:,2); 0.5*(mu_tilde2_next2-1)], 0.5*(mu_tilde2_next2+1)) ]; %%%%%%
% else
%     mu_tilde2_next2 = 0;
% end
% 
% if (chi(1))
%     cons = [cons, mu_tilde_next(1)+mu_tilde_next(2)+mu_tilde2_next1+mu_tilde2_next2 <= (2*real(W_current(:,1)'*W_next(:,1))/mu_current(2) + 2*real(W_current(:,2)'*W_next(:,2))/mu_current(2) ...
%                                                             - (norm(W_current(:,1))^2+norm(W_current(:,2))^2)/(mu_current(2)^2)*mu_next(2) + Pbs) ]; %%%%%
% else
%     cons = [cons, mu_tilde_next(1)+mu_tilde_next(2)+mu_tilde2_next1+mu_tilde2_next2 <= (2*real(W_current(:,1)'*W_next(:,1))/mu_current(2) - (norm(W_current(:,1))^2)/(mu_current(2)^2)*mu_next(2) + Pbs) ]; %%%%%
% end
                                             
%%                                             

% mu_hat_next = sdpvar(L, G, 'full', 'real');
% 
% mu_hat2_next = sdpvar(L, G, 'full', 'real');



    
if (beta(2))
    mu_p1_next = sdpvar(1, L, 'full', 'real');
    mu_p2_next = sdpvar(1, L, 'full', 'real');
    for l = 1:1:L
        cons = [cons, cone([p_next(l,1) 0.5*(mu_p1_next(l)-1)], 0.5*(mu_p1_next(l)+1)) ];
        cons = [cons, cone([p_next(l,2) 0.5*(mu_p2_next(l)-1)], 0.5*(mu_p2_next(l)+1)) ];
        cons = [cons, mu_p1_next(l)+mu_p2_next(l)<=P(l) ];
        cons = [cons, mu_p2_next(l)<=10^(-10) ];
        cons = [cons, mu_p1_next(l)>=10^(-20)];
        cons = [cons, mu_p2_next(l)>=10^(-20)];
    end
elseif (beta(1))
    mu_p1_next = sdpvar(1, L, 'full', 'real');
    mu_p2_next = sdpvar(1, L, 'full', 'real');
    for l = 1:1:L
        cons = [cons, cone([p_next(l,1) 0.5*(mu_p1_next(l)-1)], 0.5*(mu_p1_next(l)+1)) ];
        cons = [cons, cone([p_next(l,2) 0.5*(mu_p2_next(l)-1)], 0.5*(mu_p2_next(l)+1)) ];
        cons = [cons, mu_p2_next(l)+mu_p2_next(l)<=P(l) ];
        cons = [cons, mu_p1_next(l)<=10^-(10) ];
        cons = [cons, mu_p1_next(l)>=10^(-20)];
        cons = [cons, mu_p2_next(l)>=10^(-20)];
    end
else
%     disp('bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb');
    mu_p1_next = sdpvar(1, L, 'full', 'real');
    mu_p2_next = sdpvar(1, L, 'full', 'real');
    for l = 1:1:L
        cons = [cons, mu_p1_next(l)>=10^(-10)];
        
%         cons = [cons, cone([p_next(l,1) 0.5*(mu_p1_next(l)-1)], 0.5*(mu_p1_next(l)+1)) ];
        cons = [cons, cone([p_next(l,1) 0.5*(mu_p1_next(l)-mu_next(1))], 0.5*(mu_p1_next(l)+mu_next(1))) ];
        
        cons = [cons, mu_p2_next(l)>=10^(-10)];
        cons = [cons, cone([p_next(l,2) 0.5*(mu_p2_next(l)-mu_next(2))], 0.5*(mu_p2_next(l)+mu_next(2))) ];
        
%         cons = [cons, (mu_p1_next(l) + mu_p2_next(l) - 2*p_current(l,1)/mu_current(2)*p_next(l,1) + p_current(l,1)^2/(mu_current(2)^2)*mu_next(2) ) <= P(l) ];
        cons = [cons, (mu_p1_next(l) + mu_p2_next(l)) <= P(l)];
    end

end

%     cons = [cons, mu_hat_next(l,1)>=10^-(10)];
%     cons = [cons, mu_hat_next(l,2)>=10^-(10)];
%     cons = [cons, mu_hat2_next(l,1)>=10^-(10)];
%     cons = [cons, mu_hat2_next(l,2)>=10^-(10)];
%     
%     cons = [cons, cone([p_next(l,1) 0.5*(mu_hat_next(l,1)-1)], 0.5*(mu_hat_next(l,1)+1)) ];
%     
%     cons = [cons, cone([p_next(l,2) 0.5*(mu_next(2)-mu_hat_next(l,2))], 0.5*(mu_hat_next(l,2)+mu_next(2))) ];
%     
%     if (beta(2))
%         cons = [cons, cone([p_next(l,1) 0.5*(mu_next(2)-mu_hat2_next(l,1))], 0.5*(mu_hat2_next(l,1)+mu_next(2))) ]; %%%%%
%     else
%         cons = [cons, mu_hat2_next(l,1)<=10^-(5)];
%     end
%     
%     if (beta(1))
%         cons = [cons, cone([p_next(l,2) 0.5*(mu_hat2_next(l,2)-1)], 0.5*(mu_hat2_next(l,2)+1)) ]; %%%%%
%     else
%         cons = [cons, mu_hat2_next(l,2)<=10^-(5)];
%     end
%     
%     if (beta(1))
% 
%         cons = [cons, mu_hat_next(l,1)+mu_hat_next(l,2)+mu_hat2_next(l,1)+mu_hat2_next(l,2) <= (2*p_current(l,1)/mu_current(2)*p_next(l,1) + 2*p_current(l,2)/mu_current(2)*p_next(l,2) ...
%                                                                                                 - (p_current(l,1)^2+p_current(l,2)^2)/(mu_current(2)^2)*mu_next(2) + P(l) ) ];
% 
%     else
%         cons = [cons, mu_hat_next(l,1)+mu_hat_next(l,2)+mu_hat2_next(l,1)+mu_hat2_next(l,2) <= (2*p_current(l,1)/mu_current(2)*p_next(l,1) - p_current(l,1)^2/(mu_current(2)^2)*mu_next(2) + P(l) ) ];
%     end
%                        
% end




if (Init || max_min==1)
    
%     myops = sdpsettings('solver','mosek');
    myops = sdpsettings('solver','sdpt3','verbose',0);
    diagnotics = solvesdp(cons, -varrho, myops)
    
    OptimalValue = double(varrho);
    
else
    
    if (max_min==0)

    %     myops = sdpsettings('solver','mosek');
        myops = sdpsettings('solver','sdpt3','verbose',0);
        diagnotics = solvesdp(cons, -obj, myops);

        OptimalValue = double(obj);
        sumrate_up = double(sum(tempsum_up));
        sumrate_down = double(sum(tempsum_down));
    else
        
%         phi = phi_up + phi_down;
        
        %     myops = sdpsettings('solver','mosek');
        myops = sdpsettings('solver','sdpt3','verbose',0);
        diagnotics = solvesdp(cons, -phi, myops)

        OptimalValue = double(phi);
        
    end
    
end




DownlinkRate_PerUser = real(double(PerGroupDownlinkRate));



UplinkRate_PerUser = real(double(PerGroupUplinkRate));




W = double(W_next);

TraceW = trace(W*W');

p = double(p_next);
PP = p*p';

mu = double(mu_next);
vartheta = double(vartheta_next);

Theta = double(Theta_next);



end

