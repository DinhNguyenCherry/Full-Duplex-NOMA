function [ W_current ] = Get_InitialW0( G, H, Pbs, HalfDuplex )
%GET_INITIALW0 Summary of this function goes here
%   Detailed explanation goes here

% global N_tx
% global K
% global G

K = size(H,2);

tau_current = G*ones(1,G);
if (nargin<4 || HalfDuplex==0)
    N_tx = size(H,1);
    W_next = sdpvar(N_tx*K,G,'full','complex');
    HalfDuplex = 0;
else
    N_tx = floor(size(H,1)/2);
    N_rx = size(H,1) - N_tx;
    
    W_next = sdpvar((N_tx+N_rx)*K,G,'full','complex');
    tau_current = G*ones(1,G);
end
mu_next = sdpvar(1, G); %Pbs/G*ones(1,G);%



obj_W = 0;
cons_W = [];

for g = 1:1:G

	for k = 1:1:K
		

		if (HalfDuplex)
            cons_W = [cons_W, real(H(:,k,g)'*W_next(((k-1)*(N_tx+N_rx)+1):k*(N_tx+N_rx),g))>=0]; % Constraint 12b
        else
            cons_W = [cons_W, real(H(:,k,g)'*W_next(((k-1)*N_tx+1):k*N_tx,g))>=0]; % Constraint 12b
        end
		
    end
    
    cons_W = [cons_W, cone([W_next(:,g); 0.5*(mu_next(g)-tau_current(g))],0.5*(mu_next(g)+tau_current(g)))]; % Constraint 39b
	
end


cons_W = [cons_W, sum(mu_next) <= Pbs]; % Constraint 39a

myops = sdpsettings('solver','sdpt3','verbose',0);

diagnotics = solvesdp(cons_W, -obj_W, myops);

OptimalValue_W0 = double(obj_W);


W_current = double(W_next);

mu_current = double(mu_next);

end

