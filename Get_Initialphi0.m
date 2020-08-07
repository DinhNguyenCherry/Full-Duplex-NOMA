function [ phi_current ] = Get_Initialphi0( H, G_hat, W_current, p_current, HalfDuplex )
%GET_INITIALPHI0 Summary of this function goes here
%   Detailed explanation goes here

% global G
% global K
% global N_tx
% global sigma_K

G = size(W_current,2);
K = size(G_hat,2);

sigma_K = 0.01*ones(1, K);

if (nargin<5 || HalfDuplex==0)
    HalfDuplex = 0;
    N_tx = size(H,1);
else
%     global N_rx
    N_tx = floor(size(H,1)/2);
    N_rx = size(H,1)-N_tx;
end

% phi_current = zeros(K, G);

for g = 1:1:G
    
    if (HalfDuplex)
        W_all_k = reshape(W_current(:, g), N_tx+N_rx, K);
    else
        W_all_k = reshape(W_current(:, g), N_tx, K);
    end
	
	for k = 1:1:K
	
		W_without_k = W_all_k(:, [1:(k-1) (k+1):K]);
        
        temp = 0;
		
		temp = temp + (H(:,k)'*W_without_k*W_without_k'*H(:,k));
		
		temp = temp + sum((p_current(:,g).^2).*(abs(G_hat(:,k)).^2));
		
		temp = temp + sigma_K(k)^2;
        
        phi_current(k,g) = sqrt(temp);
		
	end
	
end

end

