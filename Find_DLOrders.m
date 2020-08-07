% % % The MATLAB CODE is used for the following paper: Hieu V. Nguyen, Van-Dinh Nguyen, Octavia A. Dobre, Diep N. Nguyen, Eryk Dutkiewicz, and Oh-Soon Shin, 
% % % "Joint Power Control and User Association for NOMA-Based Full-Duplex Systems,"
% % % IEEE Transactions on Communications, vol. 67, no. 11, pp. 8037-8055, Nov. 2019.


function [ alpha_next, maxSumDL ] = Find_DLOrders( Hd, G_SI, G_lk, Rate_Threshold, W_current, p_current )
%FIND_ULORDERS Summary of this function goes here
%   Detailed explanation goes here
N = size(G_SI,1);
L = size(G_lk,1);
K = floor(size(G_lk,2)/2);

PairRate = zeros(K,K);

for k = 1:1:K
    for j = 1:1:K

        PairRate(k,j) = Get_Costmat(Hd, G_lk, Rate_Threshold, W_current, p_current, k, j);
        
    end
end

PairRate

[Sol, maxSumDL]=lapjv(-PairRate);
maxSumDL = -maxSumDL;

alpha_next = zeros(K,K);
pos = Sol + [0:(K-1)]*K;
alpha_next(pos) = 1;

end

function [ Weight ] = Get_Costmat( Hd, G_lk, Rate_Threshold, W_current, p_current, idx_k, idx_j )

K = floor(size(G_lk,2)/2);
N = floor(size(W_current,1)/K/2);

sigma = 10*sqrt(10^(-174/10)*10^(7));

W_current_mat = reshape(W_current, N, 2*K);

InnerIDXinterference = [(1:idx_k-1) ((idx_k+1):(K+idx_j-1)) ((K+idx_j+1):(2*K))];

InnerUserRate = log(1+(abs(Hd(:,idx_k)'*W_current_mat(:,idx_k))^2)/ ...
                    (norm(Hd(:,idx_k)'*W_current_mat(:,InnerIDXinterference))^2 + sum((abs(G_lk(:,idx_k)).*p_current).^2) + sigma^2) );
                
% if (InnerUserRate>=Rate_Threshold)
    OuterIDXinterference = [(1:(K+idx_j-1)) ((K+idx_j+1):(2*K))];
    
    Rate2At1 = (abs(Hd(:,idx_k)'*W_current_mat(:,K+idx_j))^2)/ ...
                    (norm(Hd(:,idx_k)'*W_current_mat(:,OuterIDXinterference))^2 + sum((abs(G_lk(:,idx_k)).*p_current).^2) + sigma^2);
                
    Rate2At2 = (abs(Hd(:,K+idx_j)'*W_current_mat(:,K+idx_j))^2)/ ...
                    (norm(Hd(:,K+idx_j)'*W_current_mat(:,OuterIDXinterference))^2 + sum((abs(G_lk(:,idx_k)).*p_current).^2) + sigma^2);
                
    OuterUserRate = log(1+min(Rate2At1, Rate2At2));
    
%     if (OuterUserRate>=Rate_Threshold)
        Weight = InnerUserRate + OuterUserRate;
%     else
%         Weight = -1000;
%     end
%     
% else
%     Weight = -1000;
% end



end

