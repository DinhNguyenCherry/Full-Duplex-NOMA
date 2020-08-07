function [ D_Hd, D_Hu, D_gCCI, positionDownlinkUsers, positionUplinkUsers ] = CreateD_MultiZones( K, L, Parameters, AllCells, Order, posDLUs, posULUs )
%CREATED Summary of this function goes here
%   Detailed explanation goes here

RadiusOfCell = Parameters{1};
RadiusOfNearestUser = Parameters{2};
InnerZoneRadius = [RadiusOfNearestUser Parameters{3}];

if (nargin<nargin('CreateD_MultiZones'))
    
%% large-fading for K downlink users  - multizones

positionDownlinkUsers = [];

for iZ = 2:1:length(InnerZoneRadius)
    RadiusOfNearestUser = InnerZoneRadius(iZ-1);
    
    rvector = RadiusOfNearestUser*ones(1,K(iZ-1)) + (InnerZoneRadius(iZ)-RadiusOfNearestUser)*rand(1,K(iZ-1));
    anglevector = 2*pi*rand(1,K(iZ-1));
    
    positionDownlinkUsers = [positionDownlinkUsers; (rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];
    
end

distanceBS_To_DownlinkUsers = sqrt(sum(positionDownlinkUsers'.^2));
        
PL_downlink = GetPathloss(103.8, 20.9, distanceBS_To_DownlinkUsers/1000);

D_Hd = (diag(PL_downlink)).^(0.5)

% %% large-fading for K downlink users  - inner zone  
%     rvector = RadiusOfNearestUser*ones(1,K) + (InnerZoneRadius-RadiusOfNearestUser)*rand(1,K);
%     anglevector = 2*pi*rand(1,K);
%     
%     positionDownlinkUsers = [(rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];
%         
%     distanceBS_To_DownlinkUsers_Zone1 = sqrt(sum(positionDownlinkUsers'.^2));
%         
% 
%     
%     PL_downlink_Zone1 = GetPathloss(103.8, 20.9, distanceBS_To_DownlinkUsers_Zone1/1000);
% 
%         
%     D_H1k = (diag(PL_downlink_Zone1)).^(0.5)
%     
%     
% %% large-fading for K downlink users  - outer zone  
%     rvector = InnerZoneRadius*ones(1,K) + (RadiusOfCell-InnerZoneRadius)*rand(1,K);
%     anglevector = 2*pi*rand(1,K);
%     
%     positionDownlinkUsers = [positionDownlinkUsers; (rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];
% 
%         
%     distanceBS_To_DownlinkUsers_Zone2 = sqrt(sum(positionDownlinkUsers(K+1:2*K,:)'.^2));
%         
% 
%     PL_downlink_Zone2 = GetPathloss(103.8, 20.9, distanceBS_To_DownlinkUsers_Zone2/1000);
% 
%         
%     D_H2j = (diag(PL_downlink_Zone2)).^(0.5)    
    
%% large-fading for L uplink users  
    RadiusOfNearestUser=Parameters{2};
    rvector = RadiusOfNearestUser*ones(1,L) + (RadiusOfCell-RadiusOfNearestUser)*rand(1,L);
    anglevector = 2*pi*rand(1,L);
    
    positionUplinkUsers = [(rvector.*cos(anglevector))' (rvector.*sin(anglevector))'];

        
    distanceBS_To_UplinkUsers = sqrt(sum(positionUplinkUsers'.^2));


    PL_uplink = GetPathloss(103.8, 20.9, distanceBS_To_UplinkUsers/1000);

        
    D_Hu = (diag(PL_uplink)).^(0.5)
    
    
%% large-fading between uplink and downlink users
    
    positionDownUsers_Ltimes = kron(positionDownlinkUsers(1:sum(K),:),ones(L,1));
    positionUpUsers_2Ktimes = repmat(positionUplinkUsers(1:L,:), sum(K), 1);
    
    distance_up2KLdown = sqrt(sum((positionDownUsers_Ltimes-positionUpUsers_2Ktimes).^2 ,2));

    PL_PQup2KLdown = GetPathloss(145.4, 37.5, distance_up2KLdown/1000);
    
    D_gCCI = (reshape(PL_PQup2KLdown, L, sum(K))).^(0.5)
        

else
    
    positionDownlinkUsers = posDLUs;
    positionUplinkUsers = posULUs;
    
    disBS2DLUs = sqrt(sum(posDLUs'.^2));
    disBS2ULUs = sqrt(sum(posULUs'.^2));
    
    posDLUs_Ltimes = kron(positionDownlinkUsers,ones(L,1));
    posULUs_Ktimes = repmat(positionUplinkUsers, K, 1);
    disUp2Down = sqrt(sum((posDLUs_Ltimes-posULUs_Ktimes).^2 ,2));
    
    PL_DLUs = GetPathloss(103.8, 20.9, disBS2DLUs/1000);
    D_Hd = (PL_DLUs).^(0.5);
    D_H1k = diag(D_Hd(1:K))
    D_H2j = diag(D_Hd(K+1:2*K))
    
    PL_ULUs = GetPathloss(103.8, 20.9, disBS2ULUs/1000);  
    D_Hu = (diag(PL_ULUs)).^(0.5)
    
    PL_Up2Down = GetPathloss(145.4, 37.5, disUp2Down/1000);
    D_gCCI = (reshape(PL_Up2Down, L, K)).^(0.5)
    
end

end


function [Pathloss] = GetPathloss(firstPar, secondPar, distance)

% distance (km)

Pathloss_dB = firstPar + secondPar*log10(distance);

Pathloss = 10.^(-Pathloss_dB/10);

end

