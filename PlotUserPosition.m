close all
clear
clc

load('ForTable.mat','positionDownlinkUsers');
load('ForTable.mat','positionUplinkUsers');
load('ForTable.mat','RadiusOfCell');

% load('Rate_all_4.mat','D_H');
% load('Rate_all_4.mat','D_G_channel');
% 
% load('Rate_all_4.mat','StandardDeviation');
% load('Rate_all_4.mat','ploss');
% load('Rate_all_4.mat','RadiusOfCell');
% load('Rate_all_4.mat','RadiusOfNearestUser');


% K=size(D_H,1);
% distanceBS_To_DownlinkUsers = 1000*ones(1,K);
% 
% while (max(distanceBS_To_DownlinkUsers)>100 || min(distanceBS_To_DownlinkUsers)<10)
% 
%     Zvector = StandardDeviation*randn(1,K)
%     betavector_downlink = (diag(D_H)').^2
% 
%     distanceBS_To_DownlinkUsers = ((10.^(Zvector/10))./betavector_downlink).^(1/ploss)*RadiusOfNearestUser;
% 
% end

% anglevector = 2*pi*rand(1,K);
% positionDownlinkUsers = [(distanceBS_To_DownlinkUsers.*cos(anglevector))' (distanceBS_To_DownlinkUsers.*sin(anglevector))']


% L=size(D_G_channel,1);
% distanceBS_To_UplinkUsers = 1000*ones(1,L);
% 
% while (max(distanceBS_To_UplinkUsers)>100 || min(distanceBS_To_UplinkUsers)<10)
% 
%     Zvector = StandardDeviation*randn(1,L);
%     betavector_uplink = (diag(D_G_channel)').^2;
% 
%     distanceBS_To_UplinkUsers = ((10.^(Zvector/10))./betavector_uplink).^(1/ploss)*RadiusOfNearestUser;
% 
% end
% 
% anglevector = 2*pi*rand(1,L);
% positionUplinkUsers = [(distanceBS_To_UplinkUsers.*cos(anglevector))' (distanceBS_To_UplinkUsers.*sin(anglevector))']


%% plot

ang=0:0.01:2*pi; 
xp=RadiusOfCell*cos(ang);
yp=RadiusOfCell*sin(ang);
ax = plot(0+xp,0+yp);

hold on
scatter(0,0 , 'k^ ')

Users = size(positionDownlinkUsers,1);
for iUser = 1:1:Users
    scatter(positionDownlinkUsers(iUser,1), positionDownlinkUsers(iUser,2), 'r* ')
    text(positionDownlinkUsers(iUser,1)+4, positionDownlinkUsers(iUser,2), num2str(iUser),'color','r');
    scatter(positionUplinkUsers(iUser,1), positionUplinkUsers(iUser,2), 'b+ ')
    length(num2str(iUser))
    text(positionUplinkUsers(iUser,1)-5-3*length(num2str(iUser)), positionUplinkUsers(iUser,2), num2str(iUser),'color','b');
end

set(gca,'XTick',[-100 :20: 100])
set(gca,'YTick',[-100 :20: 100])

xlim([-105 105])
ylim([-105 105])
axis('square')

