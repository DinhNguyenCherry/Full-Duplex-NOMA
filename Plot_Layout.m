% % % The MATLAB CODE is used for the following paper: Hieu V. Nguyen, Van-Dinh Nguyen, Octavia A. Dobre, Diep N. Nguyen, Eryk Dutkiewicz, and Oh-Soon Shin, 
% % % "Joint Power Control and User Association for NOMA-Based Full-Duplex Systems,"
% % % IEEE Transactions on Communications, vol. 67, no. 11, pp. 8037-8055, Nov. 2019.

function [ done ] = Plot_Layout( RadiusOfCell, RadiusOfInnerzone, positionUplinkUsers, positionDownlinkUsers )
%PLOT_LAYOUT Summary of this function goes here
%   Detailed explanation goes here

hold on

ang=0:0.01:2*pi; 

xp=RadiusOfCell*cos(ang);
yp=RadiusOfCell*sin(ang);
ax = plot(0+xp,0+yp);

xp1=RadiusOfInnerzone*cos(ang);
yp1=RadiusOfInnerzone*sin(ang);
ax1 = plot(0+xp1,0+yp1);

scatter(0,0 , 'k^ ')
scatter(positionUplinkUsers(:,1), positionUplinkUsers(:,2), 'b+ ')
scatter(positionDownlinkUsers(:,1), positionDownlinkUsers(:,2), 'r* ')


set(gca,'XTick',[-100 :20: 100])
set(gca,'YTick',[-100 :20: 100])

xlim([-(RadiusOfCell+5) (RadiusOfCell+5)])
ylim([-(RadiusOfCell+5) (RadiusOfCell+5)])
axis('square')

done = 1;


end

