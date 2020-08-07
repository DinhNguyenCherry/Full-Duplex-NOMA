close all
clear
clc

addpath(genpath('./yalmip'));
addpath(genpath('./SeDuMi_1_3'));
addpath(genpath('./SDPT3-4.0'));

tic

%% Initialization

NumberOfRunning = 40;

filename = 'SRvsThreshold_Alg2_Users10_G4_3.mat';

global MaxIteration
MaxIteration = 50;

% global G
% global K
% global L

G = 4;
Range_G = 2;
K = 10;
L = 10;

% global Rate_Threshold
% Rate_Threshold = 0.5:0.1:1.2;
Rate_Threshold = [1.5:0.5:3];%0.5:0.1:1.2;%0.5:0.25:2.25;

Pbs_dB = 26; % dowlink power [dBm]
Pl_dB = 10 * ones(L,1);%[10; 10]; % uplink power - Pl [dBm]

if (length(Pl_dB)~=L)
	disp('Error: length of uplink power vector and L are not matched together.')
end

Pbs = 10^(Pbs_dB/10); % downlink power
P = 10.^(Pl_dB/10); % uplink power

global N_tx
N_tx = 10;
global N_rx
N_rx = 10;

rho_dB = [-75];
full_range_rho = [-90:10:-10];

% global rho
% rho = 10^(rho_dB/10);

global sigma
sigma = 0.01;
global sigma_K
sigma_K = 0.01*ones(1, K);

RadiusOfCell = 100;
RadiusOfNearestUser = 10;
StandardDeviation = 8;
ploss = 3;

Parameters = [RadiusOfCell RadiusOfNearestUser StandardDeviation ploss];

CenterCell = [0 0];

% FirstTier = 2*RadiusOfCell*[1 0; cos(pi/3) sin(pi/3); -cos(pi/3) sin(pi/3); -1 0; -cos(pi/3) -sin(pi/3); cos(pi/3) -sin(pi/3)];

Order = [1];

AllCells = [CenterCell];% FirstTier; SecondTier];

NumOfSamples1 = NumberOfRunning;
NumOfSamples2 = NumberOfRunning;

% i_G = 0;
count = 0;
i_NumOfSim = 0;
Alg1_Rate_all = zeros(NumberOfRunning,length(Rate_Threshold));
Alg2_Rate_all = zeros(NumberOfRunning,length(Rate_Threshold));
HD_Rate = zeros(NumberOfRunning,length(Rate_Threshold));
HD_Rate_All = zeros(NumberOfRunning, length(Rate_Threshold));
HD_DownRate_All = zeros(NumberOfRunning, length(Rate_Threshold));
HD_UpRate_All = zeros(NumberOfRunning, length(Rate_Threshold));

DownlinkRate_PerGroupPerUser = zeros(K, G);
UplinkRate_PerGroupPerUser = zeros(L, G);
time_group = zeros(1,G);

% [D_H, D_G_channel, D_G_hat, positionDownlinkUsers, positionUplinkUsers ] = CreateD( K, L, Parameters, AllCells, Order );

%     H_full = load('ForTable3.mat', 'H')
%     G_channel_full = load('ForTable3.mat', 'G_channel')
%     G_SI_full = load('ForTable3.mat', 'G_SI')
%     G_hat_full = load('ForTable3.mat', 'G_hat')

    load('ForTable.mat', 'D_H')
    load('ForTable.mat', 'D_G_channel')
    load('ForTable.mat', 'D_G_hat')



% load('SRvsThreshold_Alg2_G3_Users10_1.mat', 'D_H')
% load('SRvsThreshold_Alg2_G3_Users10_1.mat', 'D_G_channel')
% load('SRvsThreshold_Alg2_G3_Users10_1.mat', 'D_G_hat')

% load('Rate_all_4.mat', 'D_H')
% load('Rate_all_4.mat', 'D_G_channel')
% load('Rate_all_4.mat', 'D_G_hat')




Alg1_OptValue_Stat = [];
Alg2_OptValue_Stat = [];
OptValue_HalfDuplex_Stat = [];

% h = waitbar(0, [ str_Group ':  ' num2str(floor(0*100 / NumberOfRunning)) ' % Completed'],'Name','Percent of completion');
h = waitbar(0, 'Please wait ...','Name','Percent of completion');

if (length(Rate_Threshold)>1)

Clus = parcluster('local');
Clus.NumWorkers = 4;

poolobj = parpool(Clus, Clus.NumWorkers);
end

while (i_NumOfSim<NumberOfRunning)
    
    i_NumOfSim = i_NumOfSim + 1
    
    H = CreateChannel(1, 1, N_tx, K)*D_H;
    G_channel = CreateChannel(1, 1, N_rx, L)*D_G_channel;

    G_SI = CreateChannel(1, 1, N_tx, N_rx);
    G_hat = CreateChannel(1, 1, L, K).*D_G_hat;
    
%     H = H_full.H(:, 1:K);
%     G_channel = G_channel_full.G_channel(:, 1:L);
% 
%     G_SI = G_SI_full.G_SI;
%     G_hat = G_hat_full.G_hat(1:L, 1:K);
    
% Half Duplex Channel
%     H = CreateChannel(1, 1, N_tx+N_rx, K)*D_H;
%     G_channel = CreateChannel(1, 1, N_rx+N_tx, L)*D_G_channel;
% 
%     G_SI = zeros(N_rx+N_tx,N_rx+N_tx); %CreateChannel(1, 1, N_tx, N_rx);
%     G_hat = CreateChannel(1, 1, L, K).*D_G_hat;
    
    count = 0;
    
    percent = (i_NumOfSim-1) / (NumberOfRunning);
    waitbar(percent, h,[ num2str(floor(percent*100)) ' % Completed'])

% for g = Range_G
parfor i = 1:1:length(Rate_Threshold) % for HalfDuplex
% for i = 1:1:length(Rate_Threshold) % for HalfDuplex
% parfor i = 1:1:length(rho_dB)
    
    i_Rate_Threshold = Rate_Threshold(i);
    
    str = ['Rate Threshold = ' num2str(i_Rate_Threshold)];
    disp(['-----------------' str '-----------------']);
    
%     G = g;
    rho = 10^(rho_dB/10);
%     count = count + 1;
    
%     str_Group = ['Executing for SI = ' num2str(i_rho) ' dB'];
%     
%     disp(['-------------------' str_Group '------------------------']);
    
%     percent = ((i_NumOfSim-1)*length(rho_dB)+(count-1)) / (NumberOfRunning*length(rho_dB));
%     waitbar(percent, h,[ num2str(floor(percent*100)) ' % Completed - ' str_Group])
    
 
    
    %% Full Duplex
    
%     [D_H, D_G_channel, D_G_hat] = CreateD( K, L, Parameters, AllCells, Order );
% 
%     H = CreateChannel(1, 1, N_tx, K)*D_H;
%     G_channel = CreateChannel(1, 1, N_rx, L)*D_G_channel;
% 
%     G_SI = CreateChannel(1, 1, N_tx, N_rx);
%     G_hat = CreateChannel(1, 1, L, K).*D_G_hat;
    


%     H = CreateChannel(1, 1, N_tx, K);
%     G_channel = CreateChannel(1, 1, N_rx, L);
% 
%     G_SI = CreateChannel(1, 1, N_tx, N_rx);
%     G_hat = CreateChannel(1, 1, L, K);


    disp('------------------------------ Full Duplex Initialization ---------------------------------');
    

    % Solve w0
    [ W_current ] = Get_InitialW0( G, H, Pbs );


    % Solve p0 
    [ p_current ] = Get_Intitialp0( G, P );


    % Intialize phi0

    phi_current = Get_Initialphi0( H, G_hat, W_current, p_current);

%     [isBreak1, Alg1_OptValue, Alg1_OptValueChain, DownlinkRate, UplinkRate, time_opt] = Algorithm2(Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, rho, 1, i_Rate_Threshold);
    Alg1_OptValue = 0;
%     
%     if (isBreak1)
%         NumOfSamples1 = NumOfSamples1 - 1;
%         i_NumOfSim = i_NumOfSim - 1;
%         break;
%     end
    
    [isBreak2, Alg2_OptValue, Alg2_OptValueChain, DownlinkRate, UplinkRate, time_opt] = Algorithm2(Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, rho, 0, i_Rate_Threshold);
    
    DownlinkRate_PerGroupPerUser = DownlinkRate;
    UplinkRate_PerGroupPerUser = UplinkRate;
    time_group = time_opt;
%     Alg2_OptValue = 0;
    
%     if (isBreak2)
%         NumOfSamples2 = NumOfSamples2 - 1;
%         i_NumOfSim = i_NumOfSim - 1;
%         break;
%     end
%     
%     Alg1_Rate_all(i_NumOfSim, count) = Alg1_OptValue;
%     Alg2_Rate_all(i_NumOfSim, count) = Alg2_OptValue;
    
    Alg1_Rate_all(i_NumOfSim, i) = Alg1_OptValue;
    Alg2_Rate_all(i_NumOfSim, i) = Alg2_OptValue;

%     Alg1_Rate_all(i_NumOfSim, i+5) = Alg1_OptValue;
%     Alg2_Rate_all(i_NumOfSim, i+5) = Alg2_OptValue;
    
%     Alg1_OptValue_Stat = [Alg1_OptValue_Stat Alg1_OptValue];
%     Alg2_OptValue_Stat = [Alg2_OptValue_Stat Alg2_OptValue];
    

    

%     
%     Alg1_OptValueChain;
%     Alg2_OptValueChain;

    %% Half Duplex
    
%     if (count==1)
%     
%     disp('------------------------------ Half Duplex Initialization ---------------------------------');
%     
% %     [D_H, D_G_channel, D_G_hat] = CreateD( K, L, Parameters, AllCells, Order );
% 
%     H = CreateChannel(1, 1, N_tx+N_rx, K)*D_H;
%     G_channel = CreateChannel(1, 1, N_rx+N_tx, L)*D_G_channel;
% 
%     G_SI = zeros(N_rx+N_tx,N_rx+N_tx);%CreateChannel(1, 1, N_tx, N_rx);
%     G_hat = CreateChannel(1, 1, L, K).*D_G_hat;
    
%     [ W_current ] = Get_InitialW0(G, H, Pbs, 1 );
%     
%     [ p_current ] = Get_Intitialp0( P );
%     
%     phi_current = Get_Initialphi0( H, G_hat, W_current, p_current, 1);
% 
%     [OptValue_HalfDuplex, OptimalDownlinkRate, OptimalUplinkRate] = HalfDuplex( Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, rho, 1, 2*i_Rate_Threshold);
    OptValue_HalfDuplex = 0;
    OptimalDownlinkRate = 0;
    OptimalUplinkRate = 0;
%     
    HD_Rate_All(i_NumOfSim, i) = OptValue_HalfDuplex;
    HD_DownRate_All(i_NumOfSim, i) = OptimalDownlinkRate;
    HD_UpRate_All(i_NumOfSim, i) = OptimalUplinkRate;
%     
%     end
    
    
    
end
%     waitbar(i_NumOfSim / NumberOfRunning, h,[ str_Group ':  ' num2str(floor(i_NumOfSim*100 / NumberOfRunning)) ' % Completed'])
    
    

    
%     if (count==1)
%         HD_Rate(count) = mean(OptValue_HalfDuplex_Stat);
%         HD_Rate = repmat(HD_Rate(count), 1, length(rho_dB));
%     end
    
end

if (length(Rate_Threshold)>1)
delete(poolobj);

clear Clus
end

if (NumberOfRunning == 1)
    Alg1_Rate = Alg1_Rate_all
    Alg2_Rate = Alg2_Rate_all
    HD_Rate = HD_Rate_All
else
    Alg1_Rate = mean(Alg1_Rate_all)
    Alg2_Rate = mean(Alg2_Rate_all)
    HD_Rate = mean(HD_Rate_All)
end

close(h);

save(filename);

figure(1)

hold on
% plot(Rate_Threshold, Alg1_Rate, 'bs-', 'linewidth', 2, 'markersize',9);
plot(Rate_Threshold, Alg2_Rate, 'rs-', 'linewidth', 2, 'markersize',9);
% plot(Rate_Threshold, HD_Rate, 'k^-', 'linewidth', 2, 'markersize',9);

% legend('FD - fixed time group', 'FD - optimal time group');

time = toc
hour = 0; min=0;
if (time>3600)
    hour = floor(time/3600);
    time = mod(time,3600);
end
if (time>60)
    min = floor(time/60);
    time = mod(time,60);
end
disp(['Running Time = ' num2str(hour) ' h ' num2str(min) ' m ' num2str(time) ' s.' ]);