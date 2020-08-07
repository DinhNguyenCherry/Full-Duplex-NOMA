close all
clear
clc

addpath(genpath('./yalmip'));
addpath(genpath('./SeDuMi_1_3'));
addpath(genpath('./SDPT3-4.0'));

% addpath(genpath('./CDF'));
addpath(genpath('./CDF/SumRate'));

%get the directory of your output files:
% pathname = fileparts('./CDF/');
pathname = fileparts('./CDF/SumRate/');


tic

%% Initialization

% NumberOfRunning = 40;

Method = 1; % 0: proposed
            % 1: proposed switches to half-duplex - fixed time interval
            % 2: IEEE Access - time grouping & non HA selection
            % 3: Conventional FD method - one group & non HA selection
            %[]: an assigned mode for antennas - size 1x2
            
Reused_Channels = 0; % Reused channels is for Method 1 3
Reused_from = 0;

Order_Topos = [51:1:60];
Order_RndChannels = [1:1:5];%[1:1:3];



h = waitbar(0, 'Please wait ...','Name','Percent of completion');
            
for Order_Topo = Order_Topos          
% Order_RndChannel = 3; % First group of channels for a topology of users' positions, 1


for Order_RndChannel=Order_RndChannels
    
%     if (mod(Order_RndChannel,3)==1)
%         NumberOfRunning = 20;
%     else
%         NumberOfRunning = 40;
%     end

    NumberOfRunning = 20;




rootname = ['Method' num2str(Method) '_Topo' num2str(Order_Topo) '_RndChannels'];


filename = [ rootname num2str(Order_RndChannel) '.mat'];


savedname = fullfile(pathname, filename)

dataname = ['Method' num2str(Reused_from) '_Topo' num2str(Order_Topo) '_RndChannels'];


% global MaxIteration
% MaxIteration = 50;

N = 8; % number of antennas at BS 2N
G = 2;
K = 8;
L = 8;

error_lk = zeros(L,K);
Rate_Threshold = 0.5;

eta = 1;

Pbs_dB = 26; % dowlink power [dBm]
Pl_dB = 23 * ones(L,1);%[10; 10]; % uplink power - Pl [dBm]

if (length(Pl_dB)~=L)
	disp('Error: length of uplink power vector and L are not matched together.')
end

Pbs = 10^(Pbs_dB/10); % downlink power
P = 10.^(Pl_dB/10); % uplink power


rho_dB = -70;%[-90:10:-10];
full_range_rho = [-90:10:-10];

% rho = 10^(rho_dB/10);

sigma = 0.01;

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
Alg1_Rate_all = zeros(NumberOfRunning,length(full_range_rho));
Alg2_Rate_all = zeros(NumberOfRunning,length(full_range_rho));
Alg3_Rate_all = zeros(NumberOfRunning,length(full_range_rho));
HD_Rate = zeros(NumberOfRunning,length(full_range_rho));

if (Reused_Channels==0)
    
    if (Order_RndChannel==1)

        [D_Hd, D_Hu, D_glk, positionDownlinkUsers, positionUplinkUsers] = CreateD( K, L, Parameters, AllCells, Order );
        loadingfile = [];
%         Layout = Plot_Layout(RadiusOfCell, positionUplinkUsers, positionDownlinkUsers )

    else

        loadingfile = [ rootname num2str(1) '.mat'];

        load(loadingfile,'D_Hd');
        load(loadingfile,'D_Hu');
        load(loadingfile,'D_glk');
        load(loadingfile,'positionDownlinkUsers');
        load(loadingfile,'positionUplinkUsers');

    end
    
else
    
    loadingfile = [ dataname num2str(Order_RndChannel) '.mat'];

end

% load('Parameter_convergence1.mat','D_Hd');
% load('Parameter_convergence1.mat','D_Hu');
% load('Parameter_convergence1.mat','D_glk');
% load('Parameter_convergence1.mat','positionDownlinkUsers');
% load('Parameter_convergence1.mat','positionUplinkUsers');

% load('Convergence.mat','positionDownlinkUsers');
% load('Convergence.mat','positionUplinkUsers');
% load('Convergence.mat','D_Hd');
% load('Convergence.mat','D_Hu');
% load('Convergence.mat','D_glk');



Alg1_OptValue_Stat = [];
Alg2_OptValue_Stat = [];
OptValue_HalfDuplex_Stat = [];

UplinkRate_PerGroupPerUser1 = zeros(L,G);
DownlinkRate_PerGroupPerUser1 = zeros(K,G);

UplinkRate_PerGroupPerUser2 = zeros(L,G);
DownlinkRate_PerGroupPerUser2 = zeros(K,G);

UplinkRate_PerGroupPerUser3 = zeros(L,G);
DownlinkRate_PerGroupPerUser3 = zeros(K,G);

% h = waitbar(0, [ str_Group ':  ' num2str(floor(0*100 / NumberOfRunning)) ' % Completed'],'Name','Percent of completion');
% h = waitbar(0, 'Please wait ...','Name','Percent of completion');

    percent = ((Order_Topo-Order_Topos(1))*length(Order_RndChannels)+Order_RndChannel-Order_RndChannels(1)) / (length(Order_Topos)*length(Order_RndChannels));
    waitbar(percent, h,['Topo->' num2str(Order_Topo) ', Channel file->' num2str(Order_RndChannel) ' ## Completed: ' num2str(floor(percent*100)) ' %'])

Clus = parcluster('local');
Clus.NumWorkers = 4;

poolobj = parpool(Clus, Clus.NumWorkers);

AllChannels = cell(4,NumberOfRunning);
AllChannels_Hu = cell(1,NumberOfRunning);
AllChannels_Hd = cell(1,NumberOfRunning);
AllChannels_GSI = cell(1,NumberOfRunning);
AllChannels_Glk = cell(1,NumberOfRunning);
AllOmega = cell(1,NumberOfRunning);

Method_Proposed = zeros(1, NumberOfRunning);
Method_HD = zeros(1, NumberOfRunning);
Method_nonHAsel = zeros(1, NumberOfRunning);
Method_ConvenFD = zeros(1, NumberOfRunning);

% load('Method0_RandChannels1.mat','AllChannels');
if (Reused_Channels)
    load(loadingfile,'NumberOfRunning');
    load(loadingfile,'D_Hd');
    load(loadingfile,'D_Hu');
    load(loadingfile,'D_glk');
    load(loadingfile,'AllChannels_Hu');
    load(loadingfile,'AllChannels_Hd');
    load(loadingfile,'AllChannels_GSI');
    load(loadingfile,'AllChannels_Glk');
end

% while (i_NumOfSim<NumberOfRunning)
%     
%     i_NumOfSim = i_NumOfSim + 1

Alg1_OptValue = zeros(1,NumberOfRunning);
nonHAsel = zeros(1,NumberOfRunning);

if (Reused_Channels==0)
    
    for i_NumOfSim=1:1:NumberOfRunning
    
        AllChannels_Hu{1, i_NumOfSim} = CreateChannel(1, 1, 2*N, L)*D_Hu;
        AllChannels_Hd{1, i_NumOfSim} = CreateChannel(1, 1, 2*N, K)*D_Hd;
        AllChannels_GSI{1, i_NumOfSim} = CreateChannel(1, 1, N, N);
        AllChannels_Glk{1, i_NumOfSim} = CreateChannel(1, 1, L, K).*D_glk;
        
%         AllChannels_Hu{1, i_NumOfSim} = Hu(:,:,i_NumOfSim);
%         AllChannels_Hd{1, i_NumOfSim} = Hd(:,:,i_NumOfSim);
%         AllChannels_GSI{1, i_NumOfSim} = G_SI(:,:,i_NumOfSim);
%         AllChannels_Glk{1, i_NumOfSim} = G_lk(:,:,i_NumOfSim);
    end
end

parfor i_NumOfSim=1:1:NumberOfRunning
    
    disp(['##############  Channel Index: ' num2str(i_NumOfSim) '   ##############']);
    

%     if (Reused_Channels==0)
        
        
        Hu(:,:,i_NumOfSim) = AllChannels_Hu{1, i_NumOfSim};
        Hd(:,:,i_NumOfSim) = AllChannels_Hd{1, i_NumOfSim};
        G_SI(:,:,i_NumOfSim) = AllChannels_GSI{1, i_NumOfSim};
        G_lk(:,:,i_NumOfSim) = AllChannels_Glk{1, i_NumOfSim};
    
%         Hu(:,:,i_NumOfSim) = CreateChannel(1, 1, 2*N, L)*D_Hu;
%         Hd(:,:,i_NumOfSim) = CreateChannel(1, 1, 2*N, K)*D_Hd;
%         G_SI(:,:,i_NumOfSim) = CreateChannel(1, 1, N, N);
%         G_lk(:,:,i_NumOfSim) = CreateChannel(1, 1, L, K).*D_glk;
%         
%         AllChannels_Hu{1, i_NumOfSim} = Hu(:,:,i_NumOfSim);
%         AllChannels_Hd{1, i_NumOfSim} = Hd(:,:,i_NumOfSim);
%         AllChannels_GSI{1, i_NumOfSim} = G_SI(:,:,i_NumOfSim);
%         AllChannels_Glk{1, i_NumOfSim} = G_lk(:,:,i_NumOfSim);
        

%         AllChannels{1, i_NumOfSim} = Hu;
%         AllChannels{2, i_NumOfSim} = Hd;
%         AllChannels{3, i_NumOfSim} = G_SI;
%         AllChannels{4, i_NumOfSim} = G_lk;
        
%     else
% %         Hu = AllChannels{1, i_NumOfSim};
% %         Hd = AllChannels{2, i_NumOfSim};
% %         G_SI = AllChannels{3, i_NumOfSim};
% %         G_lk = AllChannels{4, i_NumOfSim};
%         Hu = AllChannels_Hu{1, i_NumOfSim};
%         Hd = AllChannels_Hd{1, i_NumOfSim};
%         G_SI = AllChannels_GSI{1, i_NumOfSim};
%         G_lk = AllChannels_Glk{1, i_NumOfSim};
%     end

    



    
    
%     load('Parameter_convergence1.mat','Hd');
%     load('Parameter_convergence1.mat','Hu');
%     load('Parameter_convergence1.mat','G_SI');
%     load('Parameter_convergence1.mat','G_lk');
    
%     load('Convergence.mat','Hd');
%     load('Convergence.mat','Hu');
%     load('Convergence.mat','G_SI');
%     load('Convergence.mat','G_lk');
    
    
    count = 0;
    
%     percent = (i_NumOfSim-1) / (NumberOfRunning);
%     waitbar(percent, h,[ num2str(floor(percent*100)) ' % Completed'])

% for g = Range_G
for i = length(rho_dB) % for HalfDuplex
% parfor i = 1:1:length(rho_dB)
    
    i_rho = rho_dB(i);
    

    rho = 10^(i_rho/10);
    
    iRate_Threshold = Rate_Threshold;
    
 
    
    %% Full Duplex
    

    disp('------------------------------ Full Duplex Initialization ---------------------------------');
    

    % Solve w0
%     [ W_current ] = Get_InitialW0( G, Hu, Pbs );


    % Solve p0 
%     [ p_current ] = Get_Intitialp0( P );


    % Intialize phi0

%     phi_current = Get_Initialphi0( H, G_hat, W_current, p_current);


%     [isBreak1, Alg1_OptValue, Alg1_OptValueChain, Alg1_OptValueChain_sum, UplinkRate_PerGroupPerUser1, DownlinkRate_PerGroupPerUser1, prework] = ProposedAlg(Pbs, P, Hu, Hd, G_SI, G_lk, error_lk, rho, 0, Rate_Threshold, 0);
    [isBreak1, Alg1_OptValue(i_NumOfSim), Alg1_OptValueChain, Alg1_OptValueChain_sum, UplinkRate_PerGroupPerUser1, DownlinkRate_PerGroupPerUser1, prework, Omega, AntennaMode, OptValue_permode,~,nonHAsel(i_NumOfSim),~] = ProposedAlg(Pbs, P, Hu(:,:,i_NumOfSim), Hd(:,:,i_NumOfSim), G_SI(:,:,i_NumOfSim), G_lk(:,:,i_NumOfSim), error_lk, rho, Method, iRate_Threshold, 0, eta);
%     Alg1_OptValue = 0;
%     
%     if (isBreak1)
%         NumOfSamples1 = NumOfSamples1 - 1;
%         i_NumOfSim = i_NumOfSim - 1;
%         break;
%     end
    

%     [isBreak2, Alg2_OptValue, Alg2_OptValueChain, Alg2_OptValueChain_sum, UplinkRate_PerGroupPerUser2, DownlinkRate_PerGroupPerUser2, ~] = ProposedAlg(Pbs, P, Hu, Hd, G_SI, G_lk, error_lk, rho, 0, Rate_Threshold, 1);
    Alg2_OptValue = 0;
    
%     if (isBreak2)
%         NumOfSamples2 = NumOfSamples2 - 1;
%         i_NumOfSim = i_NumOfSim - 1;
%         break;
%     end
%     
%     Alg1_Rate_all(i_NumOfSim, count) = Alg1_OptValue;
%     Alg2_Rate_all(i_NumOfSim, count) = Alg2_OptValue;

%     [isBreak3, Alg3_OptValue, Alg3_OptValueChain, Alg3_OptValueChain_sum, UplinkRate_PerGroupPerUser3, DownlinkRate_PerGroupPerUser3, ~] = ProposedAlg(Pbs, P, Hu, Hd, G_SI, G_lk, error_lk, rho, 0, Rate_Threshold, 2, eta);
    Alg3_OptValue = 0;
    
%     Alg1_Rate_all(i_NumOfSim, i) = Alg1_OptValue;
%     Alg2_Rate_all(i_NumOfSim, i) = Alg2_OptValue;%sum(sum(UplinkRate_PerGroupPerUser))+sum(sum(DownlinkRate_PerGroupPerUser));%
%     Alg3_Rate_all(i_NumOfSim, i) = Alg3_OptValue;
    
    if (Method==0)
        Method_Proposed(1, i_NumOfSim) = Alg1_OptValue(i_NumOfSim);
        AllOmega{1, i_NumOfSim} = Omega;
%         temp = nonHAsel{1,1};
        Method_nonHAsel(1, i_NumOfSim) = nonHAsel(i_NumOfSim);%temp;
    else
        if (Method==3)
            Method_ConvenFD(1, i_NumOfSim) = Alg1_OptValue(i_NumOfSim);
        else
            Method_HD(1, i_NumOfSim) = Alg1_OptValue(i_NumOfSim);
        end
    end

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
    
%     [ W_current ] = Get_InitialW0( H, Pbs, 1 );
%     
%     [ p_current ] = Get_Intitialp0( P );
%     
%     phi_current = Get_Initialphi0( H, G_hat, W_current, p_current, 1);
% 
%     [OptValue_HalfDuplex] = HalfDuplex( Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, rho, 1);
%     
%     OptValue_HalfDuplex_Stat = [OptValue_HalfDuplex_Stat OptValue_HalfDuplex];
%     
%     end
    
    
    
end
%     waitbar(i_NumOfSim / NumberOfRunning, h,[ str_Group ':  ' num2str(floor(i_NumOfSim*100 / NumberOfRunning)) ' % Completed'])
    
    

    
%     if (count==1)
%         HD_Rate(count) = mean(OptValue_HalfDuplex_Stat);
%         HD_Rate = repmat(HD_Rate(count), 1, length(rho_dB));
%     end

end

delete(poolobj);

clear Clus


% if (NumberOfRunning == 1)
%     Alg1_Rate = Alg1_Rate_all
%     Alg2_Rate = Alg2_Rate_all
%     Alg3_Rate = Alg3_Rate_all
% else
%     Alg1_Rate = mean(Alg1_Rate_all)
%     Alg2_Rate = mean(Alg2_Rate_all)
%     Alg3_Rate = mean(Alg3_Rate_all)
% end

% HD_Rate = repmat(mean(OptValue_HalfDuplex_Stat), 1, length(rho_dB));



% figure;
% 
% hold on

% plot([1:1:length(Alg1_OptValueChain)], Alg1_OptValueChain, 'r-', 'linewidth', 2, 'markersize',9);
% plot([1:1:length(nonHAsel{1,2})], nonHAsel{1,2}, 'b-', 'linewidth', 2, 'markersize',9);
% plot([1:1:length(Alg2_OptValueChain_sum)], Alg2_OptValueChain_sum, 'r-', 'linewidth', 2, 'markersize',9);
% plot([1:1:length(Alg3_OptValueChain)], Alg3_OptValueChain, 'g-', 'linewidth', 2, 'markersize',9);
% plot([1:1:length(Alg3_OptValueChain)], Alg3_OptValueChain*eta, 'g-', 'linewidth', 2, 'markersize',9);


% plot(full_range_rho, Alg1_Rate, 'bs--', 'linewidth', 2, 'markersize',9);
% plot(full_range_rho, Alg2_Rate, 'rs-', 'linewidth', 2, 'markersize',9);
% plot(rho_dB, HD_Rate, 'k^-', 'linewidth', 2, 'markersize',9);

% legend('FD - fixed time group', 'FD - optimal time group');

save(savedname, '-regexp', '^(?!h$).');

end

end

close(h);

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