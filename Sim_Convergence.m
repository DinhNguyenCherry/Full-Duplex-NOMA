close all
clear
clc

% pathCell = regexp(path, pathsep, 'split');
% if (~any(strcmpi('.\yalmip', pathCell)))
%     disp('_________________ Adding path ... _________________')
%     addpath(genpath('./yalmip'));
%     addpath(genpath('./SeDuMi_1_3'));
%     addpath(genpath('./SDPT3-4.0'));
% end
% clear pathCell

addpath(genpath('./yalmip'));
addpath(genpath('./SeDuMi_1_3'));
addpath(genpath('./SDPT3-4.0'));



pathname = fileparts('./Figures/Convergence/');
addpath(genpath('./Figures/Convergence'));

tic

%% Initialization

NumberOfRunning = 1;

Method = 5; % 1: Exshautive search
            % 2: Reinforcement and Bipartite graph
            % 3: Continuous relaxation problem (CRP)
            % 4: CRP - penalty functions
            % 5: FD - NOMA
            % 6: Conventional FD - no NOMA
            % 7: Half-Duplex
            %[]: an assigned mode for antennas - size 1x2

Methodname = {'ExSearch', 'HS_Search', 'CRP', 'CRP_PF', 'Con_FD', 'FD_noNOMA', 'HD'};       

% global MaxIteration
% MaxIteration = 50;

N = 10; % number of antennas at BS Nt = Nr = N
G = 1;
K = 6;
L = 4; 

% error_lk = zeros(L,K);

% error = delta*(SNR)^(-upsilon)
% delta = 0;
% upsilon = 0.7;



Rate_Threshold = 1;

max_min = 0; % 0: sum rate
             % 1: max-min with eta=1
             % 2: max-min with optional values of eta
eta = 1;

filename = ['Convergence_' Methodname{Method} 'MultiZones_test2.mat'];

loadingfile = [];%['Convergence_CRPsel1.mat'];%['Convergence.mat'];'Convergence_CRPsel.mat''Convergence_CRP_PFmultizones3.mat'
recalcD = 0;
RenewChannel = 1;

% load('Method0_Topo15_RndChannels1.mat','D_Hd');
% load('Method0_Topo15_RndChannels1.mat','D_Hu');
% load('Method0_Topo15_RndChannels1.mat','D_glk');
% 
% load('Method0_Topo15_RndChannels1.mat','positionDownlinkUsers');
% load('Method0_Topo15_RndChannels1.mat','positionUplinkUsers');

savedname = fullfile(pathname, filename);

Pbs_dB = 38; % dowlink power [dBm]
Pl_dB = 18 * ones(L,1);%[10; 10]; % uplink power - Pl [dBm]

if (length(Pl_dB)~=L)
	disp('Error: length of uplink power vector and L are not matched together.')
end

Pbs = 10^(Pbs_dB/10); % downlink power
P = 10.^(Pl_dB/10); % uplink power


rho_dB = -100;%[-90:10:-10];
full_range_rho = [-90:10:-10];

% rho = 10^(rho_dB/10);

sigma = sqrt(10^(-174/10)*10^(7));%0.01;

sigma_K = sqrt(10^(-174/10)*10^(7))*ones(1, K);%0.01*ones(1, K);

% Error_RD.error_CCI  = (diag(delta*(P/sigma^2).^(-upsilon))*ones(L,K)).^(0.5)
% Error_RD.error_up   = (diag(delta*(P/sigma^2).^(-upsilon))*ones(L,1)).^(0.5)
% Error_RD.error_down = (diag(delta*(Pbs./(sigma_K.^2)).^(-upsilon))*ones(K,1)).^(0.5)

% Error_RD.error_CCI  = 0*ones(L,K)
% Error_RD.error_up   = 0.01
% Error_RD.error_down = 0.01

RadiusOfCell = 100;
RadiusOfNearestUser = 10;
RadiusOfInnerzone = 50;
% StandardDeviation = 8;
% ploss = 3;

% Parameters = [RadiusOfCell RadiusOfNearestUser StandardDeviation ploss];
Parameters = [RadiusOfCell RadiusOfNearestUser RadiusOfInnerzone];

CenterCell = [0 0];

% FirstTier = 2*RadiusOfCell*[1 0; cos(pi/3) sin(pi/3); -cos(pi/3) sin(pi/3); -1 0; -cos(pi/3) -sin(pi/3); cos(pi/3) -sin(pi/3)];

Order = [1];

AllCells = [CenterCell];

NumOfSamples1 = NumberOfRunning;
NumOfSamples2 = NumberOfRunning;

% i_G = 0;
count = 0;
i_NumOfSim = 0;
Alg1_Rate_all = zeros(NumberOfRunning,length(full_range_rho));
Alg2_Rate_all = zeros(NumberOfRunning,length(full_range_rho));
Alg3_Rate_all = zeros(NumberOfRunning,length(full_range_rho));
HD_Rate = zeros(NumberOfRunning,length(full_range_rho));

if (isempty(loadingfile))
    [D_H1k, D_H2j, D_Hu, D_gCCI, positionDownlinkUsers, positionUplinkUsers] = CreateD( K, L, Parameters, AllCells, Order );
    D_Hd = diag([diag(D_H1k); diag(D_H2j)]);
else
    % load('Parameter_convergence1.mat','D_Hd');
    % load('Parameter_convergence1.mat','D_Hu');
    % load('Parameter_convergence1.mat','D_glk');
    % load('Parameter_convergence1.mat','positionDownlinkUsers');
    % load('Parameter_convergence1.mat','positionUplinkUsers');

%     load('Convergence.mat','positionDownlinkUsers');
%     load('Convergence.mat','positionUplinkUsers');
%     load('Convergence.mat','D_Hd');
%     load('Convergence.mat','D_Hu');
%     load('Convergence.mat','D_glk');

    load(loadingfile,'positionDownlinkUsers');
    load(loadingfile,'positionUplinkUsers');
    if (recalcD)
        [D_H1k, D_H2j, D_Hu, D_gCCI, ~, ~] = CreateD( K, L, Parameters, AllCells, Order, positionDownlinkUsers, positionUplinkUsers);
        D_Hd = diag([diag(D_H1k); diag(D_H2j)]);
    else
%         load(loadingfile,'D_H1k');
%         load(loadingfile,'D_H2j');
        load(loadingfile,'D_Hu');
        load(loadingfile,'D_Hd');
        load(loadingfile,'D_gCCI');
    end
end

Layout = Plot_Layout(RadiusOfCell, RadiusOfInnerzone, positionUplinkUsers, positionDownlinkUsers )

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
h = waitbar(0, 'Please wait ...','Name','Percent of completion');

% Clus = parcluster('local');
% Clus.NumWorkers = 4;

% poolobj = parpool(Clus, Clus.NumWorkers);



while (i_NumOfSim<NumberOfRunning)
    
    i_NumOfSim = i_NumOfSim + 1
    
    if (RenewChannel)
    
        Hu = 1*CreateChannel(1, 1, N, L)*D_Hu;
        Hd = 10*CreateChannel(1, 1, N, 2*K)*D_Hd;
%         Hd1 = 10*CreateChannel(1, 1, N, K)*D_H1k;
%         Hd2 = 10*CreateChannel(1, 1, N, K)*D_H2j;
%         D_Hd = diag([diag(D_H1k);diag(D_H2j)]);
%         Hd = [Hd1, Hd2];

        G_SI = 1*CreateChannel(1, 1, N, N);
        G_lk = 10*CreateChannel(1, 1, L, 2*K).*D_gCCI;
        
%         Error_RD.error_up   = sqrt(Error_RD.error_up)*diag(D_Hu);
%         Error_RD.error_down = sqrt(Error_RD.error_down)*diag(D_Hd);
%         Error_RD.error_CCI  = sqrt(Error_RD.error_CCI).*D_glk;
    
    else
    
%     load('Parameter_convergence1.mat','Hd');
%     load('Parameter_convergence1.mat','Hu');
%     load('Parameter_convergence1.mat','G_SI');
%     load('Parameter_convergence1.mat','G_lk');
    
%     load('Convergence.mat','Hd');
%     load('Convergence.mat','Hu');
%     load('Convergence.mat','G_SI');
%     load('Convergence.mat','G_lk');

        load(loadingfile,'Hd1');
        load(loadingfile,'Hd2');
        load(loadingfile,'Hd');
        load(loadingfile,'Hu');
        load(loadingfile,'G_SI');
        load(loadingfile,'G_lk');
    
    end
    
    count = 0;
    
    percent = (i_NumOfSim-1) / (NumberOfRunning);
    waitbar(percent, h,[ num2str(floor(percent*100)) ' % Completed'])

% for g = Range_G
for i = length(rho_dB) % for HalfDuplex
% parfor i = 1:1:length(rho_dB)
    
    i_rho = rho_dB(i);
    

    rho = 10^(i_rho/10);
    
    iRate_Threshold = Rate_Threshold*log(2);
    
 
    
    %% Full Duplex
    

    disp('------------------------------ Full Duplex Initialization ---------------------------------');
    

    if (Method==1)
        [isBreak1, Alg1_OptValue, Alg1_OptValueChain, global_ULRate_PerUser, global_DLRate_PerUser, global_Solution ] = ... 
        ProposedAlg_ExSearch(Pbs, P, Hu, Hd, G_SI, G_lk, rho, iRate_Threshold, max_min, eta, ['Convergence']);
    elseif (Method==2)
        [isBreak1, Alg1_OptValue, Alg1_OptValueChain, global_ULRate_PerUser, global_DLRate_PerUser, global_Solution ] = ... 
        ProposedAlg_Reinforcement(Pbs, P, Hu, Hd, G_SI, G_lk, rho, iRate_Threshold, max_min, eta, ['Convergence']);
    elseif (Method==3)
        [isBreak1, Alg1_OptValue, Alg1_OptValueChain, global_ULRate_PerUser, global_DLRate_PerUser, global_Solution ] = ... 
        ProposedAlg_RelaxProb(Pbs, P, Hu, Hd, G_SI, G_lk, rho, iRate_Threshold, max_min, eta, ['Convergence']);
    elseif (Method==4)
        [isBreak1, Alg1_OptValue, Alg1_OptValueChain, global_ULRate_PerUser, global_DLRate_PerUser, global_Solution ] = ... 
            ProposedAlg_RelaxProb_pen(Pbs, P, Hu, Hd, G_SI, G_lk, rho, iRate_Threshold, max_min, eta, ['Convergence']);
    else
        [isBreak1, Alg1_OptValue, Alg1_OptValueChain, global_ULRate_PerUser, global_DLRate_PerUser, global_Solution ] = ... 
        ProposedAlg_ExSearch(Pbs, P, Hu, Hd, G_SI, G_lk, rho, iRate_Threshold, max_min, eta, ['Convergence -' Methodname{Method}]);
    end
    
    
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
    
    Alg1_Rate_all(i_NumOfSim, i) = Alg1_OptValue;
    Alg2_Rate_all(i_NumOfSim, i) = Alg2_OptValue;%sum(sum(UplinkRate_PerGroupPerUser))+sum(sum(DownlinkRate_PerGroupPerUser));%
    Alg3_Rate_all(i_NumOfSim, i) = Alg3_OptValue;

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

% delete(poolobj);

% clear Clus
disp('------------------------ 2 problems ------------------------');

UplinkRate_PerGroupPerUser1
DownlinkRate_PerGroupPerUser1

UplinkRate_PerGroupPerUser2
DownlinkRate_PerGroupPerUser2

UplinkRate_PerGroupPerUser3
DownlinkRate_PerGroupPerUser3

disp('------------------------------------------------------------');

if (NumberOfRunning == 1)
    Alg1_Rate = Alg1_Rate_all
    Alg2_Rate = Alg2_Rate_all
    Alg3_Rate = Alg3_Rate_all
else
    Alg1_Rate = mean(Alg1_Rate_all)
    Alg2_Rate = mean(Alg2_Rate_all)
    Alg3_Rate = mean(Alg3_Rate_all)
end

% HD_Rate = repmat(mean(OptValue_HalfDuplex_Stat), 1, length(rho_dB));

close(h);

save(savedname);

%%
figure;

hold on

plot([1:1:length(Alg1_OptValueChain)], Alg1_OptValueChain, 'r-', 'linewidth', 2, 'markersize',9);


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