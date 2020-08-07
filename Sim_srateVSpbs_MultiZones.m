close all
clear
clc

pathCell = regexp(path, pathsep, 'split');
if (~any(strcmpi('.\yalmip', pathCell)))
    disp('_________________ Adding path ... _________________')
    addpath(genpath('./yalmip'));
    addpath(genpath('./SeDuMi_1_3'));
    addpath(genpath('./SDPT3-4.0'));
end
clear pathCell

pathname = fileparts('./Figures/SrateVSpbs_MultiZones_zone3/');

tic

%% Initialization


Files = 1:2;
NumberOfRunning = 10;

Method = 4; % 1: Exshautive search
            % 2: Reinforcement and Bipartite graph
            % 3: Continuous relaxation problem (CRP)
            % 4: CRP - penalty functions
            % 5: FD - NOMA
            % 6: Conventional FD - no NOMA
            % 7: Half-Duplex
            %[]: an assigned mode for antennas - size 1x2

Methodname = {'ExSearch', 'HS_Search', 'CRP', 'CRP_PF', 'Con_FD', 'FD_noNOMA', 'HD'}; 
CellSize = '_smallcell';

NoWorkers = 4;

% global MaxIteration
% MaxIteration = 50;

N = 10; % number of antennas at BS Nt = Nr = N
G = 1;
K = 4;
L = 4; 

NoUsers = [4 2 2 4];

% error_lk = zeros(L,K);

% error = delta*(SNR)^(-upsilon)
% delta = 0;
% upsilon = 0.7;

recalcD = 1;

Rate_Threshold = 1;

max_min = 0; % 0: sum rate
             % 1: max-min with eta=1
             % 2: max-min with optional values of eta
eta = 1;

zone3 = 1;

h = waitbar(0, 'Please wait ...','Name','Percent of completion');

for iFile = Files
    
    percent = (find(iFile==Files)-1) / (length(Files));
    waitbar(percent, h,[ num2str(floor(percent*100)) ' % Completed'])

filename = ['[' Methodname{Method} CellSize '] SrateVSpbs_' num2str(iFile) '.mat']
loadingfile = [];%['Convergence.mat'];

RenewChannel = 1;

% load('Method0_Topo15_RndChannels1.mat','D_Hd');
% load('Method0_Topo15_RndChannels1.mat','D_Hu');
% load('Method0_Topo15_RndChannels1.mat','D_glk');
% 
% load('Method0_Topo15_RndChannels1.mat','positionDownlinkUsers');
% load('Method0_Topo15_RndChannels1.mat','positionUplinkUsers');

savedname = fullfile(pathname, filename);

Pbs_dB = 26:4:42; % dowlink power [dBm]
Pl_dB = 18 * ones(L,1);%[10; 10]; % uplink power - Pl [dBm]

if (length(Pl_dB)~=L)
	disp('Error: length of uplink power vector and L are not matched together.')
end

Pbs = 10.^(Pbs_dB/10); % downlink power
P = 10.^(Pl_dB/10); % uplink power


rho_dB = -100;%[-90:10:-10];
full_range_rho = [-90:10:-10];

rho = 10^(rho_dB/10);

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
RadiusOfInnerzone = [30 50 70 100];
% StandardDeviation = 8;
% ploss = 3;

% Parameters = [RadiusOfCell RadiusOfNearestUser StandardDeviation ploss];
Parameters = {RadiusOfCell, RadiusOfNearestUser, RadiusOfInnerzone};

CenterCell = [0 0];

% FirstTier = 2*RadiusOfCell*[1 0; cos(pi/3) sin(pi/3); -cos(pi/3) sin(pi/3); -1 0; -cos(pi/3) -sin(pi/3); cos(pi/3) -sin(pi/3)];

Order = [1];

AllCells = [CenterCell];

NumOfSamples1 = NumberOfRunning;
NumOfSamples2 = NumberOfRunning;

% i_G = 0;
count = 0;
% i_NumOfSim = 0;
% OptSRate_all = zeros(NumberOfRunning,length(Pbs));
OptSRate_vec = zeros(1,NumberOfRunning*length(Pbs));
Alg2_Rate_all = zeros(NumberOfRunning,length(Pbs));
Alg3_Rate_all = zeros(NumberOfRunning,length(Pbs));
HD_Rate = zeros(NumberOfRunning,length(Pbs));

if (isempty(loadingfile))
    [D_Hd, D_Hu, D_gCCI, positionDownlinkUsers, positionUplinkUsers] = CreateD_MultiZones( NoUsers, L, Parameters, AllCells, Order );
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
    else
        load(loadingfile,'D_H1k');
        load(loadingfile,'D_H2j');
        load(loadingfile,'D_Hu');
        load(loadingfile,'D_gCCI');
    end
end

% Layout = Plot_Layout(RadiusOfCell, RadiusOfInnerzone, positionUplinkUsers, positionDownlinkUsers )

% Alg1_OptValue_Stat = [];
% Alg2_OptValue_Stat = [];
% OptValue_HalfDuplex_Stat = [];

% UplinkRate_PerGroupPerUser1 = zeros(L,G);
% DownlinkRate_PerGroupPerUser1 = zeros(K,G);
% 
% UplinkRate_PerGroupPerUser2 = zeros(L,G);
% DownlinkRate_PerGroupPerUser2 = zeros(K,G);
% 
% UplinkRate_PerGroupPerUser3 = zeros(L,G);
% DownlinkRate_PerGroupPerUser3 = zeros(K,G);

% h = waitbar(0, [ str_Group ':  ' num2str(floor(0*100 / NumberOfRunning)) ' % Completed'],'Name','Percent of completion');
% h = waitbar(0, 'Please wait ...','Name','Percent of completion');

AllChannels = cell(NumberOfRunning, 4);

for i_NumOfSim = 1:1:NumberOfRunning
    if (RenewChannel)
    
        Hu = 1*CreateChannel(1, 1, N, L)*D_Hu;
        Hd = 10*CreateChannel(1, 1, N, 3*K)*D_Hd;
%         Hd1 = 10*CreateChannel(1, 1, N, K)*D_H1k;
%         Hd2 = 10*CreateChannel(1, 1, N, K)*D_H2j;
%         Hd = [Hd1, Hd2];

        G_SI = 1*CreateChannel(1, 1, N, N);
        G_lk = 10*CreateChannel(1, 1, L, 3*K).*D_gCCI;
        
        AllChannels{i_NumOfSim,1} = Hu;
        AllChannels{i_NumOfSim,2} = Hd;
        AllChannels{i_NumOfSim,3} = G_SI;
        AllChannels{i_NumOfSim,4} = G_lk;
        

    
    else
    
        
        load(loadingfile,'AllChannels');
    
    end
end


Clus = parcluster('local');
Clus.NumWorkers = NoWorkers;

poolobj = parpool(Clus, Clus.NumWorkers);

i_NumOfSim = 0;

AllPoints = NumberOfRunning*length(Pbs);

% while (i_NumOfSim<NumberOfRunning)
parfor iPoint = 0:1:(AllPoints-1)
    
%     i_NumOfSim = i_NumOfSim + 1
    
    i_NumOfSim = floor(iPoint/length(Pbs))+1;
    
    Hu = AllChannels{i_NumOfSim,1};
    Hd = AllChannels{i_NumOfSim,2};
    G_SI = AllChannels{i_NumOfSim,3};
    G_lk = AllChannels{i_NumOfSim,4};
    
    
    
%     count = 0;
    
%     percent = (i_NumOfSim-1) / (NumberOfRunning);
%     waitbar(percent, h,[ num2str(floor(percent*100)) ' % Completed'])

% Clus = parcluster('local');
% Clus.NumWorkers = 3;
% 
% poolobj = parpool(Clus, Clus.NumWorkers);
    
% for g = Range_G
% parfor i = 1:1:length(Pbs) % for HalfDuplex
% parfor i = 1:1:length(rho_dB)
    
    i = mod(iPoint,length(Pbs))+1;
    iPbs = Pbs(i);
    

%     rho = 10^(i_rho/10);
    
    iRate_Threshold = Rate_Threshold*log(2);
    
    strdisp = ['(rho, RateThres, Pbs, NoChannel) = (' num2str(rho_dB) ', ' num2str(Rate_Threshold) ', ' num2str(Pbs_dB(i)) ', ' num2str(i_NumOfSim) ') - ' Methodname{Method}] %#ok<PFBNS>
    
    %% Full Duplex
    

    disp('------------------------------ Full Duplex Initialization ---------------------------------');
    

    if (Method==1)
        [isBreak1, Alg1_OptValue, Alg1_OptValueChain, global_ULRate_PerUser, global_DLRate_PerUser, global_Solution ] = ... 
        ProposedAlg_ExSearch(iPbs, P, Hu, Hd, G_SI, G_lk, rho, iRate_Threshold, max_min, eta, strdisp);
    elseif (Method==2)
        [isBreak1, Alg1_OptValue, Alg1_OptValueChain, global_ULRate_PerUser, global_DLRate_PerUser, global_Solution ] = ... 
        ProposedAlg_Reinforcement(iPbs, P, Hu, Hd, G_SI, G_lk, rho, iRate_Threshold, max_min, eta, strdisp);
    elseif (Method==3)
        [isBreak1, Alg1_OptValue, Alg1_OptValueChain, global_ULRate_PerUser, global_DLRate_PerUser, global_Solution ] = ... 
        ProposedAlg_RelaxProb(iPbs, P, Hu, Hd, G_SI, G_lk, rho, iRate_Threshold, max_min, eta, strdisp);
    elseif (Method==4)
        if (zone3 ==1)
            [isBreak1, Alg1_OptValue, Alg1_OptValueChain, global_ULRate_PerUser, global_DLRate_PerUser, global_Solution ] = ... 
            ProposedAlg_RelaxProb_pen_MultiZones(iPbs, P, Hu, Hd, G_SI, G_lk, rho, iRate_Threshold, max_min, eta, strdisp);
        else
        [isBreak1, Alg1_OptValue, Alg1_OptValueChain, global_ULRate_PerUser, global_DLRate_PerUser, global_Solution ] = ... 
            ProposedAlg_RelaxProb_pen(iPbs, P, Hu, Hd, G_SI, G_lk, rho, iRate_Threshold, max_min, eta, strdisp);
        end
    else
        [isBreak1, Alg1_OptValue, Alg1_OptValueChain, global_ULRate_PerUser, global_DLRate_PerUser, global_Solution ] = ... 
        ProposedAlg_ExSearch(iPbs, P, Hu, Hd, G_SI, G_lk, rho, iRate_Threshold, max_min, eta, strdisp);
    end
    
    
%     Alg1_OptValue = 0;
%     
%     if (isBreak1)
%         NumOfSamples1 = NumOfSamples1 - 1;
%         i_NumOfSim = i_NumOfSim - 1;
%         break;
%     end
    

%     [isBreak2, Alg2_OptValue, Alg2_OptValueChain, Alg2_OptValueChain_sum, UplinkRate_PerGroupPerUser2, DownlinkRate_PerGroupPerUser2, ~] = ProposedAlg(Pbs, P, Hu, Hd, G_SI, G_lk, error_lk, rho, 0, Rate_Threshold, 1);
%     Alg2_OptValue = 0;
    
%     if (isBreak2)
%         NumOfSamples2 = NumOfSamples2 - 1;
%         i_NumOfSim = i_NumOfSim - 1;
%         break;
%     end
%     
%     Alg1_Rate_all(i_NumOfSim, count) = Alg1_OptValue;
%     Alg2_Rate_all(i_NumOfSim, count) = Alg2_OptValue;

%     [isBreak3, Alg3_OptValue, Alg3_OptValueChain, Alg3_OptValueChain_sum, UplinkRate_PerGroupPerUser3, DownlinkRate_PerGroupPerUser3, ~] = ProposedAlg(Pbs, P, Hu, Hd, G_SI, G_lk, error_lk, rho, 0, Rate_Threshold, 2, eta);
%     Alg3_OptValue = 0;
    
    OptSRate_vec(iPoint+1) = Alg1_OptValue;
%     OptSRate_all(i_NumOfSim, i) = Alg1_OptValue;
%     Alg2_Rate_all(i_NumOfSim, i) = Alg2_OptValue;%sum(sum(UplinkRate_PerGroupPerUser))+sum(sum(DownlinkRate_PerGroupPerUser));%
%     Alg3_Rate_all(i_NumOfSim, i) = Alg3_OptValue;


    
    
    
% end

% delete(poolobj);
% 
% clear Clus

end
    

delete(poolobj);

clear Clus

    temp = reshape(OptSRate_vec,length(Pbs),NumberOfRunning);
    OptSRate_all = temp';

% disp('------------------------ 2 problems ------------------------');
% 
% UplinkRate_PerGroupPerUser1
% DownlinkRate_PerGroupPerUser1
% 
% UplinkRate_PerGroupPerUser2
% DownlinkRate_PerGroupPerUser2
% 
% UplinkRate_PerGroupPerUser3
% DownlinkRate_PerGroupPerUser3
% 
% disp('------------------------------------------------------------');

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

% close(h);

save(savedname, '-regexp', '^(?!h$).');

end

close(h);

%%
% figure;
% 
% hold on

% plot([1:1:length(Alg1_OptValueChain)], Alg1_OptValueChain, 'r-', 'linewidth', 2, 'markersize',9);


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