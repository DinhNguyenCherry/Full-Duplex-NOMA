% close all
% clear
% clc

addpath(genpath('./yalmip'));
addpath(genpath('./SeDuMi_1_3'));
addpath(genpath('./SDPT3-4.0'));

tic

%% Initialization

NumberOfRunning = 5;

global MaxIteration
MaxIteration = 20;

global G
global K
global L

G = 1;
Range_G = 1:1:3;
K = 4;
L = 4;

global Rate_Threshold
Rate_Threshold = 1;

Pbs_dB = 20; % dowlink power [dB]
Pl_dB = 10 * ones(L,1);%[10; 10]; % uplink power - Pl [dB]

if (length(Pl_dB)~=L)
	disp('Error: length of uplink power vector and L are not matched together.')
end

Pbs = 10^(Pbs_dB/10); % downlink power
P = 10.^(Pl_dB/10); % uplink power

global N_tx
N_tx = 4;
N_rx = 4;

rho_dB = -30;

global rho
rho = 10^(rho_dB/10);

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

i_G = 0;
Alg1_Rate = zeros(1,length(Range_G));
Alg2_Rate = zeros(1,length(Range_G));

for g = Range_G
    
    G = g;
    i_G = i_G + 1;
    
    str_Group = ['Executing for G = ' num2str(g)];
    
    disp(['-------------------' str_Group '------------------------']);
    
    Alg1_OptValue_Stat = [];
    Alg2_OptValue_Stat = [];
    
    h = waitbar(0, [ str_Group ':  ' num2str(floor(0*100 / NumberOfRunning)) ' % Completed'],'Name','Percent of completion');

for i_NumOfSim = 1:1:NumberOfRunning
    
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




    
%     x_current = 1/sqrt(G)*ones(K,G);
%     y_current = 1/sqrt(G)*ones(L,G);
%     tau_current = G*ones(1,G);% sdpvar(1, G);

    disp('------------------------------ Initialization ---------------------------------');
    

    % Solve w0


    [ W_current ] = Get_InitialW0( H, Pbs );
%     [ W_0, ~ ] = Get_InitialW0( H, tau_current, Pbs );

%     traceWW = trace(W_0*W_0')


    % Solve p0 

    [ p_current ] = Get_Intitialp0( P );
%     [ p_0,  ~] = Get_Intitialp0( P, tau_current );



    % Intialize phi0

    phi_current = Get_Initialphi0( H, G_hat, W_current, p_current);
%     phi_0 = Get_Initialphi0( H, G_hat, W_0, p_0);

    [isBreak1, Alg1_OptValue, Alg1_OptValueChain] = Algorithm2(Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, 1);
    [isBreak2, Alg2_OptValue, Alg2_OptValueChain] = Algorithm2(Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, 0);
    
    Alg1_OptValue_Stat = [Alg1_OptValue_Stat Alg1_OptValue];
    Alg2_OptValue_Stat = [Alg2_OptValue_Stat Alg2_OptValue];
    
    if (isBreak1)
        NumOfSamples1 = NumOfSamples1 - 1;
    end
    
    if (isBreak2)
        NumOfSamples2 = NumOfSamples2 - 1;
    end
    
    Alg1_OptValueChain;
    Alg2_OptValueChain;
    
    waitbar(i_NumOfSim / NumberOfRunning, h,[ str_Group ':  ' num2str(floor(i_NumOfSim*100 / NumberOfRunning)) ' % Completed'])
    
end
    
    close(h);
    
    Alg1_Rate(i_G) = mean(Alg1_OptValue_Stat);
    Alg2_Rate(i_G) = mean(Alg2_OptValue_Stat);

end

figure(1)

hold on
plot(Range_G, Alg1_Rate, 'b.--', 'linewidth', 1, 'markersize',7);
plot(Range_G, Alg2_Rate, 'r.-', 'linewidth', 1, 'markersize',7);

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