% % % The MATLAB CODE is used for the following paper: Hieu V. Nguyen, Van-Dinh Nguyen, Octavia A. Dobre, Diep N. Nguyen, Eryk Dutkiewicz, and Oh-Soon Shin, 
% % % "Joint Power Control and User Association for NOMA-Based Full-Duplex Systems,"
% % % IEEE Transactions on Communications, vol. 67, no. 11, pp. 8037-8055, Nov. 2019.


function [ isBreak, global_OptValue, global_OptValueChain, global_ULRate_PerUser, global_DLRate_PerUser, global_Solution ] = ProposedAlg_RelaxProb( Pbs, P, Hu, Hd, G_SI, G_lk, rho, Rate_Threshold, max_min, eta, strdisp )
%PROPOSEDALG Summary of this function goes here
%   Detailed explanation goes here

if (nargin==nargin('ProposedAlg_RelaxProb')-1)
    strdisp = [];
end



if (Rate_Threshold==0)
    Rate_Threshold = 0.001;
end

if (max_min)
    Rate_Threshold = 0;
end

Converg = strcmp(strdisp,'Convergence');
if (Converg)
    MaxIteration = 50;
else
    MaxIteration = 50;
end


N = floor(size(Hu,1)/2);
L = size(Hu,2);
K = size(Hd,2);


global_OptValueChain = [];

global_Solution = 0;

if (isempty(strdisp))
    strdisp = ['For SI = ' num2str(10*log10(rho)) ' dB & Pbs = ' num2str(10*log10(Pbs)) ' dBm'];
end

global_OptValue = -1000;

   
    
        disp(['***** Getting a feasible point ..... ' strdisp]);

        GetBreak = 1;
        
        while (GetBreak)
            [ GetBreak, OptimalValue, OptimalValue_preStep, ~, ~, W_current, p_current, addVar_current] = ...
                    GetInitialization(Pbs, P, Hu, Hd, G_SI, G_lk, rho, Rate_Threshold, max_min, eta, {1, 'DL'}, {});
            RowSum = sum(round(addVar_current{1}),1);
            ColSum = sum(round(addVar_current{1}),2);
            if (~isempty(find(RowSum<1)) || ~isempty(find(ColSum<1)))
                GetBreak = 1;
            end
        end
        GetBreak = 1;
        while (GetBreak)
            [ GetBreak, OptimalValue, OptimalValue_preStep, ~, ~, W_current, p_current, addVar_current] = ...
                    GetInitialization(Pbs, P, Hu, Hd, G_SI, G_lk, rho, Rate_Threshold, max_min, eta, {1, 'UL'}, {W_current, p_current, addVar_current});
        end
        if (~GetBreak)
            GetBreak = 1;
        
            while (GetBreak)
                  [ GetBreak, OptimalValue, OptimalValue_preStep, ~, ~, W_current, p_current, addVar_current] = ...
                    GetInitialization(Pbs, P, Hu, Hd, G_SI, G_lk, rho, Rate_Threshold, max_min, eta, {2, 'DL/UL'}, {W_current, p_current, addVar_current});
            end        
        end


        %% Loop for optimization

        if (GetBreak)

            isBreak = 1;
            OptValueChain = [];
            UplinkRate_PerUser = zeros(L,1);
            DownlinkRate_PerUser = zeros(K,2);
            
            

        else

            isBreak = 0;

            disp(['---------------------------------- ' strdisp ' - Run iterative algorithm --------------------------------']);

            n = 0;
            OptimalValue_current = OptimalValue;
            OptValueChain = [OptimalValue_preStep OptimalValue];
%             if (abs(OptimalValue_preStep)<10)
%                 OptValueChain = [OptimalValue_preStep OptimalValue];
%                 OptValueChain_rate.MinUp = [min(sum(UplinkRate_PerGroupPerUser,2))];
%                 OptValueChain_rate.MaxUp = [max(sum(UplinkRate_PerGroupPerUser,2))];
%                 OptValueChain_rate.MinDown = [min(sum(DownlinkRate_PerGroupPerUser,2))];
%                 OptValueChain_rate.MaxDown = [max(sum(DownlinkRate_PerGroupPerUser,2))];
% %                 [sum(sum(UplinkRate_PerGroupPerUser)) + sum(sum(DownlinkRate_PerGroupPerUser))];
%             else
%                 OptValueChain = [OptimalValue];
%                 OptValueChain_rate.MinUp = [min(sum(UplinkRate_PerGroupPerUser,2))];
%                 OptValueChain_rate.MaxUp = [max(sum(UplinkRate_PerGroupPerUser,2))];
%                 OptValueChain_rate.MinDown = [min(sum(DownlinkRate_PerGroupPerUser,2))];
%                 OptValueChain_rate.MaxDown = [max(sum(DownlinkRate_PerGroupPerUser,2))];
% %                 OptValueChain_sum = [sum(sum(UplinkRate_PerGroupPerUser)) + sum(sum(DownlinkRate_PerGroupPerUser))];
%             end

            while (n<MaxIteration)

                disp(['******************* ' strdisp ' --- Iteration: ' num2str(n+1) ' *********************']);
                % Solve (43)

                [ OptimalValue, UplinkRate_PerUser_next, DownlinkRate_PerUser_next, W_next, p_next, addVar_next, Status] = ... 
                Get_optSolutionPerIteration2(Pbs, P, Hu, Hd, G_SI, G_lk, rho, Rate_Threshold, W_current, p_current, addVar_current, {0, 'DL/UL'}, max_min, eta);
                                                
                % Update
                
                
                if (isFeasible(Status, addVar_next, OptimalValue_current, OptimalValue, sum(UplinkRate_PerUser_next/log(2),2), sum(DownlinkRate_PerUser_next/log(2),2), Rate_Threshold, n))


                    W_current =	W_next;

                    p_current = p_next;
                    
                    addVar_current = addVar_next;
                    
                    
                    
                    UplinkRate_PerUser = UplinkRate_PerUser_next;
                    
                    DownlinkRate_PerUser = DownlinkRate_PerUser_next;


                    % Check convergence

        %             if (abs(OptimalValue-OptimalValue_current)<10^-2)
                    if (checkConvergence(OptValueChain, OptimalValue, Converg))
                        OptValueChain = [OptValueChain repmat(OptimalValue_current, 1, MaxIteration-length(OptValueChain) )];
%                         OptValueChain_rate.MinUp = [OptValueChain_rate.MinUp repmat(OptValueChain_rate.MinUp(length(OptValueChain_rate.MinUp)), 1, MaxIteration-length(OptValueChain_rate.MinUp) )];
%                         OptValueChain_rate.MaxUp = [OptValueChain_rate.MaxUp repmat(OptValueChain_rate.MaxUp(length(OptValueChain_rate.MaxUp)), 1, MaxIteration-length(OptValueChain_rate.MaxUp) )];
%                         OptValueChain_rate.MinDown = [OptValueChain_rate.MinDown repmat(OptValueChain_rate.MinDown(length(OptValueChain_rate.MinDown)), 1, MaxIteration-length(OptValueChain_rate.MinDown) )];
%                         OptValueChain_rate.MaxDown = [OptValueChain_rate.MaxDown repmat(OptValueChain_rate.MaxDown(length(OptValueChain_rate.MaxDown)), 1, MaxIteration-length(OptValueChain_rate.MaxDown) )];
%                         OptValueChain_sum = [OptValueChain_sum repmat(OptimalValue_current, 1, MaxIteration-length(OptValueChain_sum) )];
                        break;
                    else
                        OptimalValue_current = OptimalValue
                        OptValueChain = [OptValueChain OptimalValue_current];
%                         OptValueChain_rate.MinUp = [OptValueChain_rate.MinUp min(sum(UplinkRate_PerGroupPerUser,2))];
%                         OptValueChain_rate.MaxUp = [OptValueChain_rate.MaxUp max(sum(UplinkRate_PerGroupPerUser,2))];
%                         OptValueChain_rate.MinDown = [OptValueChain_rate.MinDown min(sum(DownlinkRate_PerGroupPerUser,2))];
%                         OptValueChain_rate.MaxDown = [OptValueChain_rate.MaxDown max(sum(DownlinkRate_PerGroupPerUser,2))];
%                         OptValueChain_sum = [OptValueChain_sum (sum(sum(UplinkRate_PerGroupPerUser)) + sum(sum(DownlinkRate_PerGroupPerUser)))];
                    end

                else

                    disp('Infeasible point --> keeping the latest feasible point');
%                     OptValueChain = [OptValueChain OptimalValue_current];
%                     OptValueChain_sum = [OptValueChain_sum (sum(sum(UplinkRate_PerGroupPerUser)) + sum(sum(DownlinkRate_PerGroupPerUser)))];
                    OptValueChain = [OptValueChain repmat(OptimalValue_current, 1, MaxIteration-length(OptValueChain) )];
%                     OptValueChain_rate.MinUp = [OptValueChain_rate.MinUp repmat(OptValueChain_rate.MinUp(length(OptValueChain_rate.MinUp)), 1, MaxIteration-length(OptValueChain_rate.MinUp) )];
%                     OptValueChain_rate.MaxUp = [OptValueChain_rate.MaxUp repmat(OptValueChain_rate.MaxUp(length(OptValueChain_rate.MaxUp)), 1, MaxIteration-length(OptValueChain_rate.MaxUp) )];
%                     OptValueChain_rate.MinDown = [OptValueChain_rate.MinDown repmat(OptValueChain_rate.MinDown(length(OptValueChain_rate.MinDown)), 1, MaxIteration-length(OptValueChain_rate.MinDown) )];
%                     OptValueChain_rate.MaxDown = [OptValueChain_rate.MaxDown repmat(OptValueChain_rate.MaxDown(length(OptValueChain_rate.MaxDown)), 1, MaxIteration-length(OptValueChain_rate.MaxDown) )];
%                     OptValueChain_sum = [OptValueChain_sum repmat(OptimalValue_current, 1, MaxIteration-length(OptValueChain_sum) )];
                    break;

                end

                n = n + 1;

            end
            
            disp('-------------------- Converging point ----------------------------');
            
            UplinkRate_PerUser = UplinkRate_PerUser/log(2)
            DownlinkRate_PerUser = DownlinkRate_PerUser/log(2)
            
            
            
            
            UplinkPower_PerUser = sum(p_current.^2,2);
            Sum_Downlink_power = real(trace(W_current*W_current'));
                                
            OptValueChain = OptValueChain/log(2)
            
%             OptValueChain_rate.MinUp = OptValueChain_rate.MinUp/log(2);
%             OptValueChain_rate.MaxUp = OptValueChain_rate.MaxUp/log(2);
%             OptValueChain_rate.MinDown = OptValueChain_rate.MinDown/log(2);
%             OptValueChain_rate.MaxDown = OptValueChain_rate.MaxDown/log(2);
            
%             OptValueChain_sum = OptValueChain_sum/log(2);
            
            
        end

        disp('------------------------------------------------------------------');
            
        if (~isBreak && global_OptValue<OptimalValue_current)

            
            global_OptValue = OptimalValue_current;
            
            global_Solution = {W_current, p_current, addVar_current};

            global_DLRate_PerUser = DownlinkRate_PerUser;

%             global_Sum_Downlink_power = real(trace(W_current*diag(time_current)*W_current'));
%             global_Sum_Downlink_power = real(trace(W_current*W_current'));

            
            global_ULRate_PerUser = UplinkRate_PerUser;

%             global_UplinkPower_PerGroupPerUser = p_current.^2;
%             global_UplinkPower_PerUser = global_UplinkPower_PerGroupPerUser*time_current';
%             global_UplinkPower_PerUser = sum(global_UplinkPower_PerGroupPerUser,2);

            global_OptValueChain = OptValueChain;
%             if (abs(OptValueChain_rate.MinUp(end)-OptValueChain_rate.MaxUp(end))<10-1 && abs(OptValueChain_rate.MinDown(end)-OptValueChain_rate.MaxDown(end))<10^-1)
%             global_OptValueChain_rate = OptValueChain_rate;
%             end
%             global_OptValueChain_sum = OptValueChain_sum;

        end


        if (global_OptValue>0)

            disp([' ################### FINAL RESULT -- ' strdisp ' #############################']);

            global_OptValue = global_OptValue/log(2)
            
            global_ULRate_PerUser
            global_DLRate_PerUser           
            
            global_OptValueChain
            disp(' ##############################################################');
        end

end

function [check] = isFeasible(Status, addVar_next, OptimalValue_current, OptimalValue, UplinkRate, DownlinkRate, Threshold, Step)

% check if the current point is feasible and the next point is infeasible

checkStatus = 1;
if (~isempty(findstr('Infeasible', Status)))
    checkStatus = 0;
    disp('Out of feasible region --> Not pass range check');
end

checkAssignVar = 1;
if (sum(sum(addVar_next{1}<0))>0 && sum(sum(addVar_next{2}<0))>0)
    checkAssignVar = 0;
    disp('Assignment Variables --> Not pass range check');
end

checkObjVal = 1;

if ((OptimalValue_current>0 && OptimalValue<OptimalValue_current && Step>0) || (round(OptimalValue_current,20)==0) || (OptimalValue_current>1000))
    checkObjVal = 0;
    disp('OptimalValue --> Not pass range check');
    OptimalValue_current
    OptimalValue
end

% check = ((Step<20)||(checkStatus && checkAssignVar && checkObjVal));
check = (checkStatus && checkAssignVar && checkObjVal);


end


function [check] = checkConvergence(OptValueChain, OptimalValue, Converg)

check = 0;
% if (OptimalValue>80)
%     check=1;
% end
if (length(OptValueChain)>10)
    if (abs(OptimalValue-OptValueChain(length(OptValueChain))) < 10^-3)
        check = 1;
    end
    if (abs(OptimalValue-OptValueChain(length(OptValueChain)-5)) < 0.01)
        check = 1;
    end
    if (~Converg)
    if (abs(OptimalValue-OptValueChain(length(OptValueChain)))/abs(OptimalValue) < 0.01)
        check = 1;
    end
    end
end

end

