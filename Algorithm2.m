% % % The MATLAB CODE is used for the following paper: Hieu V. Nguyen, Van-Dinh Nguyen, Octavia A. Dobre, Diep N. Nguyen, Eryk Dutkiewicz, and Oh-Soon Shin, 
% % % "Joint Power Control and User Association for NOMA-Based Full-Duplex Systems,"
% % % IEEE Transactions on Communications, vol. 67, no. 11, pp. 8037-8055, Nov. 2019.


function [ GetBreak, OptimalValue_final, OptValueChain, DownlinkRate_PerGroupPerUser, UplinkRate_PerGroupPerUser, time_current ] = Algorithm2( Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, rho, Fixed_timegroup_assignment, Rate_Threshold )
%ALGORITHM2 Summary of this function goes here
%   Detailed explanation goes here

% global MaxIteration
MaxIteration = 50;


strdisp = ['For SI = ' num2str(10*log10(rho)) ' dB'];

if (Fixed_timegroup_assignment)
    strdisp = ['Fixed time -- ' strdisp];
else
    strdisp = ['Opt. time -- ' strdisp];
end

disp(['***** Getting a feasible point ..... ' strdisp]);

[ GetBreak, OptimalValue, OptimalValue_preStep, DownlinkRate_PerGroupPerUser, UplinkRate_PerGroupPerUser, RDown_current, RThDown_current, RUp_current, RTh_current, W_current, p_current, phi_current, time_current, w_tilde_current, p_bar_current, alpha_bar_current, alpha, beta_bar_current, beta] = GetInitialization(Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, rho, Fixed_timegroup_assignment, Rate_Threshold);




%% Loop for optimization

if (GetBreak)
    
    GetBreak = 1;
    OptimalValue_final = 0;
    OptValueChain = [];
    
else
    
    GetBreak = 0;

    disp(['---------------------------------- ' strdisp ' - Run iterative algorithm --------------------------------']);

    n = 0;
    OptimalValue_current = OptimalValue;
    if (abs(OptimalValue_preStep)<10)
        OptValueChain = [OptimalValue_preStep OptimalValue];
    else
        OptValueChain = [OptimalValue];
    end

    %     RDown_current = 0;
    %     RThDown_current = zeros(K,G);
    %     alpha_bar_current = zeros(K,G);
    %     mu_current = zeros(K,G);
    % 
    %     for g = 1:1:G
    %         for k = 1:1:K
    %             PerGroupPerUserRate = log2(1 + real(H(:,k)'*W_current((k-1)*N_tx+1:k*N_tx,g))^2/phi_current(k,g)^2);
    %             alpha_bar_current(k,g) = sqrt(alpha_current(k,g)*PerGroupPerUserRate);
    %     %         mu_current(k,g) = norm(W_current((k-1)*N_tx+1:k*N_tx,g))^2/alpha_current(k,g);
    %             mu_current(k,g) = Pbs/K;
    %             RDown_current = RDown_current + time_current(g)*alpha_bar_current(k,g)^2;
    %     %         RDown_current = RDown_current + time_current(g)*PerGroupPerUserRate;
    %             RThDown_current(k,g) = sqrt(time_current(g)*alpha_bar_current(k,g)^2);
    %         end
    %     end
    % 
    %     RUp_current = 0;
    %     RTh_current = zeros(L,G);
    %     beta_bar_current = zeros(L,G);
    % 
    %     for g = 1:1:G
    %         for l = 1:1:L
    %             PerGroupPerUserRate = log2(1 + GetSINR(G_channel(:,l:L), G_SI, W_current(:,g), p_current(l:L,g)));
    %             beta_bar_current(l,g) = sqrt(beta_current(l,g)*PerGroupPerUserRate);
    %             RUp_current = RUp_current + time_current(g)*beta_bar_current(l,g)^2;
    %             RTh_current(l,g) = sqrt(time_current(g)*beta_bar_current(l,g)^2);
    %         end
    %     end
    % 
    %     w_tilde_current  = sum(W_current.^2);
    %     p_bar_current = repmat(P, 1, G)/G*inv(diag(time_current));

    while (n<MaxIteration)

        disp(['******************* ' strdisp ' --- Iteration: ' num2str(n+1) ' *********************']);
        % Solve (45)

    % 	[W_next, p_next, phi_next, eta_next, z_next, z_bar_next, kapa_next, upsilon_next, upsilon_bar_next, alpha_next, beta_next, time_next, tau_next, OptimalValue] = Get_optSolutionPerIteration(Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, eta_current, z_current, z_bar_current, kapa_current, upsilon_current, upsilon_bar_current, alpha_current, beta_current, time_current, tau_current, mu_current, theta_current);
    % 	[W_next, p_next, phi_next, z_next, z_bar_next, upsilon_next, upsilon_bar_next, time_next, x_next, y_next, OptimalValue] = Get_optSolutionPerIteration2(Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, z_current, z_bar_current, upsilon_current, upsilon_bar_current, alpha_current, beta_current, time_current, x_current, y_current, mu_current, theta_current);
        [ OptimalValue, DownlinkRate_PerGroupPerUser_next, UplinkRate_PerGroupPerUser_next, RDown_next, RThDown_next, RUp_next, RTh_next, W_next, p_next, phi_next, time_next, w_tilde_next, p_bar_next, alpha_bar_next, alpha_next, beta_bar_next, beta_next] = Get_optSolutionPerIteration4(Pbs, P, H, G_channel, G_SI, G_hat, W_current, p_current, phi_current, time_current, w_tilde_current, p_bar_current, RDown_current, RThDown_current, RUp_current, RTh_current, alpha_bar_current, beta_bar_current, rho, Fixed_timegroup_assignment, Rate_Threshold);

        % Update
        
        if (isFeasible(alpha, beta, time_current, alpha_next, beta_next, time_next, OptimalValue_current, OptimalValue))

            RDown_current = RDown_next;

            RThDown_current = RThDown_next;

            alpha_bar_current = alpha_bar_next;

    %         alpha_current = alpha_next;

        %         mu_current = mu_next;

            alpha = alpha_next;

            RUp_current = RUp_next;

            RTh_current = RTh_next;

            beta_bar_current = beta_bar_next;

            beta = beta_next;

            W_current =	W_next;

            Sum_Downlink_power = real(trace(W_current*diag(time_current)*W_current'));

            DownlinkRate_PerGroupPerUser = DownlinkRate_PerGroupPerUser_next*diag(time_next);

            DownlinkRate_PerUser = sum(DownlinkRate_PerGroupPerUser,2);

            DownlinkPower_PerUserPerGroup = DownlinkPower_Analysis(W_current, time_current);

            UplinkRate_PerGroupPerUser = UplinkRate_PerGroupPerUser_next*diag(time_next);

            UplinkRate_PerUser = sum(UplinkRate_PerGroupPerUser,2);


            p_current = p_next;

            phi_current = phi_next;

            time_current = time_next;

            w_tilde_current = w_tilde_next;

            p_bar_current = p_bar_next;


            % Check convergence

%             if (abs(OptimalValue-OptimalValue_current)<10^-2)
            if (checkConvergence(OptValueChain, OptimalValue))
                break;
            else
                OptimalValue_current = OptimalValue
                OptValueChain = [OptValueChain OptimalValue_current];
            end
            
        else
            
            disp('Infeasible point --> keeping the latest feasible point');
            OptValueChain = [OptValueChain OptimalValue_current];
            break;
        
        end
        
        n = n + 1;

    end

    disp([' ################### FINAL RESULT -- ' strdisp ' #############################']);
    alpha = round(alpha)
    beta = round(beta)
    time_current
    OptimalValue_final = 1/log(2)*OptimalValue_current
    OptValueChain = OptValueChain/log(2);
    DownlinkRate_PerGroupPerUser = DownlinkRate_PerGroupPerUser/log(2)
    DownlinkRate_PerUser = sum(DownlinkRate_PerGroupPerUser,2)
    Sum_Downlink_power = real(trace(W_current*diag(time_current)*W_current'))
    UplinkRate_PerGroupPerUser = UplinkRate_PerGroupPerUser/log(2)
    UplinkRate_PerUser = sum(UplinkRate_PerGroupPerUser,2)
    UplinkPower_PerGroupPerUser = p_current.^2
    UplinkPower_PerUser = UplinkPower_PerGroupPerUser*time_current'
    OptValueChain
    disp(' ##############################################################');

end

end

function [check] = isFeasible(alpha, beta, time, alpha_next, beta_next, time_next, OptimalValue_current, OptimalValue)

% check if the current point is feasible and the next point is infeasible

check_alpha = [round(alpha)>1] + [round(alpha)<0];
check_beta = [round(beta)>1] + [round(beta)<0];

check_time = ~(round(sum(time),3)==1.000);

check1 = [(sum(sum(check_alpha)) + sum(sum(check_beta)) + check_time) == 0]; % check1 = 1 means that current point is feasible

check_alpha_next = [round(alpha_next)>1] + [round(alpha_next)<0];
check_beta_next = [round(beta_next)>1] + [round(beta_next)<0];

check_time_next = ~(round(sum(time_next),3)==1.000);

check2 = [(sum(sum(check_alpha_next)) + sum(sum(check_beta_next)) + check_time_next) == 0]; % check2 = 1 means that next point is feasible

check3 = 1;

if ((OptimalValue_current>0 && OptimalValue<0) || (round(OptimalValue_current,5)==0) || (OptimalValue_current>1000))
    check3 = 0;
end

check = (~(check1 && ~check2) && check3);

end

function [check] = checkConvergence(OptValueChain, OptimalValue)

check = 0;

if (length(OptValueChain)>10)
    if (abs(OptimalValue-OptValueChain(length(OptValueChain))) < 10^-5)
        check = 1;
    end
    if (abs(OptimalValue-OptValueChain(length(OptValueChain)-5)) < 0.1)
        check = 1;
    end
end

end
