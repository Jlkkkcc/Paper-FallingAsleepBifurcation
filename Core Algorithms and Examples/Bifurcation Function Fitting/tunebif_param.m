function [params_tuned,rsq_init,rsq_final,dd,x_ini,iffail]= tunebif_param(params_ref,xx_smooth,tvec)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: JL

% Fine-tune initial parameteres of the harvesting model for falling asleep
% in 30-minutes dynamics

% The algorithm aligns maximum dip to t=0;

% The algorithm also fine tunes a bit the K values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(params_ref)
    params_ref = [6,6.482,0.8,0.34];
end

% params_ref = [6,6.482,0.8,0.34];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Loops for fine-tune K
max_loop_K = 10;
nlp_K = 0;
iffail = 0;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Loops for fine-tune speed
tune_step_r = 0.02;
tune_step_m = 0.005;
max_loop = 200;

params_adjust = params_ref;
params_prev = params_adjust;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start loop

while nlp_K < max_loop_K

    K_actual = mean(xx_smooth(1:100))+0.25*nlp_K;
    K_ref = params_adjust(2);

    % Initial adjust
    r_new = params_ref(1)*K_ref/K_actual;

    if r_new >20    % Add upper constraint
        r_new = 20;
    end
    if r_new <1     % Add upper constraint
        r_new = 1;
    end

    params_adjust = [r_new,K_actual,params_adjust(3),params_adjust(4)];

    x_ini = [mean(xx_smooth(1:20));1];
    [~,dd] = ode45(@(t,x) harvest(t,x,params_adjust),tvec,x_ini);

    % figure
    % % plot(t,dd(:,1),'-o',t,dd(:,2),'-.')
    % plot(t,dd(:,1),'-o')
    % hold on
    % plot(t,xx_smooth)
    if size(dd,1)<length(xx_smooth)

        if nlp_K == 0
            iffail = 1;
            rsq_final = NaN;
            rsq_init = NaN;
            break
        else
            params_adjust = params_prev;
            break
        end
    end

    rsq_init = rsquare(xx_smooth,dd(:,1));

    % rsq_init

    n_loop = 0;

    %%%%%%%%%%%%%%%%
    % Loop tune speed
    tzero = find(tvec==0);
    dd_diff = diff(dd(:,1));
    min_diff_idx = find(dd_diff == min(dd_diff))+1;
    % Avoid an initial dip due to large K
    trimcut = 0;
    if min_diff_idx<20
        dd_diff = dd_diff(21:end);
        trimcut = 20;
        min_diff_idx = find(dd_diff == min(dd_diff))+1+trimcut;
    end
    if min(dd_diff)>-0.003
        min_diff_idx = length(dd_diff)+1+trimcut;
    end

    while abs(min_diff_idx-tzero)>10 %Check difference in time gaps
        params_ini_Beforeloop = params_adjust;
        if min_diff_idx>tzero % Positive difference

            r_new = params_adjust(1)-tune_step_r;
            params_adjust(1) = r_new;

        elseif min_diff_idx<tzero

            r_new = params_adjust(1)+tune_step_r;
            params_adjust(1) = r_new;

        end
        if r_new<=0    % Ensure that r_new is larger than 0 and tune m
            r_new = 0.5;
            params_adjust(1) = r_new;

            if min_diff_idx>tzero % Positive difference

                m_new = params_adjust(4)+tune_step_m;
                params_adjust(4) = m_new;

            elseif min_diff_idx<tzero

                m_new = params_adjust(4)+tune_step_m;
                params_adjust(4) = m_new;

            end

        end
        [~,dd] = ode45(@(t,x) harvest(t,x,params_adjust),tvec,x_ini);
        if size(dd,1)<length(xx_smooth)
            params_adjust = params_ini_Beforeloop;
            break
        end
%         rsq_new = rsquare(xx_smooth,dd(:,1));
        dd_diff_new = diff(dd(:,1));
        min_diff_idx = find(dd_diff_new == min(dd_diff_new))+1;
        % Avoid an initial dip due to large K
        trimcut = 0;
        if min_diff_idx<20
            dd_diff_new = dd_diff_new(21:end);
            trimcut = 20;
            min_diff_idx = find(dd_diff_new == min(dd_diff_new))+1+trimcut;
        end
        if min(dd_diff_new)>-0.003
            min_diff_idx = length(dd_diff_new)+1+trimcut;
        end

        % rsq_diff = rsq_new-rsq_prev;
        % rsq_prev = rsq_new;
        %
        % if rsq_diff<rsq_change  %No improvement
        %
        %     rsq_notimprove = rsq_notimprove+1;
        % end
        n_loop = n_loop+1;
        % params_prev = params_adjust;
        if n_loop>max_loop % Maximum loop
            break
        end
    end

    [~,dd] = ode45(@(t,x) harvest(t,x,params_adjust),tvec,x_ini);
    % figure
    % plot(t,dd(:,1),'-o')
    % hold on
    % plot(t,xx_smooth)
    rsq_final = rsquare(xx_smooth,dd(:,1));

    % See if K improves things
    if nlp_K == 0
        rsq_K_prev = rsq_final;
        params_prev = params_adjust;
    else
        rsq_Kdiff = rsq_final - rsq_K_prev;
        if rsq_Kdiff<0   % Worser r2: change back
            params_adjust = params_prev;
        else
            rsq_K_prev = rsq_final;
            params_prev = params_adjust;
        end
    end
    nlp_K = nlp_K+1;
    if nlp_K>max_loop_K
        break
    end

end

if iffail
    params_tuned = params_ref;
else
    params_tuned = params_adjust;
end
[~,dd] = ode45(@(t,x) harvest(t,x,params_tuned),tvec,x_ini);

rsq_final = rsquare(xx_smooth,dd(:,1));





























