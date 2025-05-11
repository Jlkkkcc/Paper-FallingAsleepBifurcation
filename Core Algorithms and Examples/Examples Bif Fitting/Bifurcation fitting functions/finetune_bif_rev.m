function [params_optim,rsq_init,mse_init,rsq_final,mse_final,dd,iffail] = finetune_bif_rev(params_ini,xx,tvec,xini,nloop)

% Performs further tuning of the bifurcation parameters of the harvesting
% model for 30-minute and less sleep-onset S-variable dynamics;

% Set default
rsq_init = NaN;
rsq_final = NaN;
mse_final = NaN;
mse_init = NaN;
params_optim = params_ini;
iffail = 0;
dd = [];

% Initial r-squared
[~,dd] = ode45(@(t,x) harvest(t,x,params_ini),tvec,xini);
rsq_init = rsquare(xx,dd(:,1));
mse_init = mean((dd(:,1) - xx').^2);

cr_index_init = rsq_init-mse_init;    % A novel criterion index; The larger the better
% flag = (cr_index<0);     % Negative situations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimisation
% By default K, h is assumed to be not changing

r_ini = params_ini(1);
m_ini = params_ini(4);
h_ini = params_ini(3);

rtune_step = 0.02;
mtune_step = 0.005;
htune_step = 0.001;

if isempty(nloop)
    nloop = 250;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start from r tuning;
rsq_inc_r = NaN(1,nloop);
mse_inc_r = NaN(1,nloop);
params_now = params_ini;
% Positive direction
for nlp = 1:nloop

    r_now = r_ini+nlp*rtune_step;
    params_now(1) = r_now;
    [~,ddnow] = ode45(@(t,x) harvest(t,x,params_now),tvec,xini);
    if size(ddnow,1)<length(xx)    %Error message
        continue
    end
    rsq_now = rsquare(xx,ddnow(:,1));
    mse_now = mean((ddnow(:,1) - xx').^2);

    rsq_inc_r(nlp) = rsq_now;
    mse_inc_r(nlp) = mse_now;

end
clear cr_index_inc
cr_index_inc = rsq_inc_r-mse_inc_r;

rsq_dec_r = NaN(1,nloop);
mse_dec_r = NaN(1,nloop);
params_now = params_ini;
% Negative direction
for nlp = 1:nloop

    r_now = r_ini-nlp*rtune_step;
    if r_now <=0   % No negative r
        continue
    end
    params_now(1) = r_now;
    [~,ddnow] = ode45(@(t,x) harvest(t,x,params_now),tvec,xini);
    if size(ddnow,1)<length(xx)    %Error message
        continue
    end
    rsq_now = rsquare(xx,ddnow(:,1));
    mse_now = mean((ddnow(:,1) - xx').^2);

    rsq_dec_r(nlp) = rsq_now;
    mse_dec_r(nlp) = mse_now;

end
clear cr_index_dec
cr_index_dec = rsq_dec_r-mse_dec_r;
clear cr_index_r
cr_index_r = [cr_index_inc,cr_index_dec];

% Find best r
clear rbest
% if sum(cr_index_r>0)>0    % Any positive r2
    best_cridx = find(cr_index_r== max(cr_index_r),1);
    max_crnow = max(cr_index_r);
    if max_crnow<cr_index_init    % Worse than initial
        rbest = r_ini;
    else
        if best_cridx<=nloop   % Increasing direction
            rbest = r_ini+best_cridx*rtune_step;
        elseif best_cridx>nloop   % Decreasing direction
            rbest = r_ini-(best_cridx-nloop)*rtune_step;
        end
    end
% end

if isempty(rbest)
    disp('Initial taken, r-optim failed')
    iffail = 1;
    rsq_final = rsq_init;
    mse_final = mse_init;
    return
end
params_optim(1) = rbest;
[~,dd] = ode45(@(t,x) harvest(t,x,params_optim),tvec,xini);
rsq_final = rsquare(xx,dd(:,1));
mse_final = mean((dd(:,1) - xx').^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rsq_final > 0.9  % No need to tune further m
    return
end
cr_index_new = rsq_final - mse_final;    %Update to previous

% Continue to tune m if necessary
rsq_inc_m = NaN(1,nloop);
mse_inc_m = NaN(1,nloop);
params_now = params_optim;
% Positive direction
for nlp = 1:nloop

    m_now = m_ini+nlp*mtune_step;
    params_now(4) = m_now;
    [~,ddnow] = ode45(@(t,x) harvest(t,x,params_now),tvec,xini);
    if size(ddnow,1)<length(xx)    %Error message
        continue
    end
    rsq_now = rsquare(xx,ddnow(:,1));
    mse_now = mean((ddnow(:,1) - xx').^2);

    rsq_inc_m(nlp) = rsq_now;
    mse_inc_m(nlp) = mse_now;

end
clear cr_index_inc
cr_index_inc = rsq_inc_m - mse_inc_m;

rsq_dec_m = NaN(1,nloop);
mse_dec_m = NaN(1,nloop);
params_now = params_optim;
% Negative direction
for nlp = 1:nloop

    m_now = m_ini-nlp*mtune_step;
    if m_now <=0   % No negative r
        continue
    end
    params_now(4) = m_now;
    [~,ddnow] = ode45(@(t,x) harvest(t,x,params_now),tvec,xini);
    if size(ddnow,1)<length(xx)    %Error message
        continue
    end
    rsq_now = rsquare(xx,ddnow(:,1));
    mse_now = mean((ddnow(:,1) - xx').^2);

    rsq_dec_m(nlp) = rsq_now;
    mse_dec_m(nlp) = mse_now;

end
clear cr_index_dec
cr_index_dec = rsq_dec_m - mse_dec_m;
clear cr_index_r
cr_index_r = [cr_index_inc,cr_index_dec];

% Find best m
clear mbest
% if sum(cr_index_r>0)>0    % Any positive r2
    best_cridx = find(cr_index_r== max(cr_index_r),1);
    max_crnow = max(cr_index_r);
    if max_crnow<cr_index_new    % Worse than initial
        mbest = m_ini;
    else
        if best_cridx<=nloop   % Increasing direction
            mbest = m_ini+best_cridx*mtune_step;
        elseif best_cridx>nloop   % Decreasing direction
            mbest = m_ini-(best_cridx-nloop)*mtune_step;
        end
    end
% end

params_optim(4) = mbest;
[~,dd] = ode45(@(t,x) harvest(t,x,params_optim),tvec,xini);
rsq_final = rsquare(xx,dd(:,1));
mse_final = mean((dd(:,1) - xx').^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Further tune h if necessary

if rsq_final > 0.9  % No need to tune further m
    return
end
cr_index_new = rsq_final - mse_final;    %Update to previous

% Continue to tune h if necessary
rsq_inc_h = NaN(1,nloop);
mse_inc_h = NaN(1,nloop);
params_now = params_optim;
% Positive direction
for nlp = 1:nloop

    h_now = h_ini+nlp*htune_step;
    params_now(3) = h_now;
    [~,ddnow] = ode45(@(t,x) harvest(t,x,params_now),tvec,xini);
    if size(ddnow,1)<length(xx)    %Error message
        continue
    end
    rsq_now = rsquare(xx,ddnow(:,1));
    mse_now = mean((ddnow(:,1) - xx').^2);

    rsq_inc_h(nlp) = rsq_now;
    mse_inc_h(nlp) = mse_now;

end
clear cr_index_inc
cr_index_inc = rsq_inc_h - mse_inc_h;

rsq_dec_h = NaN(1,nloop);
mse_dec_h = NaN(1,nloop);
params_now = params_optim;
% Negative direction
for nlp = 1:nloop

    h_now = h_ini-nlp*htune_step;
    if h_now <=0   % No negative h
        continue
    end
    params_now(3) = h_now;
    [~,ddnow] = ode45(@(t,x) harvest(t,x,params_now),tvec,xini);
    if size(ddnow,1)<length(xx)    %Error message
        continue
    end
    rsq_now = rsquare(xx,ddnow(:,1));
    mse_now = mean((ddnow(:,1) - xx').^2);

    rsq_dec_h(nlp) = rsq_now;
    mse_dec_h(nlp) = mse_now;

end
clear cr_index_dec
cr_index_dec = rsq_dec_h - mse_dec_h;
clear cr_index_r
cr_index_r = [cr_index_inc,cr_index_dec];

% Find best m
clear hbest
% if sum(cr_index_r>0)>0    % Any positive r2
    best_cridx = find(cr_index_r== max(cr_index_r),1);
    max_crnow = max(cr_index_r);
    if max_crnow<cr_index_new    % Worse than initial
        hbest = h_ini;
    else
        if best_cridx<=nloop   % Increasing direction
            hbest = h_ini+best_cridx*htune_step;
        elseif best_cridx>nloop   % Decreasing direction
            hbest = h_ini-(best_cridx-nloop)*htune_step;
        end
    end
% end

params_optim(3) = hbest;
[~,dd] = ode45(@(t,x) harvest(t,x,params_optim),tvec,xini);
rsq_final = rsquare(xx,dd(:,1));
mse_final = mean((dd(:,1) - xx').^2);






