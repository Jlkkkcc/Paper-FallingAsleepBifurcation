%% Individual fitting modifications Mar2024

% Revised without optimisation

% Author: JL

%% Load individual traces

clear all
clc
load TS_SleepDistIndividuals_indivft_ftadd.mat

%% Pre-process  / Revised pipeline for fitting

max_timeOnset = 30*60;     % Maximum time to sleep onset in seconds
max_epc_beforeOnset = floor(max_timeOnset/jmp_len)-1;

xoverFcn = @(T, Y) MyEventFunction(T, Y);
options = odeset('RelTol',1e-15,'AbsTol',1e-15,'Events',xoverFcn);

ft_touse = [15,46,47,53,2,9,33,41,54];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Storage of results
len_allsbjs_perft = [];
rsq_all_perft = [];
minparam_all_perft = [];
rmse_all_perft = [];
fitsuccess_perft = [];

mdl_ts_perft = [];
c_ts_perft = [];
xini_opt_all_perft = [];
sdynam_ts_perft = [];
tvec_all_perft = [];

rst_preopt_allsbj_perft = [];
rst_afteropt_allsbj_perft = [];

iftrim_perft = [];
Trimmed_epcs_perft = [];

for nft = 1:length(ft_touse)

    ftid = ft_touse(nft);
    clear ts_sbj_nomed
    ts_sbj_nomed = ts_perft_allsbj{nft};

    % Control parameters
    ts_touse = ts_sbj_nomed;
    ifsmooth = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Storage of results
    len_allsbjs = NaN(num_sbjs,1);
    rsq_all = NaN(num_sbjs,1);

    minparam_all = [];
    rmse_all = NaN(num_sbjs,1);

    fitsuccess = NaN(num_sbjs,1);

    mdl_ts = [];
    c_ts = [];
    xini_opt_all = [];
    sdynam_ts = [];
    tvec_all = [];

    rst_preopt_allsbj = NaN(num_sbjs,3);
    rst_afteropt_allsbj = NaN(num_sbjs,3);

    iftrim = NaN(num_sbjs,1);
    Trimmed_epcs = NaN(num_sbjs,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If failed optimisation setting
    rsq_th = 0.9;
    ifskip = NaN(num_sbjs,1);
    iffail_opt = NaN(num_sbjs,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fittings
    for nsbj = 1:num_sbjs

        % % If already fitted
        % if ismember(nsbj,sbjs_fitted)
        %     continue
        % end

        clear ts_now xx_smooth tvec_now dd dd_2 rsq_init rsq_final mse_preopt mse_final rsq_final2 rsq_alltrials params_tuned iffail

        if ifelig(nsbj)        % Eligible for fittings

            ts_now = ts_touse{nsbj}.ftdist_noart;
            len_thissbj = length(ts_now);
            len_allsbjs(nsbj) = len_thissbj;

            % Whether length of onset larger than maximum required
            if (len_thissbj-epc_postasleep) > max_epc_beforeOnset
                gap_this = len_thissbj-epc_postasleep - max_epc_beforeOnset;
                len_now1 = max_epc_beforeOnset + epc_postasleep;
                ts_now = ts_now(gap_this+1:end);
                tvec_now = -(len_now1-201)/20:0.05:10;
            else
                len_now1 = len_thissbj;
                tvec_now = -(len_now1-201)/20:0.05:10;
            end

            % Smooth data
            if ifsmooth
                xx_smooth = smoothdata(ts_now,2,'movmedian',[20,0]);
            else
                xx_smooth = ts_now;
            end
            xmin = min(xx_smooth);
            xx_smooth = xx_smooth - xmin;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Try Fitting

            %%%%%%%%%%%%%%%%%%
            % Cut out tails % Trim uniformly 20 epochs if so
            threshold = mean(xx_smooth(20:60)); % Define the artifact threshold
            artifact_tail = (sum(xx_smooth(1:20)>(threshold*1.5)));
            if artifact_tail
                iftrim(nsbj) = 1;
                Trimmed_epcs(nsbj) = 20;
                xx_smooth = xx_smooth(21:end);
                tvec_now = tvec_now(21:end);

            end

            %%%%%%%%%%%%%%%%%%
            % Step 1: Fine-tune the initial params
            [params_tuned,rsq_init,rsq_final,dd,x_ini_tuned,iffail] = tunebif_param([],xx_smooth,tvec_now);
            if iffail             % Complete failure this subject
                fitsuccess(nsbj) = 0;
                disp(['Subject ' num2str(nsbj) ' completely failed'])
                continue
            end
            if rsq_final > rsq_th        % if already good enough skip optimisation
                fitsuccess(nsbj) = 1;
                ifskip(nsbj) = 1;

                sdynam_ts{nsbj} = xx_smooth;

                rsq_all(nsbj) = rsq_final;
                minparam_all{nsbj} = params_tuned;
                rmse_all(nsbj) = mean((dd(:,1) - xx_smooth').^2);
                rst_preopt_allsbj(nsbj,:) = [rsq_final,rmse_all(nsbj),rsq_final-rmse_all(nsbj)];
                rst_afteropt_allsbj(nsbj,:) = [rsq_final,rmse_all(nsbj),rsq_final-rmse_all(nsbj)];

                xini_opt_all{nsbj} = x_ini_tuned;
                mdl_ts{nsbj} = dd(:,1);
                c_ts{nsbj} = dd(:,2);
                tvec_all{nsbj} = tvec_now;
                disp(['Subject ' num2str(nsbj) ' initial tune best'])
                continue
            end

            %%%%%%%%%%%%%%%%%%
            % Step 2: Further fine tuning without optimisation
            [params_optim,rsq_preopt,mse_preopt,rsq_final2,mse_final,dd_2,iffail1] = finetune_bif_rev(params_tuned,xx_smooth,tvec_now,x_ini_tuned,[]);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % If optimisation failed
            iffail_opt(nsbj) = iffail1;
            if iffail1
                disp(['Subject ' num2str(nsbj) ' optimisation failed, taking initial'])
            end

            rst_preopt_allsbj(nsbj,:) = [rsq_preopt,mse_preopt,rsq_preopt-mse_preopt];
            rst_afteropt_allsbj(nsbj,:) = [rsq_final2,mse_final,rsq_preopt-mse_final];

            if rst_afteropt_allsbj(nsbj,3)<rst_preopt_allsbj(nsbj,3)
                error('Check function')
            end

            xini_opt_all{nsbj} = x_ini_tuned;
            rsq_all(nsbj) = rsq_final2;
            minparam_all{nsbj} = params_optim;
            rmse_all(nsbj) = mse_final;

            mdl_ts{nsbj} = dd_2(:,1);
            c_ts{nsbj} = dd_2(:,2);

            sdynam_ts{nsbj} = xx_smooth;
            tvec_all{nsbj} = tvec_now;

            fitsuccess(nsbj) = 1;
            disp(['Subject ' num2str(nsbj) ' successfully fitted with optimisation'])


        end

    end

    % Save per feature
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Storage of results
    len_allsbjs_perft{nft} = len_allsbjs;
    rsq_all_perft{nft} = rsq_all;
    minparam_all_perft{nft} = minparam_all;
    rmse_all_perft{nft} = rmse_all;
    fitsuccess_perft{nft} = fitsuccess;

    mdl_ts_perft{nft} = mdl_ts;
    c_ts_perft{nft} = c_ts;
    xini_opt_all_perft{nft} = xini_opt_all;
    sdynam_ts_perft{nft} = sdynam_ts;
    tvec_all_perft{nft} = tvec_all;

    rst_preopt_allsbj_perft{nft} = rst_preopt_allsbj;
    rst_afteropt_allsbj_perft{nft} = rst_afteropt_allsbj;

    iftrim_perft{nft} = iftrim;
    Trimmed_epcs_perft{nft} = Trimmed_epcs;

end

save('MdlFit_Sdynamic_allIndivs_perft.mat')

%% Computation of tipping points and reshape into tables

c_critc = NaN(num_sbjs,1);
t_critc = NaN(num_sbjs,1);
c3sol_all = [];
xsol_allt_allsbj = [];
ifproblem = NaN(num_sbjs,1);    % Error in C, meaning the fitting was not well constraint
ifNoC = NaN(num_sbjs,1);
ifTransC = NaN(num_sbjs,1);

epcs_score_pre_max = floor(max_timeOnset/score_size);
epcs_post_score_now = 20;
max_epcs_slp = epcs_score_pre_max+epcs_post_score_now;

savedir_pre = '/mnt/Storage_JL/JL Analysis/UK Dementia Research Institute Dropbox/Junheng Li/Data proc and analysis/Analysis MESA/Dynamical and state analysis/Dynamical systems and Early warning signals/Revision Paper Dec24/Indiv Bif fitting Ftadd Jan25/Indiv Bifdiagrams';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Only does it when there is a valid R2

for nsbj = 1:num_sbjs
    
    if rsq_all(nsbj) > 0   % Valid R2
        
        ifNoC(nsbj) = 0;
        ifproblem(nsbj) = 0;
        ifTransC(nsbj) = 0;

        clear c param_min xsol_allt tvec_this
        c = c_ts{nsbj};
        param_min = minparam_all{nsbj};
        tvec_this = tvec_all{nsbj};

        r = param_min(1);
        K = param_min(2);
        h= param_min(3);
        m = param_min(4);

        % x_theo = zeros(length(c));
        c_3sol = [];

        f=figure('Visible','off');
        hold on
        xsol_allt = [];

        syms x

        for idxc = 1:length(c)

            % Roots of the theoretical equation when dx/dt = 0;
            %     x = roots([1/10, -1, (c(idxc)+1/10), -1,0]);
            eqn = r*x*(1-x/K) - c(idxc)*x^2 / (x^2+h^2);
            xsol_all = solve(eqn==0);
            xsol_all = double(xsol_all);

            xsol = [];
            for jjj = 1:length(xsol_all)
                if isreal(xsol_all(jjj))
                    if abs(xsol_all(jjj)) ~= 0
                        xsol = [xsol,xsol_all(jjj)];
                    end
                end
            end
            % Find 3 solution points
            if length(xsol) == 3
                c_3sol = [c_3sol,c(idxc)];
            end

            scatter(c(idxc)*ones(1,length(xsol)),xsol,2,'k')

            xsol_allt{idxc} = xsol;

        end
        box off
        xlabel('Control variable')
        ylabel('System state')
        set(gca,'FontSize', 12)
        set(gca,'TickDir','out')
        set(gca,'ticklength',3*get(gca,'ticklength'))
        set(gca,'lineWidth',2)

        xsol_allt_allsbj{nsbj} = xsol_allt;
        c3sol_all{nsbj} = c_3sol;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check critical values
        if ~isempty(c_3sol) && (length(c_3sol)<length(c) )  % No solution or full solution

            idx_critic = find(c == c_3sol(end));   % Critical point
            t_critc(nsbj) = tvec_this(idx_critic);
            c_critc(nsbj) = c(idx_critic);
            disp(['Subject',num2str(nsbj),' Critical time ', num2str(t_critc(nsbj)),' min'])
        else
            if length(c_3sol) == length(c)      % Failure situation
                ifproblem(nsbj) = 1;
                disp(['Subject',num2str(nsbj),' Bad bif diagram'])
            end
            % Run transcritical detections
            sfitnow = mdl_ts{nsbj};
            sfitnow = sfitnow(21:end);
            diff_sfitnow = diff(sfitnow);
            min_deriv_idx = find(diff_sfitnow == min(diff_sfitnow));
            idx_critic = min_deriv_idx+1+20;

            % No Critical point at all
            if abs(min(diff_sfitnow))<0.005 || (idx_critic >= length(sfitnow))
                ifNoC(nsbj) = 1;
                disp(['Subject',num2str(nsbj),' No Critical point detected at all'])
                t_critc(nsbj) = NaN;
                c_critc(nsbj) = NaN;
            else 
                ifTransC(nsbj) = 1;    % Transcritical case
                t_critc(nsbj) = tvec_this(idx_critic);
                c_critc(nsbj) = c(idx_critic);
                disp(['Subject',num2str(nsbj),' Transcritical, taking minimum time ', num2str(t_critc(nsbj)),'min'])
            end
        end
        % Construct the title with line breaks
        graphTitle = ['Time Critical: ',num2str(t_critc(nsbj)),' min'];
        % Set the title
        title(graphTitle, 'Interpreter', 'none');
        box off
        savedir = [savedir_pre,'/Sbj',num2str(nsbj)];
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
        figname_svg = fullfile(savedir,['BifDiagram_SbjNo',num2str(nsbj),'.svg']);
        figname1 = fullfile(savedir,['BifDiagram_SbjNo',num2str(nsbj),'.fig']);
        set(f,'CreateFcn','set(gcbo,"visible","on")');
        savefig(f,figname1)
        saveas(f,figname_svg)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Make some further plots
        clear s_raw sfit
        s_raw = sdynam_ts{nsbj};
        sfit = mdl_ts{nsbj};
        f =figure('Visible','off');
        plot(tvec_this,s_raw,'LineWidth',1.5,'Color','k')
        hold on
        plot(tvec_this,sfit,'LineWidth',3,'Color','b')
        set(gca,'FontSize', 12)
        set(gca,'TickDir','out')
        set(gca,'ticklength',2*get(gca,'ticklength'))
        set(gca,'lineWidth',2)
        ylabel('s')
        xlabel('Time (min)')
        line([0,0],ylim,'LineStyle','--','LineWidth',2,'Color','r');
        % Construct the title with line breaks
        graphTitle = ['R2: ',num2str(rsq_all(nsbj))];
        % Set the title
        title(graphTitle, 'Interpreter', 'none');
        set(f,'CreateFcn','set(gcbo,"visible","on")');
        box off
        savedir = [savedir_pre,'/Sbj',num2str(nsbj)];
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
        figname_svg = fullfile(savedir,['Mdlfit_SbjNo',num2str(nsbj),'.svg']);
        figname1 = fullfile(savedir,['Mdlfit_SbjNo',num2str(nsbj),'.fig']);
        savefig(f,figname1)
        saveas(f,figname_svg)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot hypnograms
        scrs_this = oriscore_asleep{nsbj};
        scrs_this = scrs_this(1:end-1);
        if length(scrs_this)>max_epcs_slp
            scrs_this_now = scrs_this(end-epcs_post_score_now-epcs_score_pre_max+1:end);
        else
            scrs_this_now = scrs_this(1:end);
        end
        tvec_hypno = score_size:score_size:length(scrs_this_now)*score_size;
        tvec_hypno = (tvec_hypno-(length(scrs_this_now)-epcs_post_score_now)*score_size)/60;
        f=figure('Visible','off');
        plot(tvec_hypno,scrs_this_now,'LineWidth',2)
        ylim([0,3])
        xlabel('Time (min)')
        ylabel('Stages')
        set(gca,'FontSize', 12)
        set(gca,'TickDir','out')
        set(gca,'ticklength',2*get(gca,'ticklength'))
        set(gca,'lineWidth',2)
        set(f,'CreateFcn','set(gcbo,"visible","on")');
        box off
        savedir = [savedir_pre,'/Sbj',num2str(nsbj)];
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
        figname_svg = fullfile(savedir,['Hypnogram',num2str(nsbj),'.svg']);
        figname1 = fullfile(savedir,['Hypnogram',num2str(nsbj),'.fig']);
        savefig(f,figname1)
        saveas(f,figname_svg)

        clear f
        close all
    end

end

save('MdlFit_Sdynamic_allIndivs_WithTipping.mat')

%% Reshape into tables

clear tbl_results
tbl_results = [];

varNames = ["SbjID","Raw-S","Smoothed-S","R-squared","Optimal-Params","Min-RMSE","Fit-Suceess","Sleep-stages-original","Sleep-stages-trimmed","Expanded_Hyp","Time-to-asleep-min","Fitted-Model-S","Fitted-ControlP","Time-vector","Initial-x","IfTrim","TrimmedEpcs","IfSkip","IfFailOpt","T-Critic","C_Critic","C_3sol","Xsol_all","IfBifFail","IfNoCritical","IfTransCritical"];
epcs_score_pre_max = floor(max_timeOnset/score_size);
epcs_post_score_now = 20;
max_epcs_slp = epcs_score_pre_max+epcs_post_score_now;
post_gap_now = 1;

ftepc_len = 6;
ftepc_jmp = 3;
score_size = 30;

for nsbj = 1:num_sbjs
    
    if ifelig(nsbj)        % Eligible for fittings

        ts_now = ts_touse{nsbj}.ftdist_noart;
        len_thissbj = length(ts_now);

        xx_smooth = sdynam_ts{nsbj};
        tvec_now = tvec_all{nsbj};
        
        % sleep scores
        scrs_this = oriscore_asleep{nsbj};
        scrs_this = scrs_this(1:end-1);
        if length(scrs_this)>max_epcs_slp
            scrs_this_now = scrs_this(end-epcs_post_score_now-epcs_score_pre_max+1:end);
        else
            scrs_this_now = scrs_this(1:end);
        end
        timeonset_this = ((epcs_to_asleep(nsbj)-1)*score_size) /60;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Expanded hynopgram
        hyp_expand = zeros(length(xx_smooth),1);
        trimepcs = Trimmed_epcs(nsbj);
        if isnan(trimepcs)
            trimepcs = 0;
        end
        for ii = 1:length(xx_smooth)
            scr_idx = floor(((ii-1+trimepcs)*ftepc_jmp +ftepc_jmp)/(score_size+0.001))+1;
            if scr_idx>length(scrs_this_now)
                scr_idx = length(scrs_this_now);
            end
            hyp_expand(ii) = scrs_this_now(scr_idx);
        end

        if isempty(tbl_results)
            tbl_results = table(nsbj,{ts_now},{xx_smooth},rsq_all(nsbj),{minparam_all{nsbj}},rmse_all(nsbj),fitsuccess(nsbj),{scrs_this},{scrs_this_now},{hyp_expand},timeonset_this,{mdl_ts{nsbj}},{c_ts{nsbj}},{tvec_now},{xini_opt_all{nsbj}}, ...
                iftrim(nsbj),Trimmed_epcs(nsbj),ifskip(nsbj),iffail_opt(nsbj),t_critc(nsbj),c_critc(nsbj),{c3sol_all{nsbj}},{xsol_allt_allsbj{nsbj}},ifproblem(nsbj),ifNoC(nsbj),ifTransC(nsbj),'VariableNames',varNames);
        else
            tbl_new = table(nsbj,{ts_now},{xx_smooth},rsq_all(nsbj),{minparam_all{nsbj}},rmse_all(nsbj),fitsuccess(nsbj),{scrs_this},{scrs_this_now},{hyp_expand},timeonset_this,{mdl_ts{nsbj}},{c_ts{nsbj}},{tvec_now},{xini_opt_all{nsbj}}, ...
                iftrim(nsbj),Trimmed_epcs(nsbj),ifskip(nsbj),iffail_opt(nsbj),t_critc(nsbj),c_critc(nsbj),{c3sol_all{nsbj}},{xsol_allt_allsbj{nsbj}},ifproblem(nsbj),ifNoC(nsbj),ifTransC(nsbj),'VariableNames',varNames);
            tbl_results = [tbl_results;tbl_new];
        end

    end

end

save('FitResults_Full.mat','tbl_results')

%% Functions ODE


function [VALUE, ISTERMINAL, DIRECTION] = MyEventFunction(T,Y)
%The event function stops when VALUE == 0 and
%ISTERMINAL==1
%a. Define the timeout in seconds
TimeOut = 10;
%
%b. The solver runs until this VALUE is negative (does not change the sign)
    VALUE = toc-TimeOut;
%c. The function should terminate the execution, so
ISTERMINAL = 1;
%d. The direction does not matter
DIRECTION = 0;
end

