%% EWS evaluations for sleep distance dyanmics (Critical slowing down) on each individual

% Author: JL

%% Load data

clear all
clc
load 'TS_SleepDistIndividuals.mat'

%% Group-level calculation of early warning signals (autocorrelation and variance)

ts_sbj = ts_sbj_nomed;

T_pval = table;
T_tau = table;
T_bds = table;

% EWS parameters
ews_win = 100;   % 5-min running window
ews_surr = 1000;
jmp = epc_len*(1-ovlap);       % Time jump (considering overlap)

for sbj = 1:num_sbjs
    if ifelig(sbj)
        
        clear x t t_sleep
        varname_pval = {'Subject No.'};
        varname_tau = {'Subject No.'};
        varname_bds = {'Subject No.'};

        val_pval = [sbj];
        val_tau = [sbj];
        val_bds = [sbj];

        ft_name = 'DistOnset';

        x = ts_sbj{sbj}.ftdist_noart;

        t_sleep = epcs_to_asleep(sbj)*score_size;
        t = ((stp_ftepcs_all_filt{sbj}-1)/Fs) - t_sleep;
        % t = [t;t(end)+3];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Patch NaNs
        if sum(isnan(x))>0
            x = medfilt_nan(x,4);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EWS analysis
        [ews_out,tau] = ews_sleep_paper(x,t,ews_win,1);
        [pval] = ews_pval_paper(x,t,ews_win,1,ews_surr,tau);

        ews_out_all{sbj}.(ft_name) = ews_out;

        fdnames = fieldnames(ews_out);
        fdnames = strcat(fdnames,['_',ft_name]);
        varname_pval = [varname_pval;fdnames];
        varname_tau = [varname_tau;fdnames];
        val_pval = [val_pval,transpose(cell2mat(struct2cell(pval)))];
        val_tau = [val_tau,transpose(cell2mat(struct2cell(tau)))];

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Table construction
        if isempty(T_pval)
            T_pval = array2table(val_pval,'VariableNames',varname_pval);
            T_tau = array2table(val_tau,'VariableNames',varname_tau);
            T_bds = array2table(val_bds,'VariableNames',varname_bds);


        else
            T_pval_new = array2table(val_pval,'VariableNames',varname_pval);
            T_tau_new = array2table(val_tau,'VariableNames',varname_tau);
            T_bds_new = array2table(val_bds,'VariableNames',varname_bds);

            T_pval = [T_pval;T_pval_new];
            T_tau = [T_tau;T_tau_new];
            T_bds = [T_bds;T_bds_new];

        end
        disp(['Subject ',num2str(sbj),' finished'])
    end

end
% 
% clear ft_sbj
save('EWS_Individuals.mat')

%% Further processing for plots

ews_to_avg = {'AR1','StD'};
% ft_to_avg = 'SP_Summaries_welch_rect_area_5_1';
ft_to_avg = 'DistOnset';

ews_out_all1 =  ews_out_all(~cellfun(@isempty, ews_out_all));

sbjsamp = 1;

max_ews_len = length(ews_out_all{sbjsamp}.(ft_to_avg).(ews_to_avg{1}));           % Maximum EWS time-series length
N_valid = size(T_pval,1);

% Average sleep scores
maxlents = length(ts_sbj_nomed{sbjsamp}.ftdist_noart);

t_sleep = (epcs_to_asleep(sbjsamp)-1)*score_size;
t = ((stp_ftepcs_all_filt{sbjsamp}-1)/Fs) - t_sleep;
% t = [t;t(end)+3];
tews = t(end-max_ews_len+1:end);

ews_all = [];

for news = 1:length(ews_to_avg)

    kk = 0;

    figure
    hold on
    allewsts = NaN(N_valid,max_ews_len);
    avg_score = NaN(N_valid,maxlents);
    cnt = 0;

    for sbj = 1:length(ews_out_all)
        
        if ~isempty(ews_out_all{sbj})
            cnt = cnt+1;
            if ifpostdepthsmall(cnt)
                continue
            end
            if length((ews_out_all{sbj}.(ft_to_avg).(ews_to_avg{news}))) <max_ews_len
                continue
            end

            kk = kk +1;
            ews_tsnow = ews_out_all{sbj}.(ft_to_avg).(ews_to_avg{news});
            len_ewsnow = length(ews_tsnow);

            allewsts(kk,max_ews_len-len_ewsnow+1:end) = ews_tsnow;

            oriscoresbj = oriscore_asleep{sbj};
            idxnow = stp_ftepcs_all_filt{sbj};
            % idxnow = [idxnow;idxnow(end)+768];
            stpepcs = floor(idxnow/(Fs*score_size))+1;
            avg_score(kk,end-maxlents+1:end) = oriscoresbj(stpepcs);

            plot(tews,allewsts(kk,:),'color', [.5 .5 .5 .2])


        end

    end

    ews_all{news} = allewsts;

    plot(tews,mean(allewsts,1,'omitnan'),'k','LineWidth',2)
    xlabel('Time (s)')
    ylabel(ews_to_avg{news})

    box off
    set(gca,'TickDir','out')
    set(gca,'ticklength',2*get(gca,'ticklength'))
    set(gca,'lineWidth',2)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add statistical testings

base_per = 3*60;         % 5 minutes baseline to compare
base_epcs = floor((base_per-epc_len)/jmp)+1;

pvec_ews = [];
for news = 1:length(ews_to_avg)
    curewsmat = ews_all{news};
    pvec_ews_now = tvtest(curewsmat,base_epcs,'left');
    pvec_ews{news} = pvec_ews_now;
end

avg_score_cutews = avg_score(:,end-length(tews)+1:end);
pvec_slp = tvtest(avg_score_cutews,base_epcs,'left');

save('Figure2f.mat',"pvec_slp",'tews','avg_score','ews_all',"pvec_ews",'base_epcs','kk')







































