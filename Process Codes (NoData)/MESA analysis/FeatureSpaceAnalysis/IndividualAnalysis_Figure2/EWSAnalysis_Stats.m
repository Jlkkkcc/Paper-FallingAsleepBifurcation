%% EWS average on all subjects and relevant testings


%%
% Data loading

clear all
clc

load EWS_Individuals.mat
load PostWakeExclu_Idx.mat    % Excluding markers for post-sleep depth

%% Reshape as tables to test statistics of increasing with time

ews_to_avg = 'AR1';
ft_to_avg = 'DistOnset';

sbjsamp = 1;
ews_out_all1 = ews_out_all;
ews_out_all =  ews_out_all(~cellfun(@isempty, ews_out_all));

max_ews_len = length(ews_out_all{sbjsamp}.(ft_to_avg).(ews_to_avg));           % Maximum EWS time-series length
kk = 0;
N_valid = size(T_pval,1);
allewsts = NaN(N_valid,max_ews_len);

t_sleep = epcs_to_asleep(sbjsamp)*score_size;
t = ((stp_ftepcs_all_filt{sbjsamp}-1)/Fs) - t_sleep;
t = [t;t(end)+3];

allews = [];
allsbjid = [];
alltimeidx = [];

cnt = 0;

max_ews = length(ews_out_all{sbjsamp}.(ft_to_avg).(ews_to_avg));

for sbj = 1:length(ews_out_all) 

    if ifpostdepthsmall(sbj)  % Exclude participants with too small POST DEPTH
        continue
    end

    if ~isempty(ews_out_all{sbj})

        if length((ews_out_all{sbj}.(ft_to_avg).(ews_to_avg))) <max_ews
            continue
        end
        
        kk = kk +1;
        ews_tsnow = ews_out_all{sbj}.(ft_to_avg).(ews_to_avg);
        len_ewsnow = length(ews_tsnow);
        
        allewsts(kk,max_ews_len-len_ewsnow+1:end) = ews_tsnow;

        for jjj = 1:len_ewsnow
            cnt = cnt+1;
            allews(cnt) = ews_tsnow(jjj);
            allsbjid(cnt) = sbj;
            alltimeidx(cnt) = t(jjj+100-1)/60;

        end

    end

end

alltimeidx = alltimeidx - min(alltimeidx);  % Turning into positive

Tfinal = table((allews'),allsbjid',alltimeidx','VariableNames',{'EWS','Sbjidx','Time'});


%% GLME: EWS index increasing with time testing

[~,TF] = rmoutliers(Tfinal.EWS);
Tclean = Tfinal(~TF,:);

glme_clean = fitlme(Tclean,'EWS ~ 1 + Time + (Time|Sbjidx)');

%% Repeat again with Variance

ews_to_avg = 'StD';
ft_to_avg = 'DistOnset';

sbjsamp = 1;
ews_out_all1 = ews_out_all;
ews_out_all =  ews_out_all(~cellfun(@isempty, ews_out_all));

max_ews_len = length(ews_out_all{sbjsamp}.(ft_to_avg).(ews_to_avg));           % Maximum EWS time-series length
kk = 0;
N_valid = size(T_pval,1);
allewsts = NaN(N_valid,max_ews_len);

t_sleep = epcs_to_asleep(sbjsamp)*score_size;
t = ((stp_ftepcs_all_filt{sbjsamp}-1)/Fs) - t_sleep;
t = [t;t(end)+3];

allews = [];
allsbjid = [];
alltimeidx = [];

cnt = 0;

max_ews = length(ews_out_all{sbjsamp}.(ft_to_avg).(ews_to_avg));

for sbj = 1:length(ews_out_all) 

    if ifpostdepthsmall(sbj)  % Exclude participants with too small POST DEPTH
        continue
    end

    if ~isempty(ews_out_all{sbj})

        if length((ews_out_all{sbj}.(ft_to_avg).(ews_to_avg))) <max_ews
            continue
        end
        
        kk = kk +1;
        ews_tsnow = ews_out_all{sbj}.(ft_to_avg).(ews_to_avg);
        len_ewsnow = length(ews_tsnow);
        
        allewsts(kk,max_ews_len-len_ewsnow+1:end) = ews_tsnow;

        for jjj = 1:len_ewsnow
            cnt = cnt+1;
            allews(cnt) = ews_tsnow(jjj);
            allsbjid(cnt) = sbj;
            alltimeidx(cnt) = t(jjj+100-1)/60;

        end

    end

end

alltimeidx = alltimeidx - min(alltimeidx);  % Turning into positive

Tfinal = table((allews'),allsbjid',alltimeidx','VariableNames',{'EWS','Sbjidx','Time'});

[~,TF] = rmoutliers(Tfinal.EWS);
Tclean = Tfinal(~TF,:);

glme_clean = fitlme(Tclean,'EWS ~ 1 + Time + (Time|Sbjidx)');



