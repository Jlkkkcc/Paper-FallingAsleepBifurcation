%% Adding features evaluation

% Junheng Li

%% Data load

clear all
clc

load fteval_EWSallsbj_20Dec.mat

%% Add features to original structure data

epc_len = 6;            % Epoch length for feature evaluation
ovlap = 0.5;              % Overlap between continous epochs
post_sleep_per = 20;    % Post sleep period for feature evaluation

epc_len_data = epc_len*Fs;
ovlap_data = floor(epc_len_data*ovlap);
rm_data = epc_len_data-ovlap_data;

%%%%%%%%%
% New storage

channel_mkr = NaN(num_sbjs,num_chs);
ft_sbj_patch = ft_sbj;   % Directly patch in new features
ft_sbj_new = [];     % Store only new features
num_newft = 8;
num_ftori = 47;

save_num = 100;
cnt_save = 2;

%%
for i = 1:num_sbjs

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Subject removal previously
    if ismember(i,discard_idx)
        continue
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluate back to channel removals
    channel_mkr_this = [];
    eeg_this = eeg_asleepper_all{i};
    numch_avail = length(eeg_this.ch);
    if numch_avail == 1
        channel_mkr_this = [1,0,0];
    elseif numch_avail == 2
        if isempty(eeg_this.ch{1})
            channel_mkr_this = [0,1,0];
        else
            channel_mkr_this = [1,1,0];
        end
    elseif numch_avail == 3
        channel_mkr_this = [1-isempty(eeg_this.ch{1}),1-isempty(eeg_this.ch{2}),1];
    end
    channel_mkr(i,:) = channel_mkr_this;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % New feature evaluation per channel

    num_epcs_total = num_ftepcs_all(i);
    stp_ftepcs = stp_ftepcs_all{i};
    for jj = 1:num_epcs_total

        % Repair and patch matrices firstly
        ft_chraw = ft_sbj_patch{i}.ck{jj}.ft_ch;
        ft_chpatch = NaN(num_chs,num_ftori+num_newft);
        ft_chnew = NaN(num_chs,num_newft);
        % Loop channel
        stp = stp_ftepcs(jj);

        for ch = 1:length(channels)
            if ~channel_mkr_this(ch)   % Check if this channel is deleted
                % Check
                if size(ft_chraw,1)>=ch && sum(ft_chraw(ch,:)~=0)>0
                    error('Check')
                end
                continue
            end
            eeg_now = eeg_this.ch{ch};
            ft_chpatch(ch,1:num_ftori) = ft_chraw(ch,:);
            data_unit = eeg_now(stp:stp+epc_len_data-1);

            if max(data_unit) == 0.5 && i>=82
                if sum(~isnan(ft_chpatch(ch,:)))>0
                    error('check 2')
                end
            else
                [ftt,ft_dscrp_new] = ftadd(data_unit,Fs);
                
                ftt1 = struct2cell(ftt);
                ftt1 = transpose(cell2mat(ftt1));
                ft_chpatch(ch,num_ftori+1:end) = ftt1;
                ft_chnew(ch,:) = ftt1;
            end
           
        end

        ft_sbj_patch{i}.ck{jj}.ft_ch = ft_chpatch;
        ft_sbj_new{i}.ck{jj}.ft_ch = ft_chnew;


    end


    % 
    % for ch = 1:length(channels)
    % 
    %     if ~channel_mkr_this(ch)   % Check if this channel is deleted
    %         continue
    %     end
    % 
    %     eeg_now = eeg_this.ch{ch};
    %     num_epcs_total = num_ftepcs_all(i);
    %     stp_ftepcs = stp_ftepcs_all{i};
    % 
    %     for jj = 1:num_epcs_total
    %         stp = stp_ftepcs(jj);
    %         data_unit = eeg_ck(stp:stp+epc_len_data-1);     % Current data chunk
    %         if max(data_unit) == 0.5
    %             ft_chnow = ft_sbj_patch{i}.ck{jj}.ft_ch(ch,:);
    % 
    %             ft_sbj_patch{i}.ck{jj}.ft_ch(ch,end+1:end+num_newft) = NaN(1,num_newft);
    %             ft_sbj_new{i}.ck{jj}.ft_ch(ch,:) = NaN(1,num_newft);
    %         else
    %             [ftt,ft_dscrp_new] = ftadd(data_unit,Fs);
    % 
    %             ftt1 = struct2cell(ftt);
    % 
    %             ft_sbj{i}.ck{jj}.ft_ch(ch,:) = transpose(cell2mat(ftt1));
    %         end
    % 
    %     end
    % 
    % end
    % stp_ftepcs_all{i} = stp_ftepcs;
    
    disp(['Subject No.',num2str(i),' successfully patched'])
    
    if floor(i/save_num)>cnt_save
        save(['fteval_allsbj_temp',num2str(i),'_light.mat'],'ft_sbj_new',"ft_dscrp_new","ft_sbj_patch")
        cnt_save = cnt_save +1;
    end
end


save('fteval_allsbj_26Dec24_all.mat','-v7.3')
clear eeg_asleepper_all eeg_this
save('fteval_allsbj_26Dec24_light.mat','-v7.3')














