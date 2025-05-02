%% Demo MESA checks

% Author: JL

%% 

clear all

clc

load MESA_DemoInfo.mat

%% Loading all subject info

% The subjects are now exclused based on three criteria: bad quality, all
% channels discarded, and did not start from awake baseline

num_parts = 52;
n_part = 40;

discard_info_all = [];
time_to_sleep_all = [];
time_to_sleep_all_correct = [];
diff_list_all = [];

sleep_info_correct_all = [];
epcs_exc_1 = [];

for i = 1:num_parts
    load(['Subject_info_part_',num2str(i),'.mat'])
    
    discard_info_all = [discard_info_all;discard_info];
    time_to_sleep_all = [time_to_sleep_all;time_to_sleep];

    % Additional info
    diff_list_all = [diff_list_all,diff_list + (i-1)*n_part];
    time_to_sleep_all_correct = [time_to_sleep_all_correct;time_to_sleep_correct];

    epcs_exc_1 = [epcs_exc_1;sleep_info(:,1)];
    
    sleep_info_correct_all = [sleep_info_correct_all;sleep_info_oriscore_nonan];

end
discard_info_sum = sum(discard_info_all,2);
  % This saves the length of the initial period cut (initial artefactual periods for each participant)
save('epcsexc.mat','epcs_exc_1')  

%% Sleep onset time distribution and divisions

idx_disc = (discard_info_sum == 3);
% idx_disc = (time_to_sleep_all == 0);
time_to_sleep_all_correct(idx_disc) = NaN;
time_to_sleep_all_correct(time_to_sleep_all_correct<3*60) = NaN;


time_to_sleep_per = prctile(time_to_sleep_all_correct,[20,80]);
time_to_sleep_per = time_to_sleep_per/60;

Age(idx_disc) = NaN;
Gender(idx_disc) = NaN;

figure
histogram(time_to_sleep_all_correct/60,[3:3:90])
xlabel('Time to sleep (min)')
ylabel('Number of participants')
line([time_to_sleep_per(1),time_to_sleep_per(1)],[0,60],'Color','r','LineWidth',2,'LineStyle','--')
line([time_to_sleep_per(2),time_to_sleep_per(2)],[0,60],'Color','r','LineWidth',2,'LineStyle','--')
box off
xticks([0:10:90])

set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

%% Age final restricted sleep-onset latency

idx_disc = (discard_info_sum == 3);
time_to_sleep_all_correct(idx_disc) = NaN;
time_to_sleep_all_correct(time_to_sleep_all_correct<3*60) = NaN;

idxdisc_all = isnan(time_to_sleep_all_correct);

Age(idxdisc_all) = NaN;
Gender(idxdisc_all) = NaN;

bin_edges = [54:1:95];

figure
histogram(Age(Gender==0),bin_edges)
xlabel('Age (years)')
ylabel('Number of subjects')
title('Female')
ylim([0,40])
set(gca,'TickDir','out');
set(gca,'linewidth',2)
set(gca,'FontSize', 12)
box off
xticks(55:5:95)

figure
histogram(Age(Gender==1),bin_edges)
xlabel('Age (years)')
ylabel('Number of subjects')
title('Male')
set(gca,'TickDir','out');
set(gca,'linewidth',2)
set(gca,'FontSize', 12)
box off
xticks(55:5:95)

ylim([0,40])

ages_female = Age(Gender==0);
ages_male = Age(Gender==1);
disp('Males:')
nanmean(ages_male)
nanstd(ages_male)
disp('Female:')
nanmean(ages_female)
nanstd(ages_female)
disp('Test:')
grp = [repmat({'Males'},[length(ages_male),1]);repmat({'Females'},[length(ages_female),1])];
p1 = kruskalwallis([ages_male;ages_female],grp);


%% Included participants

sbj_remained = find((1-idxdisc_all)==1);
save("Sbjs_filtered_time.mat","sbj_remained")






