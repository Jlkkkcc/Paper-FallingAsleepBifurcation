
%% Adding post-hoc results towards the prediction result table

% Author: Junheng Li

%% Data loading

clear all
clc
load 'uniqueTableAddedParticipants1.mat'
load FittingStatsTable_AllParticipants_New_WithTC.mat

%% Table repair % This part correctes the table due to some empty-cell issues / Only run once and then comment out this section
% 
% empstr = struct('downwardCrossingTime', [], 'upwardCrossingTime', [], 'timeSpentBelow', [], 'isRealCrossing', [], 'isFinalPrediction', [], 'ifpositive', [], 'sleepStageDownwardCrossing', [], 'sleepStageUpwardCrossing', [], 'sleepStagesDuring',[], 'ifN1N2', []);
% 
% kk = uniqueTable.CrossingDetailsCell{121,1};
% kk = [{kk{1:3}},{empstr},{kk{4:end}}];
% 
% uniqueTable.CrossingDetailsCell{121,1} = kk;

% empstr = struct('downwardCrossingTime', [], 'upwardCrossingTime', [], 'timeSpentBelow', [], 'isRealCrossing', [], 'isFinalPrediction', [], 'ifpositive', [], 'sleepStageDownwardCrossing', [], 'sleepStageUpwardCrossing', [], 'sleepStagesDuring',[], 'ifN1N2', []);
% 
% kk = uniqueTable.CrossingDetailsCell{123,1};
% kk = [{kk{1:4}},{empstr},{kk{5:end}}];
% 
% uniqueTable.CrossingDetailsCell{123,1} = kk;

%% Add into the prediction table the Post-hoc model evaluation details

num_train = size(uniqueTable,1);

idxtrains = uniqueTable.Patient_ID;
nightmodel = uniqueTable.Nights_for_Modelling;
maxnights = uniqueTable.Max_Number_Nights_Available;
% sestim = uniqueTable.Estimated_S_variables;
testnig_idxall = uniqueTable.TestingNights;

rawfitidx = fittingStatsTable_new.Patient_ID;
rawnightidx = fittingStatsTable_new.Night_Number;
rawmdl_tcrt = fittingStatsTable_new.Model_CrticT;
laten_raw = fittingStatsTable_new.Time2Sleep;
ifTransC = fittingStatsTable_new.ifTransC;
ifNoC_Raw = fittingStatsTable_new.ifNoC3;
R2fit_Raw = fittingStatsTable_new.R_square;

% Adding crossing time evaluations
scrtic_Raw = fittingStatsTable_new.SVariableCritical;
tcross_Raw = fittingStatsTable_new.FinalCrossingTime;

tcrtc_alltrains = [];
tcrtc_alltrains_vec = [];

lat_alltrains = [];
lat_alltrains_vec = [];

ifTransC_alltrains = [];
ifTransC_alltrains_vec = [];

ifNoC_alltrains = [];
ifNoC_alltrains_vec = [];

R2raw_alltrains = [];
R2raw_alltrains_vec = [];

SCraw_alltrains = [];
SCraw_alltrains_vec = [];

TCraw_alltrains = [];
TCraw_alltrains_vec = [];

for ii = 1:num_train

    cross_thissbj = uniqueTable.CrossingDetailsCell{ii};
    if isempty(cross_thissbj)
        continue
    end
    
    thissbj = idxtrains{ii};
    % This training subject
    thisnightmdl = str2num(nightmodel{ii});
    testnight_thistrain = testnig_idxall{ii};
    % thistrain_tests = sestim{ii};
    num_test = length(testnight_thistrain);    %Number of testing nights
    
    % The original nights this subject has
    rawfit_idx_thissbj = contains(rawfitidx,thissbj);
    % rawfit_idx_thissbj = find(rawfit_idx_thissbj==1);
    rawfit_nightsidx_thissbj = ismember(rawnightidx,testnight_thistrain);
    idxtestn_Raw = (rawfit_idx_thissbj & rawfit_nightsidx_thissbj);

    cross_thissbj = uniqueTable.CrossingDetailsCell{ii};

    if sum(idxtestn_Raw)>num_test
        error('Check')
    end

    if length(cross_thissbj)~= num_test
        error('Check')
    end
    
    %%%%%%%%%
    % The model critical values for them
    tcrtc_alltrains{ii} = rawmdl_tcrt(idxtestn_Raw);
    tcrtc_alltrains_vec = [tcrtc_alltrains_vec;rawmdl_tcrt(idxtestn_Raw)];

    lat_alltrains{ii} = laten_raw(idxtestn_Raw);
    lat_alltrains_vec = [lat_alltrains_vec;laten_raw(idxtestn_Raw)];

    ifTransC_alltrains{ii} = ifTransC(idxtestn_Raw);
    ifTransC_alltrains_vec = [ifTransC_alltrains_vec;ifTransC(idxtestn_Raw)];

    ifNoC_alltrains{ii} = ifNoC_Raw(idxtestn_Raw);
    ifNoC_alltrains_vec = [ifNoC_alltrains_vec; ifNoC_Raw(idxtestn_Raw)];

    R2raw_alltrains{ii} = R2fit_Raw(idxtestn_Raw);
    R2raw_alltrains_vec = [R2raw_alltrains_vec;R2fit_Raw(idxtestn_Raw)];

    SCraw_alltrains{ii} = scrtic_Raw(idxtestn_Raw);
    SCraw_alltrains_vec = [SCraw_alltrains_vec;scrtic_Raw(idxtestn_Raw)];
    TCraw_alltrains{ii} = tcross_Raw(idxtestn_Raw);
    TCraw_alltrains_vec = [TCraw_alltrains_vec;tcross_Raw(idxtestn_Raw)];

end

uniqueTable_new = uniqueTable;
uniqueTable_new = addvars(uniqueTable_new,tcrtc_alltrains',lat_alltrains',SCraw_alltrains',TCraw_alltrains','NewVariableNames',{'Model_TCrtic','Slp_Latency','SCritic_Raw','TCross_Raw'});

save('uniqueTableAddedParticipants_JLRev.mat','uniqueTable_new',"TCraw_alltrains_vec","SCraw_alltrains_vec","R2raw_alltrains_vec","ifNoC_alltrains_vec","ifTransC_alltrains_vec",'tcrtc_alltrains_vec','lat_alltrains_vec')




