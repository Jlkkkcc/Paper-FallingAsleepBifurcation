%% Summary for prediction outcomes

% Junheng Li

%% Data loading

clear all
clc
% load uniqueTableUpdatedFluctuationIdx.mat
load uniqueTableAddedParticipants_JLRev.mat
uniqueTable = uniqueTable_new;

%% Extraction of the critical threshold crossing details

num_train = size(uniqueTable,1);

max_timeOnset = 30*60;

ftepc_len = 6;
ftepc_jmp = 3;
score_size = 30;

epcs_score_pre_max = floor(max_timeOnset/score_size)-1;
epcs_post_score_now = 20;
max_epcs_slp = epcs_score_pre_max+epcs_post_score_now;
post_gap_now = 1;

% Group level
situations_all = [];    % Situation: 0: Nothing detected; (complete failure)   1: Only positive detected; 2: A valid negative detection
numFA_all = [];
predany_all = [];

FA_n1n2_all = [];
FA_non12_all = [];
mycheck_FA_n1n2_all = [];
mycheck_FA_non12_all = [];

N1N2_NoFA_alln_alltrain = [];
non1n2_noFA_alln_alltrain = [];

hypnofluct_alln_alltrain = [];
FA_rate_alltrain = NaN(num_train,1);
hypnofluct_alltrain = NaN(num_train,1);

ifNotConsis_all = [];

hypnopre_alln_alltrain = [];

tbl_summary = [];
varNames = ["ID","Night_For_Model","MaxNights","Situation","numFA_Allnights","All_Pred"...
    ,"FA_N1N2","FA_NoN1N2","FA_N1N2_My","FA_NoN1N2_My","IfConsis"...
    , "N1N2_NoFA","Wake_NoFA","Original_Cross_Info","FA_Idx","FA_toDiscard","MeanFA_Rate","HypFluct_AllNight","MeanHypFluct"];
ids_all = uniqueTable.Patient_ID;
nightformodel = uniqueTable.Nights_for_Modelling;
maxnights = uniqueTable.Max_Number_Nights_Available;

% Load data for each training
for ii =1:num_train
    
    cross_thisTrain = uniqueTable.CrossingDetailsCell{ii};   

    hypno_thisTrain = uniqueTable.Night_Hypnograms{ii};

    numpred = length(cross_thisTrain);       % Number of prediction nights
    situation_this = NaN(numpred,1);    
    numFA_this = NaN(numpred,1);  
    pred_any_this = NaN(numpred,1);

    FA_withN1N2 = NaN(numpred,1);
    FA_noN1N2 = NaN(numpred,1);
    myowncheckFA_withN1N2 = NaN(numpred,1);
    myowncheckFA_noN1N2 = NaN(numpred,1);
    ifNotConsis = NaN(numpred,1);

    N1N2_NoFA_alln = NaN(numpred,1);
    non1n2_noFA_alln = NaN(numpred,1);

    ifN1N2_mycheck_alln = cell(numpred,1);
    FA_idx_alln = cell(numpred,1);
    exclude_FA_alln = cell(numpred,1);

    fluct_hyp_alln = NaN(numpred,1);
    hypre_alln = NaN(numpred,epcs_score_pre_max);

    for nnp = 1:numpred

        cross_thisn = cross_thisTrain{nnp};
        numcross = length(cross_thisn);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Shape hypnogram
        scrs_this = hypno_thisTrain{nnp};
        if scrs_this(length(scrs_this) - 20) == 2
            scrs_this = scrs_this(1:end - 1);
        end
        % scrs_this = scrs_this(1:end-1);
        tvec_hyp = (1:length(scrs_this)).*score_size;
        tvec_hyp = tvec_hyp-(length(scrs_this)-epcs_post_score_now+1)*score_size;
        tvec_hyp = tvec_hyp/60;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Direct hypnogram metrics
        scrs_preon = scrs_this(1:end-20);
        fluct = double(diff(scrs_preon)>0);
        [start, ~, ~] = ZeroOnesCount(fluct);
        diffscrs = diff(scrs_preon);
        fluct_hyp_alln(nnp) = sum(diffscrs(start));   % Fluctuation 

        hypre_alln(nnp,end-length(scrs_preon)+1:end) = scrs_preon;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check crossings
        if numcross == 1    % Special cases for one crossing
            % Check whether it's empty
            if isempty(cross_thisn.downwardCrossingTime)
                numFA_this(nnp) = 0;
                situation_this(nnp) = 0;
                continue
            end
            % Otherwise there is one detection
            if ~(cross_thisn.isRealCrossing)    % Not valid detection 
                situation_this(nnp) = 0;
                numFA_this(nnp) = 0;
                continue
            end
        end
        
        % Extract all crossings
        downwardtimesthis = [cross_thisn.downwardCrossingTime];
        upwarrdtimes = [cross_thisn.upwardCrossingTime];
        upwarrdtimes(isnan(upwarrdtimes)) = 11;     % No actual upward crossing back, make maximum
        timespentbelow = upwarrdtimes - downwardtimesthis;
        isrealcross = [timespentbelow>=1];
        isfinalpred = [cross_thisn.isFinalPrediction];
        ifpos = [cross_thisn.ifpositive];
        ifN1N2 = [cross_thisn.ifN1N2];

        ifrealPred = sum(isrealcross & [cross_thisn.isFinalPrediction] & ~ifpos);
        numFA_current = sum(~[cross_thisn.isFinalPrediction] & [isrealcross] & ~ifpos );
        numposreal = sum([isrealcross] & [cross_thisn.ifpositive] );    %Number of real post predictions
        ifpospredonly = (1-ifrealPred) * (numposreal>0);

        numFA_this(nnp) = numFA_current;

        if ifrealPred>1
            error('Must be wrong')
        end
        if logical(ifrealPred)    % There is a prediction
            situation_this(nnp) = 2;
            idxreal = find(isfinalpred==1);
            pred_any_this(nnp) = downwardtimesthis(idxreal);   % Real prediction time
        end
        if logical(ifpospredonly)       % Only positive predictions available
            situation_this(nnp) = 1;
            idxallrealpos = ([isrealcross] & [cross_thisn.ifpositive]);
            idxpos = find(idxallrealpos==1,1);   % First one
            pred_any_this(nnp) = downwardtimesthis(idxpos);
        end
        if ~logical(ifrealPred) && ~logical(ifpospredonly)
            situation_this(nnp) = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Deal with false alarms

        % Only search pre-asleep period, and before the real prediction (if
        % exists)

        if ~logical(ifrealPred) 
            idxdown_finalpred = find(tvec_hyp>0,1);
        else
            idxdown_finalpred = find(tvec_hyp<=pred_any_this(nnp),1,'last');
        end

        scrs_pre = scrs_this(1:idxdown_finalpred);
        tvec_pre = tvec_hyp(1:idxdown_finalpred);
        scrs_this_new = scrs_pre;
        % myowncheckFA_N1N2 = NaN;

        clear FA_idx ifN1N2_mycheck exclude_FA 

        if numFA_current>0        % There are false alarms
            
            FA_idx = (isrealcross & ~isfinalpred & ~ifpos);
            FA_idx = find(FA_idx==1);

            if length(FA_idx)~=numFA_current
                error('check')
            end
            
            % Record numbers
            FA_withN1N2(nnp) = sum(ifN1N2(FA_idx),'omitmissing');
            FA_noN1N2(nnp) = numFA_current - FA_withN1N2(nnp);
            
            % Check time for all FA
            downtime_FA = downwardtimesthis(FA_idx);
            uptime_FA = upwarrdtimes(FA_idx);
            ifN1N2_mycheck = NaN(numFA_current,1);
            exclude_FA = NaN(numFA_current,1);

            for nfa = 1:numFA_current

                if uptime_FA(nfa)>0
                    error('Must be wrong')
                end
                if downtime_FA(nfa)>0
                    error('must be wrong')
                end

                idxdown = find(tvec_pre<=downtime_FA(nfa),1,'last');
                idxup = find(tvec_pre>=uptime_FA(nfa),1);
                if isempty(idxup)
                    idxup = length(tvec_pre);
                end
                if logical(ifrealPred)    % Check if too close to real prediction
                    if abs(uptime_FA(nfa)-pred_any_this(nnp))<0.5
                        idxup = idxup+1;
                        exclude_FA(nfa) = 1;   % This FA can be ignored
                    end
                end
                % Find FA period / Check FA overlap on my side

                ifN1N2_mycheck(nfa) = (sum(scrs_pre(idxdown:min(idxup,length(scrs_pre))))>0);
                scrs_this_new(idxdown:idxup-1) = NaN;
                
            end
            myowncheckFA_withN1N2(nnp) = sum(ifN1N2_mycheck,'omitnan');
            myowncheckFA_noN1N2(nnp) = numFA_current - sum(ifN1N2_mycheck,'omitnan');
           
            % Check record
            if sum(ifN1N2(FA_idx)' ~= logical(ifN1N2_mycheck))>0
                warning('Something wrong in FA original')
                ifNotConsis(nnp) = 1;
            end
            ifN1N2_mycheck_alln{nnp} = ifN1N2_mycheck;
            FA_idx_alln{nnp} = FA_idx;
            exclude_FA_alln{nnp} = exclude_FA;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check false negatives: there is a N1N2 transition but no FA
        % prior to asleep
        scrs_this_new = scrs_this_new(1:end-1);      % Exclude the last epoch as t=0
        % scrs_this_new(isnan(scrs_this_new)) = [];
        [n1n2_noFA,non1n2_noFA] = testhyp(scrs_this_new);

        non1n2_noFA_alln(nnp) = non1n2_noFA;
        N1N2_NoFA_alln(nnp) = n1n2_noFA;

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A general FA rate
    FArate = mean(numFA_this,'omitnan');
    hypfluct_meann = mean(fluct_hyp_alln,'omitnan');
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Record all results for all training nights
    if isempty(tbl_summary)
        tbl_summary = table({ids_all{ii}},nightformodel(ii),maxnights(ii),{situation_this}, ...
            {numFA_this},{pred_any_this},{FA_withN1N2},{FA_noN1N2},{myowncheckFA_withN1N2},{myowncheckFA_noN1N2}, ...
            {ifNotConsis},{N1N2_NoFA_alln},{non1n2_noFA_alln},{cross_thisTrain},{FA_idx_alln},{exclude_FA_alln}, ...
            FArate,{fluct_hyp_alln},hypfluct_meann,'VariableNames',varNames);
    else
        tbl_new = table({ids_all{ii}},nightformodel(ii),maxnights(ii),{situation_this}, ...
            {numFA_this},{pred_any_this},{FA_withN1N2},{FA_noN1N2},{myowncheckFA_withN1N2},{myowncheckFA_noN1N2}, ...
            {ifNotConsis},{N1N2_NoFA_alln},{non1n2_noFA_alln},{cross_thisTrain},{FA_idx_alln},{exclude_FA_alln}, ...
            FArate,{fluct_hyp_alln},hypfluct_meann,'VariableNames',varNames);
        tbl_summary = [tbl_summary;tbl_new];
    end
    
    %%%%%%%%%%%%%%%%%%%%
    % Concatenate lists directly
    situations_all = [situations_all;situation_this];
    predany_all = [predany_all;pred_any_this];
    numFA_all = [numFA_all;numFA_this];

    FA_n1n2_all = [FA_n1n2_all;FA_withN1N2];
    FA_non12_all = [FA_non12_all;FA_noN1N2];
    mycheck_FA_n1n2_all = [mycheck_FA_n1n2_all;myowncheckFA_withN1N2];
    mycheck_FA_non12_all = [mycheck_FA_non12_all;myowncheckFA_noN1N2];

    N1N2_NoFA_alln_alltrain = [N1N2_NoFA_alln_alltrain;N1N2_NoFA_alln];
    non1n2_noFA_alln_alltrain = [non1n2_noFA_alln_alltrain;non1n2_noFA_alln];

    hypnofluct_alln_alltrain = [hypnofluct_alln_alltrain;fluct_hyp_alln];
    FA_rate_alltrain(ii) = FArate;
    hypnofluct_alltrain(ii) = hypfluct_meann;

    ifNotConsis_all = [ifNotConsis_all;ifNotConsis];

    hypnopre_alln_alltrain = [hypnopre_alln_alltrain;hypre_alln];

end


%% Check stats of predictions (And failures)

idxmask = (lat_alltrains_vec<10 | lat_alltrains_vec>90);

situations_all(idxmask) = NaN;

good_pred = (situations_all ==2);
mean(predany_all(good_pred))
std(predany_all(good_pred))

bad_pred = (situations_all ==1);
mean(predany_all(bad_pred))
std(predany_all(bad_pred))

disp('Failed entirely')
sum(situations_all==0)

%% Distribution of Post-hoc tipping points (Raw model) (Figure 4E)

% Raw model:
sum(tcrtc_alltrains_vec<0)

idxposmdl = (lat_alltrains_vec>=10 & lat_alltrains_vec<=90) & (tcrtc_alltrains_vec<0);   % Positive model predictions

% After filtering raw failure: post-hoc prediction failures
sum(situations_all(idxposmdl)==1)
sum(situations_all(idxposmdl)==0)

figure
binedges = [-20:0.5:0];
% histogram(predany_all(idxmask),binedges,'Normalization','probability')
% hold on
histogram(tcrtc_alltrains_vec(idxposmdl),binedges,'Normalization','probability')
xlabel('Time (min)')
ylabel('Probability')
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

%% Prediction errors (Figure 4Eii)

idxallinc = (lat_alltrains_vec>=10 & lat_alltrains_vec<=90) & (tcrtc_alltrains_vec<0) & (situations_all==2); 

tdiff = predany_all(idxallinc) - tcrtc_alltrains_vec(idxallinc);

figure
bed = [-15:0.5:15];
histogram(tdiff,bed,'Normalization','probability')
xlabel('Time (min)')
ylabel('Probability')
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
tdmed = median(tdiff,'omitmissing');
line([tdmed,tdmed],ylim,'LineStyle','--','LineWidth',2,'Color','r');

mean(predany_all(idxallinc))
std(predany_all(idxallinc))

mean(tdiff)
std(tdiff)


%% Early Predictinos and relationship to hypnograms (Supplementary Figure)
% 
% numFA_all(isnan(numFA_all)) = 0;
% 
% [r,p] = corr(numFA_all,hypnofluct_alln_alltrain)

hypnofluct_alltrain(isnan(FA_rate_alltrain)) = [];
FA_rate_alltrain(isnan(FA_rate_alltrain)) = [];
[r1,p1] = corr(hypnofluct_alltrain,FA_rate_alltrain)

figure
scatter(hypnofluct_alltrain,FA_rate_alltrain,'filled')
xlim([0.9,4])
ylim([-0.1,3.5])
p = polyfit(hypnofluct_alltrain,FA_rate_alltrain,1);
hold on
% xlim([-0.5,25])
xl = xlim;
dx = diff(xl)/100;
xvec = xl(1):dx:xl(2);
plot(xvec,xvec*p(1)+p(2),'r','LineWidth',2)

box off

xlabel('Fluctuation')
ylabel('False alarm rate')
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)


%% Histogram of errorneous early predictions

figure
histogram(categorical(numFA_all(~idxmask)),'BarWidth',0.6)
box off
xlabel('False alarm')
ylabel('No. of nights')
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

mean(numFA_all(~idxmask),'omitmissing')
std(numFA_all(~idxmask),'omitmissing')



