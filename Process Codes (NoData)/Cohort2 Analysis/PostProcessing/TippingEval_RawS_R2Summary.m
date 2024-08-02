
%% Tipping point evaluations - Raw S dynamics

% Junheng Li

%% 

clear all
clc
load FittingStatsTable_AllParticipants_New.mat

% Clean up exclusions
tblnew = fittingStatsTable;
tblnew([252:259,279,281,283,284],:) = [];   % Exclusion

fittingStatsTable = tblnew;

%% Calculations

numnight = size(fittingStatsTable,1);

t_critc = NaN(numnight,1);
c_critc = NaN(numnight,1);
ifNoC = NaN(numnight,1);
iftransC = NaN(numnight,1);
ifmiss = fittingStatsTable.ifMissingData;

for ii = 1:numnight

    if ifmiss(ii)  % Empty one
        continue
    end

    tvec_this = fittingStatsTable.UpdatedTimeVec{ii};
    c_3sol = fittingStatsTable.C_3_Solutions{ii};
    c = fittingStatsTable.Control_Parameter{ii};

    if ~isempty(c_3sol) && (length(c_3sol)<length(c) )  % No solution or full solution

        idx_critic = find(c == c_3sol(end));   % Critical point
        t_critc(ii) = tvec_this(idx_critic);
        c_critc(ii) = c(idx_critic);
        disp(['Case ',num2str(ii),' Critical time ', num2str(t_critc(ii)),' min'])
    else
        if length(c_3sol) == length(c)      % Failure situation
            disp(['Case ',num2str(ii),' Bad bif diagram'])
        end
        % Run transcritical detections
        sfitnow = fittingStatsTable.Fitted_Bifurcation{ii};
        sfitnow = sfitnow(21:end);
        diff_sfitnow = diff(sfitnow);
        min_deriv_idx = find(diff_sfitnow == min(diff_sfitnow));
        idx_critic = min_deriv_idx+1+20;

        % No Critical point at all
        if abs(min(diff_sfitnow))<0.005 || (idx_critic >= length(sfitnow))
            ifNoC(ii) = 1;
            disp(['Case ',num2str(ii),' No Critical point detected at all'])
        else
            ifTransC(ii) = 1;    % Transcritical case
            t_critc(ii) = tvec_this(idx_critic);
            c_critc(ii) = c(idx_critic);
            disp(['Case ',num2str(ii),' Transcritical, taking minimum time ', num2str(t_critc(ii)),'min'])
        end
    end


end

fittingStatsTable_new = fittingStatsTable;
fittingStatsTable_new = addvars(fittingStatsTable_new,t_critc,ifTransC','NewVariableNames',{'Model_CrticT','ifTransC'});

save("FittingStatsTable_AllParticipants_Cleaned_WithTC.mat",'fittingStatsTable_new')

%% R2 summary with only negative Critical tipping point as successful modelling

r2all = fittingStatsTable_new.R_square;
latall = fittingStatsTable_new.Time2Sleep;

TC_Mdl = fittingStatsTable_new.Model_CrticT;

mask_lat = (latall>=10 & latall<=90);
sum(TC_Mdl(mask_lat)>=0)

mask_all = (latall>=10 & latall<=90) & (TC_Mdl<0);

figure
edgesbin = [0:0.05:1];
pd = fitdist(r2all(mask_all),'Kernel');
pdfbin = 0:0.01:1;
ypdf = pdf(pd,pdfbin);

histogram(r2all(mask_all),edgesbin)
% histfit(tcrtc_new,30,'kernel')
hold on
plot(pdfbin,ypdf*(0.05*length(r2all(mask_all))),'LineWidth',2,'Color','r')
box off
set(gca,'FontSize', 12)
set(gca,'TickDir','out')
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)
xlabel('R2')
ylabel('Count')
xlim([0,1])
r2med = median(r2all(mask_all),'omitmissing');
line([r2med,r2med],ylim,'LineStyle','--','LineWidth',2,'Color','r');

length(r2all(mask_all))
mean(r2all(mask_all))
std(r2all(mask_all))

length(TC_Mdl(mask_all))
mean(TC_Mdl(mask_all))
std(TC_Mdl(mask_all))


