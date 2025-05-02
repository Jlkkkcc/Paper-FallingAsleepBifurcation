%% Supplementary Figure 5

% JL

clear all
clc
load SuppFig5.mat

%% Plots all
% all subjects - SuppFig5a
sbjuniqs = unique(night_patidx);
numsbj = length(sbjuniqs);

mkrs = {"diamon","o","^","pentagram","v","square"}.';
maxlen1 = length(mkrs);
colorsall = distinguishable_colors(10);
colorsall = colorsall([2,1,4,5,9,10],:);

figure
hold on
clcidx = 1;
for ii = 1:numsbj

    idxthissbj = (night_patidx == sbjuniqs(ii));
    mkridx = floor((ii-1)/maxlen1)+1;
    % if mkridx == 2 || mkridx==6 ||  mkridx == 3
        scatter(Y2all(idxthissbj,1),Y2all(idxthissbj,2),60,colorsall(clcidx,:),"filled",mkrs{mkridx})
    % else
    %     scatter(Y2all(idxthissbj,1),Y2all(idxthissbj,2),60,colorsall(clcidx,:),mkrs{mkridx})
    % end
    clcidx = clcidx+1;
    if clcidx>maxlen1
        clcidx = 1;
    end

end
box off
set(gca, 'FontSize', 12);
set(gca, 'TickDir', 'out');
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

xlabel('TSNE-1')
ylabel('TSNE-2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consistent subject - SuppFig5b
sbjtest  =11;
night_thissbj = find(ismember(night_patidx,sbjtest));
Y2thissbj = Y2all(night_thissbj,:);

figure
scatter(Y2all(:,1),Y2all(:,2),'blue','filled')
hold on
scatter(Y2all(night_thissbj,1),Y2all(night_thissbj,2),'red','filled')
box off
set(gca, 'FontSize', 12);
set(gca, 'TickDir', 'out');
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

xlabel('TSNE-1')
ylabel('TSNE-2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inconsistent subject - SuppFig5c
sbjtest  = 9;
night_thissbj = find(ismember(night_patidx,sbjtest));
Y2thissbj = Y2all(night_thissbj,:);

figure
scatter(Y2all(:,1),Y2all(:,2),'blue','filled')
hold on
scatter(Y2all(night_thissbj,1),Y2all(night_thissbj,2),'red','filled')
box off
set(gca, 'FontSize', 12);
set(gca, 'TickDir', 'out');
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

xlabel('TSNE-1')
ylabel('TSNE-2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Partly consistent subject - SuppFig5d
sbjtest  = 30;
night_thissbj = find(ismember(night_patidx,sbjtest));
Y2thissbj = Y2all(night_thissbj,:);

figure
scatter(Y2all(:,1),Y2all(:,2),'blue','filled')
hold on
scatter(Y2all(night_thissbj,1),Y2all(night_thissbj,2),'red','filled')
box off
set(gca, 'FontSize', 12);
set(gca, 'TickDir', 'out');
set(gca,'ticklength',2*get(gca,'ticklength'))
set(gca,'lineWidth',2)

xlabel('TSNE-1')
ylabel('TSNE-2')













