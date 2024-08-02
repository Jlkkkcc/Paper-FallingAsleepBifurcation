%% EWS plot function
%
% This function visualises EWS measurements along with its time-series and
% p-values, and Kendall's tau if relevant
%
% Author: Junheng Li

%% Code

function ews_plot_paper(ts,time,win, ews,tau,pval,tchange)

% Inputs:
% ts: Original time series signal
% time: The corresponding time indexes of 'ts'
% win: Window size for EWS calculation
% tau: Kendall's tau for EWS evaluation
% pval: pvalues for Kendall's tau illustrating significance of EWS
%
% Output:
% Figure plot of original time series and EWS measurements, all time
% aligned;

f = figure;
fdnames = fieldnames(ews);
num_ews = length(fdnames);

subplotnum = num_ews+1;
numcol = ceil(subplotnum/2)+1;
f.Position(3:4) = [700,numcol*220];

time_base = 60;

subplot(numcol,2,1)
plot(time,ts)
xlim([floor(min(time)/time_base)*time_base,ceil(max(time)/time_base)*time_base]);
xticks(min(xlim):time_base*5:max(xlim))
line([tchange,tchange],[min(ylim),max(ylim)],'Color','r')
% line([N1_time,N1_time],[min(ylim),max(ylim)],'Color','g')
xlabel('Time (s)')
box off
set(gca,'TickDir','out')
    set(gca,'ticklength',2*get(gca,'ticklength'))
title('Original time series')

for i = 1:subplotnum-1
    
    ews_tsnow = ews.(fdnames{i});
    tau_now = tau.(fdnames{i});
    pval_now = pval.(fdnames{i});
    
    if size(ews_tsnow)>1
        ews_tsnow = transpose(ews_tsnow);
    end
    ews_tsnow = [NaN(1,win-1),ews_tsnow];            % Pad the length due to window effects
    
    ax = subplot(numcol,2,i+1);
    plot(time, ews_tsnow)
    box off
%     set(gca,'FontSize', 16)
    set(gca,'TickDir','out')
    set(gca,'ticklength',2*get(gca,'ticklength'))
%     if contains(fdnames{i},'AR')
%         ylim([0,1.2])
%     end

%     set(gca,'lineWidth',2)
%     xlabel('Time (s)')
    line([tchange,tchange],[min(ylim),max(ylim)],'Color','r')
    xlim([floor(min(time)/time_base)*time_base,ceil(max(time)/time_base)*time_base]);
    xticks(min(xlim):time_base*5:max(xlim))
%     line([N1_time,N1_time],[min(ylim),max(ylim)],'Color','g')
    title(fdnames{i})
    txt = ['Tau=',num2str(tau_now,2), ';p-value=',num2str(pval_now,2)];
%     xlim = ax.XLim;
%     ylim = ax.YLim;
%     text(0.02*diff(xlim),0.90*ylim(2),txt)
    text(min(xlim), max(ylim)*0.95,txt)

end


