%% EWS surrogate testing
%
% Author: Junheng Li
%
% This script is a follow-up step following the ews_sleep function, to
% generate surrogate time-series data for statistical testing;
% 
% Methods for surrogate generation:
%  (see Ebisuzaki W. (1997) ;
% A method to estimate the statistical significance of a correlation when 
% the data are serially correlated. J Climate, 10, 2147-2153).
% The p-value for TVAR is built-in with the model;

%% Code script

function [pval] = ews_pval_paper(ts,time, win, ifdetrend, nsurr ,tau_ori)

% Inputs:
% ts: Original time series;
% time: Relevant time index for 'ts';
% nsurr: Number of surrogates for the statistical testing
% tau_ori: Original Kendall's tau correlation values of the EWS measured
% from real signal;
% Cntr: which version of function to use for EWS
%
% Outputs:
% pval: A structure containing same fields of 'tau_ori', returning a
% surrogate-tested p-value for each of the EWS measurements

fields = fieldnames(tau_ori);
numfield = length(fields);

ts_surr = ebisuzaki(ts,nsurr);              % Surrogate data generating
tau_surr = cell(nsurr,1);


% Surrogate EWS measurements
parfor i = 1:nsurr
    
    [~,tau_surr{i}] = ews_sleep_paper(ts_surr(:,i),time,win, ifdetrend);
%     disp('---------------------------------------')
%     disp(['Surrogate data EWS measurements: NO.', num2str(i),' successful'])

end

for nfd = 1:numfield
    fdname = fields{nfd};
    
    tau_thisews_ori = tau_ori.(fdname);
    tau_surr_thisews = zeros(nsurr,1);
    for j = 1:nsurr
        tau_surr_thisews(j) = tau_surr{j}.(fdname);
    end
    
    pval.(fdname) = (sum(tau_surr_thisews(:)>tau_thisews_ori)+1)/nsurr;      % Surrogate p-values
    
end

end






