%% Time-varying statistical testings

function pvec = tvtest(datamat, basewin,tail)

% This function is aimed to do statistical testing over time compared to
% the baseline period specified;

% Specifically for EWS and sleep stages we test for right tail
%
% Author: JL

% Input:
% datamat: NxT matrix, N is the number of samples, T is observations over
% time;
% basewin: baseline period defined (as number of data points)
% tail: a string indicating which tail to test

% N = size(datamat,1);

T = size(datamat,2);

if T<= basewin
    error('Baseline period longer than entire time series')
end

base_data = datamat(:,1:basewin);
% base_data = base_data(:);

% Normalise to the second direction
base_data_avg = mean(base_data,2,'omitnan');

pvec = NaN(T-basewin,1);

% Test over windows
for ii = basewin+1:T
    
    curdata = datamat(:,ii);

    [~,pvec(ii-basewin)] = ttest2(base_data_avg,curdata,'Tail',tail);

end

























