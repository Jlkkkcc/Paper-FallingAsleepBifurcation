%% Median filtering to replace NaN in sleep feature time-series
%
% Author: Junheng Li
%
% The main aim of this function is to replace the NaN values in the feature
% time-series (where NaN indicate artefact epochs) using a median filtering
% with a certain moving-median window
%
%% Code

function yout = medfilt_nan(y,win)
%
% Input:
% y: The time-series for median filtering
% win: Window size for moving-median
% In general , it is assumed that the maximum number of continuous NaNs in
% y should be less than 'win', otherwise error will be reported
%
% Output:
% yout: the output time-series with NaN padded with median of neighbouring
% window

len = length(y);
if len<= win
    error('Window size larger than time series length')
end

yout = y;
half_win = floor(win/2);
% For datapoints before
for i = 1:half_win
    curr_win = y(1:win);
    if sum(1-isnan(curr_win)) == 0
        error('Continuous NaNs greater than window size')
    end
    if isnan(y(i))
        yout(i) = nanmedian(curr_win);
    end
end

for i = half_win+1:len-half_win
    if isnan(y(i))
        curr_win = y(i-half_win+1:i+half_win);
        if sum(1-isnan(curr_win)) == 0
            error('Continuous NaNs greater than window size')
        end
        yout(i) = nanmedian(curr_win);
        
    end
end

for i = len-half_win+1:len
    if isnan(y(i))
        curr_win = y(len-win+1:len);
        if sum(1-isnan(curr_win)) == 0
            error('Continuous NaNs greater than window size')
        end
        yout(i) = nanmedian(curr_win);
        
    end
end



