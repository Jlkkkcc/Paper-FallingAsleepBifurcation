%% Order of integration test
% Author: Junheng Li
%
% The general idea is to test how many order differencing can achieve a
% stationary time-series
% Stationarity is measured via KPSS and ADF test (both should be satisfied)

%% 
function ord = ord_int(y)

% Input
% y: the input time-series to test
%
% Output:
% ord: the order of integration I(d)

h_kpss = kpsstest(y);
h_adf = adftest(y);

if h_kpss == 0 && h_adf == 1
    ord = 0;
else
    ord = 0;
    while h_kpss == 1 || h_adf ==0
        ord = ord +1;
        ts_diff_ord = diff(y,ord);
        h_kpss = kpsstest(ts_diff_ord);
        h_adf = adftest(ts_diff_ord);

    end
end





