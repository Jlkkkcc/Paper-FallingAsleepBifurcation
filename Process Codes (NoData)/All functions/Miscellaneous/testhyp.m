function [n1n2_noFA,non1n2_noFA] = testhyp(hyp)

% Quick test of hypnogram for N1N2 transitions and relevant FA

% JL

% Test continuous stage 0s

hyp_zero = (hyp==0);
[~, ~, k1] = ZeroOnesCount(hyp_zero);

non1n2_noFA = k1;

hyp_stages = (hyp>0);
[start, len, k2] = ZeroOnesCount(hyp_stages);
k2new = k2;

hyp_nan = isnan(hyp);
hyp_nan_lefts = [hyp_nan(2:end),0];
hyp_nan_rights = [0,hyp_nan(1:end-1)];

% Detect overlap modified
ovlaps_total = 0;
for ii = 1:k2new
    idxsta = max(1,start(ii)-1);
    idxstop = min(start(ii)+len(ii)-1,length(hyp_nan));
    if (hyp_nan(idxsta)==1)  || (hyp_nan(idxstop) == 1)
        ovlaps_total = ovlaps_total+1;
    end

end
% 
% 
% ovlaps_total = sum(hyp_nan_lefts.*hyp_stages) + sum(hyp_nan_rights.*hyp_stages);

n1n2_noFA = k2 - ovlaps_total;
% if n1n2_noFA<0
%     n1n2_noFA = 0;
% end





















