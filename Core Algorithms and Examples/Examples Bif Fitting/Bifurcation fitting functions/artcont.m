function [idx_st, num_cont] = artcont(x)

% This is a simple function identifying the continous artefact points
% (starting indexes) and continous lengths;
%
% The input x should be a logical input, indicating whether each epoch is
% artefact (1) or not (0);
%
% Author: Junheng Li

f = find(diff([0,x,0]==1));
idx_st = f(1:2:end-1);  % Start indices
num_cont = f(2:2:end)-idx_st;  % Consecutive ones’ counts