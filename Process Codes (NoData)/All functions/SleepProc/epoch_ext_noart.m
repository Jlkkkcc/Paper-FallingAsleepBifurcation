function [epcs, ifgood] = epoch_ext_noart(data, Fs, epc_len, num_epc,overlap)

% This function is used to extract artefact-free epochs from a period of
% data input
%
% Author: Junheng Li

% Input:
% data: the data vector for epoch extraction
% Fs: sampling rate of the input data
% epc_len: epoch length in data points
% num_epc: Number of epochs for output
% overlap: percentage of overlap between continuous epochs (control for artefact
% detection)
% ifplot: whether to plot the epochs chosen along with artefact detection
%
% Output:
% epcs: A N-by-M matrix, where N is number of epochs, and M is length in
% data points

len = length(data);
ifgood = 1;
epcs = [];

if epc_len*num_epc > len
    error('Too short input data!')
end

% Maximum possible number of epochs from input data
jump = floor((1-overlap)*epc_len);
num_epcs_max = floor((len-epc_len)/jump)+1;
epcs_all = zeros(num_epcs_max,epc_len);

% Divide into epochs
for i = 1:num_epcs_max
    stp = 1+(i-1)*jump;
    epcs_all(i,:) = data(stp:stp+epc_len-1);
end

% % Epoch-level stats
% epc_rms = rms(epcs_all,2);    % Compute the root mean square for each epoch
% 
% amp_th = 3*nanmean(epc_rms);      % Amplitude threshold defined as 3 times the mean RMS
% 
% isart = zeros(num_epcs_max,1);    % Logic indicater of artefacts
% 
% % Artefact detection
% for i = 1:num_epcs_max
%     if max(abs(epcs_all(i,:))) > amp_th
%         isart(i) = 1;
%         if ifplot
%             stp = 1+(i-1)*jump;
%             plot(t(stp:stp+epc_len-1),epcs_all(i,:),'r')
%         end
%     end
% end
% 
% % Choose artefact-free epochs
% 
% idx_art_free = find(isart == 0);
% if isempty(idx_art_free)
%     ifgood = 0;
%     return
%     
% end
% if length(idx_art_free)<num_epc
%     ifgood = 0;
%     return
% end

idx_rand = randperm(num_epcs_max,num_epc);
% idx_to_use = idx_art_free(idx_rand);

epcs = epcs_all(idx_rand,:);

