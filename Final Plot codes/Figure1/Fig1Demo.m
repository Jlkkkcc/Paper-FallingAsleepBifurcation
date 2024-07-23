
%% Figure1 Demo EEG plots

% Authors: Junheng Li / Tianyu Wei

% Exemplary EEG data / epochs during falling asleep
%
% Raw data:
% eeg: Raw EEG data chunks of a falling asleep process
% Fs: the sampling rate;
% scores: Corresponding sleep scores in 30 s size epochs

%% 

clear all
clc
load DataFig1.mat

%% Epochs for specific sleep stages

% Epoch length in seconds
epc_len = 6;

% Entire EEG
S_total = eeg;

% Select some representative epoch for each sleep stage period
i1st = 60*Fs;
S_1 = data(i1st+1:i1st+epc_len*Fs);

i2st = 330*Fs;
S_2 = data(i2st+1:i2st+epc_len*Fs);

i3st = 600*Fs;
S_3 = data(i3st+1:i3st+epc_len*Fs);


%% plot with fixed ylim

% Three exemplary epochs plots
figure
plot(S_1);
ylim([-0.06,0.06])
axis off
figure
plot(S_2);
ylim([-0.06,0.06])
axis off
figure
plot(S_3);
ylim([-0.06,0.06])
axis off

%% Total EEG period with time markers

% Total with time markers of the epochs taken
figure
plot(S_total);
ylim([-0.06,0.06])
axis off
hold on
yl = ylim;
line([i1st,i1st],ylim,'LineStyle','--', 'Color','red')
line([i1st+epc_len*Fs,i1st+epc_len*Fs],ylim,'LineStyle','--', 'Color','green')
line([i2st,i2st],ylim,'LineStyle','--','Color','red')
line([i2st+epc_len*Fs,i2st+epc_len*Fs],ylim,'LineStyle','--', 'Color','green')
line([i3st,i3st],ylim,'LineStyle','--','Color','red')
line([i3st+epc_len*Fs,i3st+epc_len*Fs],ylim,'LineStyle','--', 'Color','green')


