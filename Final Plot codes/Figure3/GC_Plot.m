%% Granger-Causality test final results (Figure 3G)

% See analysis in the process codes

clear all
clc
load Figure3G.mat

G = digraph(adjmat_GCtest_corBon,ft_names);
p = plot(G,'MarkerSize',10,'LineWidth',1.5)







