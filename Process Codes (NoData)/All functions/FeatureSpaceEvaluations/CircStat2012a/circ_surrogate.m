function pval = circ_surrogate(ph)

% Surrogate data testing for phase vectors

ph_sur = circ_unifrnd(length(ph));

pval = circ_cmtest(ph,ph_sur);