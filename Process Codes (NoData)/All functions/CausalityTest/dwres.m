function dw = dwres(res)

% Author: Junheng Li
% Durbin-Watson test for serial autocorrelation of Vector-autoregression
% models
%
% The DW test stats output should be close to '2' to be defined as no
% significant serial autocorrelation;
% The MATLAB's function does not currently support test foe VAR model.
%
% DW=n−1∑i=1(ri+1−ri)2n∑i=1r2i,

% Input:
% res: should be a NxM matrix, where N is the number of features
% (variables) and M is the number of observations (Time points)
%
% Output:
% dw: the DW test stats for each variable (Nx1 vector)

nobs = size(res,2);
nv = size(res,1);

dw = zeros(nv,1);

for nnv = 1:nv    % For each variable residual
    nom = 0;
    denom = 0;
    for ii = 1:nobs - 1
        nom = nom + (res(nnv,ii))^2;
        denom = denom + (res(nnv,ii+1) - res(nnv,ii))^2;
    end
    dw(nnv) = denom/nom;

end




