function theta = circ_unifrnd( varargin )
%CIRC_UNIFRND random sample from circular uniform distribution
%
%  theta=CIRC_UNIFRND( [M N ...] ) returns a matrix of size [M N ...]
%  with randomly sampled angles from the circular uniform distribution.
%
%  theta=CIRC_UNIFRND( M, N, ...) same as previous syntax.
%

%  Copyright 2005-2008 Fabian Kloosterman

if nargin==1 && isscalar(varargin{1})
    n = {varargin{1} 1};
else
    n = varargin;
end

theta = unifrnd(0,2*pi, n{:});
